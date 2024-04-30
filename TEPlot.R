# load required packages
library(tidyverse)
library(plyranges)
library(ggplot2)
library(Biostrings)
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome")

# Load in repeat information from EarlGrey
eg_gff <- read_gff('data/icerya_purchasi.filteredRepeats.divergence.gff') %>%
  mutate(type = as.character(type)) %>%
  mutate(type = sub("-.*", "", type)) %>%
  filter(!is.na(KIMURA80)) %>%
  as_tibble() %>%
  mutate(repeat_number = row_number()) %>%
  as_granges()
# Convert our data into tibble for data manipulation
eg_gff_tbl <- as_tibble(eg_gff)
# Correct the subfamily 'Penelope' (within PLE) as it is miscategorised as LINE by RepeatMasker
# Rename Penelope to move to group PLE 
eg_gff <- eg_gff_tbl %>%
mutate(type = sub('LINE/Penelope','PLE/Penelope', type))
# Get number of repeats characterised by EarlGrey
nrow(eg_gff) # 1078807
# Estimate percentage of repeat elements across entire genome
sum(eg_gff$width) #596663952


# Load in genome sequence of Icerya purchasi (for GC content estimation)
genome_seq <- Biostrings::readDNAStringSet("data/genome/GCA_952773005.1_ihIcePurc2.1_genomic.fna")
# Rename sequence names and remove 'Icerya...'
names(genome_seq) <- sub(' .*','', names(genome_seq))
# Fix up scaffold length for 3 scaffolds (one for each chromosome) 
repeat_tbl <- tibble(seqnames = names(genome_seq), scaffold_length = width(genome_seq))
# Select columns we want to interpret
eg_gff_tbl_selected <- eg_gff_tbl %>% dplyr::select(seqnames, start, end, type, ID)
# Change names of columns to characters
# Remove mitochondria chromosome 'OX731682.1'
eg_gff_tbl_mutate <- eg_gff_tbl_selected %>% mutate(
  seqnames = as.character(seqnames),
  type = as.character(type)
) %>%
  filter(seqnames != 'OX731682.1')


##################################################################################################

### SPLITTING GENOME INTO BINS

# Put sequence into bins of length 1,000,000
# Assign 1,000,000 to variable 'bin_size' for simplicity in later functions
bin_size = 1e6

# Following function makes a column with a list of all bin names:
# Make empty table
bins <- tibble()
genome_index <- tibble(seqnames = names(genome_seq), contig_width = width(genome_seq))
# Loop through
for(i in seq_along(genome_index$seqnames)){
  #number each bin. Combine this number onto the end of bin name
  new_tbl <- tibble(seqnames = genome_index$seqnames[i], bin_no = 1:ceiling(genome_index$contig_width[i]/bin_size)) %>%
    mutate(bin_name = paste0(seqnames, "_", bin_no))
  #combine this new column into our empty table
  bins <- rbind(bins, new_tbl)
}

# This function removes mitochondria chromosome and selects for bin_name only
bins <- bins %>%
  filter(seqnames !="OX731682.1") %>%
  dplyr::select(bin_name)
# Generate column for start and end of each bin
eg_gff_bins <- eg_gff_tbl_mutate %>%
  mutate(
    bin_start = floor(start/bin_size)*bin_size + 1,
    bin_end = ceiling(start/bin_size)*bin_size,
  )


##################################################################################################

### CORRECT FOR ANY OVERLAPPING REGIONS

# Check for sequences that don't overlap bins
# Repeats that don't overlap bins will be less or equal to bin end 
# Filter for repeats that are not overlapping:
eg_gff_not_overlapping <- eg_gff_bins %>% filter(end<=bin_end)

# Check for sequences that do overlap bins
# Repeat end will be more than the bin end
eg_gff_overlapping <- eg_gff_bins %>% filter(end>bin_end)

# Split overlapping regions into left and right side of the bin edge
# Left
left_overlapping <- eg_gff_overlapping %>% 
  mutate(end = bin_end)
# Right
right_overlapping = eg_gff_overlapping %>%
  mutate(start = bin_end +1)
# Combine left and right overlap back together
eg_gff_tbl_binned <- rbind(eg_gff_not_overlapping, left_overlapping) %>%
  rbind(right_overlapping)


# Creates dataframe that numbers each bin 
eg_gff_tbl_binned_for_counting <- eg_gff_tbl_binned %>%
  mutate(bin_no = bin_end/bin_size, 
         bin_name = paste0(seqnames, "_", bin_no))
# Check number of repeats characterised by EarlGrey after sorting into bins
nrow(eg_gff_tbl_binned_for_counting) # 1079421
# Sort bins in ascending order of chromosome number and start codon number
eg_gff_tbl_binned_arranged <- eg_gff_tbl_binned %>% 
  arrange(seqnames, start)


# Remove '/' from the names of repeat types
# Group repeated data into the columns: seqnames, type and bin_start. 
# Generate new column of total length of each bin (bin_total)
# Output is three columns for bin_name, type and bin_total
eg_gff_tbl_binned_grouped <- eg_gff_tbl_binned_arranged %>% 
  mutate(type = sub("/.*", "", type)) %>%
  group_by(seqnames, type, bin_start) %>% 
  mutate(bin_total = sum(end - start +1)) %>%
  ungroup() %>%
  mutate(bin_name = paste0(seqnames, "_", bin_end/bin_size)) %>% #assign each bin to a name
  dplyr::select(bin_name, type, bin_total) %>%
  base::unique() 

# Collect all data with scaffold length 
# This code corrects our data by replacing NA output as numerical value of 0  
# and fixes length of last bin to its true length (as truncated in original graphical output)
eg_gff_tbl_binned_corrected <- eg_gff_tbl_binned_grouped %>%
  pivot_wider(names_from = type, values_from = bin_total) %>% #selects for column name 'type' frequency against each bin (bin name). Assigns NA which we can input into original data
  pivot_longer(!bin_name, names_to = "type", values_to = "bin_total") %>% #input data back into original data frame, so NA is now there 
  mutate(bin_total = ifelse(is.na(bin_total), 0, bin_total)) %>% #ifelse function assigns NA as numerical value zero, if not NA keep numerical value already there 
  mutate(seqnames = sub("_.*", "", bin_name),
         bin_no = as.double(sub(".*_", "", bin_name)),
         bin_end = bin_no*bin_size,
         bin_start = bin_end-bin_size+1) %>%
  inner_join(repeat_tbl) %>% #fixing length of last bin to true length
  mutate(bin_end = ifelse(bin_end < scaffold_length, bin_end, scaffold_length),
         bin_width = bin_end-bin_start+1) %>%
  mutate(percentage = 100 * bin_total/bin_width) #percentage of that particular type of TE, for that bin 


# This function shows last 6 things at end and start of dataframe
# As a check to confirm we have generated the correct set of data
tail(eg_gff_tbl_binned_corrected)
head(eg_gff_tbl_binned_corrected)


##################################################################################################

###  PLOTTING OUR DATA
# we now want to plot distribution of repeats across genomic bins
# we will plot distribution of repeats across four groups:
# (1) All Repeat elements
# Then will plot for subfamilies within repeat types:
# (2) Subfamilies of DNA transposons
# (3) Subfamilies of LINE elements
# (4) Subfamilies of LTR  


# (1) All Repeat Elements
# Select and group repeated data (duplicated rows)
eg_gff_tbl_binned_grouped <- eg_gff_tbl_binned_corrected %>% 
  dplyr::select(seqnames, type, bin_start, bin_total)
# Remove duplicated rows
for_plotting_transposons <- base::unique(eg_gff_tbl_binned_grouped) %>%
  mutate(seqnames = sub('OX731680.1','Chromosome 1', seqnames)) %>%
  mutate(seqnames = sub('OX731681.1','Chromosome 2', seqnames))
# Plot data
plot_transposons <- ggplot(for_plotting_transposons, mapping = aes(x = bin_start, y = bin_total, fill = type)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~seqnames, ncol = 1) 
# Label x and y axis, change font size and theme
plot_transposons <- plot_transposons + xlab('Chromosome Position') + ylab('Frequency') 
# Change text size and fix legend
plot_transposons <- plot_transposons + theme_bw() + theme(text = element_text(size=15), legend.background = element_rect(color = "azure4", linetype = "solid"))
# fix up labels in the legend
plot_transposons <- plot_transposons + scale_fill_discrete(name = "Repeat Type", labels = c('DNA transposons', 'LINE', 'LTR', 'PLE',
                                         'RC', 'SINE', 'Unknown'))
plot_transposons

# (2) DNA Transposons
# Filter for 'DNA' in type column from all our data 
eg_gff_tbl_binned_DNA <- eg_gff_tbl_binned_arranged %>% 
  filter(grepl('DNA', type)) %>%
  mutate(type = sub('.*/','', type)) %>%
  mutate(type = sub("-.*", '', type)) %>%
  group_by(seqnames, type, bin_start) %>% 
  mutate(bin_total = sum(end - start +1)) %>%
  ungroup()
# Get number of DNA transposons characterised by EarlGrey
nrow(eg_gff_tbl_binned_DNA) # 170889
# Select and group repeated data (duplicated rows)
eg_gff_tbl_grouped_DNA <- eg_gff_tbl_binned_DNA %>% 
  dplyr::select(seqnames, type, bin_start, bin_total) %>%
  mutate()
# Remove duplicated rows
for_plotting_DNA = base::unique(eg_gff_tbl_grouped_DNA) %>% 
  mutate(seqnames = sub('OX731680.1','Chromosome 1', seqnames)) %>%
  mutate(seqnames = sub('OX731681.1','Chromosome 2', seqnames))
# Plot data for DNA
plot_DNA <- ggplot(for_plotting_DNA, mapping = aes(x = bin_start, y = bin_total, fill = type)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~seqnames, ncol = 1) + 
  labs(fill = 'Subfamily')
# Label x and y axis, change font size and theme
plot_DNA <- plot_DNA + xlab('Chromosome Position') + ylab('Frequency')
# Change text size and fix legend
plot_DNA <- plot_DNA + theme_bw() + theme(text = element_text(size=15), 
                                          legend.background = element_rect(color = "azure4", linetype = "solid"))
# Final Plot for DNA transposons
plot_DNA


# 2.3) LINE elements
eg_gff_tbl_binned_line = eg_gff_tbl_binned_arranged %>% 
  filter(grepl('LINE', type)) %>%
  mutate(type = sub('.*/','', type)) %>%
  mutate(type = sub("-.*", '', type)) %>%
  group_by(seqnames, type, bin_start) %>% 
  mutate(bin_total = sum(end - start +1)) %>%
  ungroup()
# Get number of DNA transposons characterised by EarlGrey
nrow(eg_gff_tbl_binned_line) # 402608
# Select and group repeated data (duplicated rows)
eg_gff_tbl_grouped_line = eg_gff_tbl_binned_line %>% 
  dplyr::select(seqnames, type, bin_start, bin_total)
# Remove duplicated rows
for_plotting_line = base::unique(eg_gff_tbl_grouped_line) %>%
  mutate(seqnames = sub('OX731680.1','Chromosome 1', seqnames)) %>%
  mutate(seqnames = sub('OX731681.1','Chromosome 2', seqnames))
# Plot data
plot_line <- ggplot(for_plotting_line, mapping = aes(x = bin_start, y = bin_total, fill = type)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~seqnames, ncol = 1) + 
  labs(fill = 'Subfamily')
#label x and y axis, change font size and theme
plot_line <- plot_line + xlab('Chromosome Position') + ylab('Frequency')
# Change text size and fix legend
plot_line <- plot_line + theme_bw() + theme(text = element_text(size=15),
                                            legend.background = element_rect(color = "azure4", linetype = "solid"))
# Final Plot 
plot_line


# (3) LTR 
# Filter for 'LTR' in type column from all our data 
eg_gff_tbl_binned_LTR <- eg_gff_tbl_binned_arranged %>% 
  filter(grepl('LTR', type)) %>%
  mutate(type = sub('.*/','', type)) %>%
  mutate(type = sub("-.*", '', type)) %>%
  group_by(seqnames, type, bin_start) %>% 
  mutate(bin_total = sum(end - start +1)) %>%
  ungroup()
# Get number of DNA transposons characterised by EarlGrey
nrow(eg_gff_tbl_binned_LTR) # 11535
# Select and group repeated data (duplicated rows)
eg_gff_tbl_grouped_LTR = eg_gff_tbl_binned_LTR %>% 
  dplyr::select(seqnames, type, bin_start, bin_total)
# Remove duplicated rows
for_plotting_LTR = base::unique(eg_gff_tbl_grouped_LTR) %>%
  mutate(seqnames = sub('OX731680.1','Chromosome 1', seqnames)) %>%
  mutate(seqnames = sub('OX731681.1','Chromosome 2', seqnames))
# Plot data
plot_LTR <- ggplot(for_plotting_LTR, mapping = aes(x = bin_start, y = bin_total, fill = type)) + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~seqnames, ncol = 1) +
  labs(fill = 'Subfamily')
# Label x and y axis, change font size and theme
plot_LTR <- plot_LTR + xlab('Chromosome Position') + ylab('Frequency')
# Change text size and fix legend
plot_LTR <- plot_LTR + theme_bw() + theme_bw() + theme(text = element_text(size=15),
                                                       legend.background = element_rect(color = "azure4", linetype = "solid"))
plot_LTR







##################################################################################################

### ESTIMATING GC CONTENT

# Make granges to get entire genome sequence 
# Align genome sequence into bins 
genome_binned_ranges <- eg_gff_tbl_binned_corrected %>%
  dplyr::select(seqnames, bin_start, bin_end, bin_name) %>%
  base::unique() %>%
  dplyr::rename(start = bin_start, end = bin_end) %>%
  as_granges()

# Calculation for GC content
GC_content <- tibble(seqnames = character(), bin_name = character(), gc_count = integer(), at_count = integer(), n_count = integer())

for(i in seq_along(genome_binned_ranges$bin_name)){
  temp_ranges <- genome_binned_ranges[i]
  temp_string <- as.character(BSgenome::getSeq(genome_seq, temp_ranges))
  temp_counts <- table(strsplit(temp_string, split="")) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Var1 = as.character(Var1))
  temp_counts_tbl <- tibble(bin_name = genome_binned_ranges$bin_name[i],
                            gc_count = sum(temp_counts[temp_counts$Var1 %in% c('G', 'g', 'C', 'c'),]$Freq),
                            at_count = sum(temp_counts[temp_counts$Var1 %in% c('A', 'a', 'T', 't'),]$Freq),
                            n_count = sum(temp_counts[temp_counts$Var1 %in% c('N', 'n'),]$Freq)
                            )
  GC_content <- rbind(GC_content, temp_counts_tbl)
}

GC_content <- cbind(GC_content, (GC_content$gc_count)/(GC_content$at_count+GC_content$gc_count)*100) %>%
  as_tibble()

colnames(GC_content)[5] <- 'GC_content' 


## Plot GC content across genomic bins 

# This function adds seqnames, start and end to our dataframe
# So the plot for distribution of GC content is consistent with repeat plots (see above)
GC_content_grouped <- GC_content %>% 
  mutate(seqnames = sub("_.*", "", bin_name)) %>%
  mutate(bin_no = as.double(sub(".*_", "", bin_name)),
         bin_no = as.double(sub(".*_", "", bin_name)),
         bin_end = bin_no*bin_size,
         bin_start = bin_end-bin_size+1) %>%
  mutate(bin_name = sub("OX73168..1_", "", bin_name)) %>%
  mutate(bin_name = as.integer(bin_name)) %>%
  dplyr::select(seqnames, bin_name, bin_start, gc_count, n_count, GC_content)
# rename seqnames to Chromosome names
for_plotting_GC <- GC_content_grouped %>%
  mutate(seqnames = sub('OX731680.1','Chromosome 1', seqnames)) %>%
  mutate(seqnames = sub('OX731681.1','Chromosome 2', seqnames))
# plotting
plot_GC <- ggplot(for_plotting_GC, mapping = aes(x = bin_start, y = GC_content, fill = GC_content)) + 
    geom_histogram(stat = "identity", width = 1e6) +
    facet_wrap(~seqnames, ncol = 1) + 
    labs(fill = "GC content")
plot_GC <- plot_GC + scale_fill_gradient(low="midnightblue", high="slateblue1")
#label x and y axis, change font size and theme
plot_GC <- plot_GC + xlab('Chromosome Position') + ylab('GC content')
plot_GC <- plot_GC + theme_bw() + theme(text = element_text(size=15))
plot_GC



##################################################################################################

### ESTIMATING REPEAT COUNT
# We want to estimate the distribution of repeat density along the genome
# I will estimate type Density as the COUNT of each repeat type in each bin 
# Counts estimated for each repeat type: DNA, LINE, LTR, PLE, RC, SINE, Unknown

# This function converts factors into characters
# Output: table with the columns bin_name, type and count
eg_gff_tbl_binned_for_counting <- eg_gff_tbl_binned_for_counting %>%
  mutate(type = sub(pattern = '-.*', replacement = '', x = type))
# Make our data frame wide
eg_gff_tbl_counted_wide <- base::table(eg_gff_tbl_binned_for_counting$bin_name, eg_gff_tbl_binned_for_counting$type) %>%
  base::as.data.frame() %>%
  as_tibble() %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2)) %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  dplyr::rename(bin_name = Var1) 
# Make our dataframe long
eg_gff_tbl_counted_long <- left_join(x = bins, y = eg_gff_tbl_counted_wide) %>%
  pivot_longer(!bin_name, names_to = 'type', values_to = 'Count') 
# Double check there is no NA
filter(.data = eg_gff_tbl_counted_long, is.na(Count))

# Remove subfamilies by removing '/' from type name
# Remove subfamilies so we can merge the same type together in one bin and so the columns are consistent in size 
eg_gff_tbl_counted_filtered <- eg_gff_tbl_counted_long %>%
  mutate(type = sub(pattern = '/.*', replacement = '', x = type)) %>%
  group_by(bin_name, type) %>%
  summarise(Total = sum(Count)) 


# (1) DNA
DNA_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('DNA', type)) 

# (2) LINE
LINE_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('LINE', type))

# (3) LTR
LTR_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('LTR', type)) 

# (4) PLE
PLE_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('PLE', type)) 

# (5) RC
RC_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('RC', type)) 

# (6) SINE
SINE_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('SINE', type)) 

# (7) Unknown
Unknown_counted <- eg_gff_tbl_counted_filtered %>%
  filter(grepl('Unknown', type)) 

    
#### get total number of repeat types characterised by EarlGrey when split into genomic bins 
sum(DNA_counted$Total) #170889 
sum(LINE_counted$Total) #402608
sum(LTR_counted$Total) #11535
sum(PLE_counted$Total) #612
sum(RC_counted$Total) #61589
sum(SINE_counted$Total) #6754
sum(Unknown_counted$Total) #425434
# sum of the total of all repeat types 
sum(170889 + 402608 + 11535 + 612 + 61589 + 6754 + 425434) #1079421
nrow(eg_gff_tbl_binned_for_counting) #1079421





