![Image the University of Edinburgh logo.](https://www.ed.ac.uk/sites/all/themes/uoe/assets/logo.png)
# Honours Research Project 
This page contain all custom script created in `R` for my Genetics Honours Project.   
My project investigates Androdioecy and characterisation of repeat elements contribution to genome size in the Cottony Cushion scale insect (_Icerya purchasi_).

# Contents
1. Annotation Resources 
2. 'Transposabe Element Annotation' script
3. 'Genomic Annotation' script
4. 'Divergence of Repeats' script
5. 'Locating Repeats in Gene Regions' script
6. Programmes used

# Annotation Resources
- **EarlGrey** is a Transposable Elements (TE) annotation pipeline which I used to identify, curate and annotate TEs across the whole genome assembly of I. _purchasi_, available at: https://github.com/TobyBaril/EarlGrey. The EarlGrey add-on '_EarGreyDivergenceCalc_' was included to estimate the divergence of repeats (quantified as the Kimura80 distance). 
     
- **AUGUSTUS** is a gene prediction programme which I used to classify genes and genomic across the whole genome assembly of I. _purchasi_.
  

## '<ins>Transposable Element Annotation</ins>' script
This `R` script takes the gff output from EarlGrey, and annotates repeat content across genomic bins of 100 kb. Estimations of GC content is also contained within this script. 


## '<ins>Genomic Annotation</ins>' script  
This `R` script takes the gff output of AUGUSTUS, and annotates gene content across genomic bins of 100 kb. This script also provides estimations of the gene count, gene percentage and range of introns, exons and regulatory regions.


## '<ins>Divergence of Repeats</ins>' script
This `R` script takes the gff output from 'EarlGrey' to plot divergence of different repeat types across the whole-genome assembly.  

## '<ins>Locating Repeats in Gene Regions</ins>' script
This `R` script takes both the gff from 'EarlGrey' and the gff from AUGUSTUS, to map where repeats lie within genomic regions.    
This script must be run last on ´R´, as requires all code from three previous scripts: (i) 'Transposabe Element Annotation' script, (ii) 'Genomic Annotation' script and (iii) 'Divergence of Repeats' script. 
