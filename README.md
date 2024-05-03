![Image the University of Edinburgh logo.](https://www.logo.wine/a/logo/University_of_Edinburgh/University_of_Edinburgh-Logo.wine.svg)
# Honours Research Project 
This page contains all custom script created in `R` for my Genetics Honours Project. This script was used to annotate repeat elements and genic content across the whole-genome assembly of the Cottony Cushion Scale Insect, _Icerya purchasi_. 
# Contents
1.  'EarlGrey' Transposabe Element Annotation
2.  'AUGUSTUS' gene annotation
3.   Divergence of Repeats
4.   Locating Repeats in Gene Regions


## 'EarlGrey' Transposable Element Annotation
**EarlGrey** is a Transposable Elements (TE) annotation pipeline, available at: https://github.com/TobyBaril/EarlGrey.  
EarlGrey is used to identify, curate and annotate TEs across the input genome of _Icerya purchasi_.    
   
The `R` script in this file, takes the gff output of 'EarlGrey' to annotate and plot repeat content across genomic bins of 100 kb. Estimations of GC content is also contained within this script. 


## 'AUGUSTUS' gene annotation
**AUGUSTUS** is a gene prediction programme. 
The ´R` script in this file, takes the gff output of the AUGUSTUS gene prediciton programme, and annotates gene content across genomic bins of 100 kb. This script also provides estimations of the gene count, gene percentage and range of introns, exons and regulatory regions.

## Divergence of Repeats
This script takes and plots the estimations of K2P distance as estimated by 'EarlGrey'. 

## Locating Repeats in Gene Regions 
This script determines where repeats lie within genomic regions. 
