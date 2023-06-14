# Nuclear-retention-in-response-to-hypoxia-in-plants
In this workflow is used to analyse the differentially expressed genes and alternative splicing events, the following steps were followed: 

![project flowchart drawio (2)](https://github.com/de-Laat/Nuclear-retention-in-response-to-hypoxia-in-plants/assets/127954517/e38d2e1a-9b6a-4a3e-8a1e-1245fcfb46c9)




The SRA files were downloaded from the NCBI database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122804).

![afbeelding](https://github.com/de-Laat/Nuclear-retention-in-response-to-hypoxia-in-plants/assets/127954517/83ea38df-7b28-46c1-9d13-e70d6289f97c)

# Pre-processing
## fastq-dump
The SRA files are converted to FASTQ files using fastq-dump

## Trimming data
The FASTQ files are trimmed using FASTX-toolkit

## Fastqc
The quality of trimmed FASTQ files is checked using fast QC

# Alternative splicing analysis
## alligning kallisto
the files are alligned to RTD3 Fasta files using kallisto

## alligning kallisto
the alligned reads are counted using kallisto

## SUPPA2
Splicing events are detected by comparing the counted files to RTD3 GFF reference files

## extract <= 0.05 p-val dpsi files (suppa)
Extract significant splicing events from dpsi files

## Made in excel
A coverage plot is made from the output of Suppa

## Performed in excel
Methods of alternative splicing are identified

# DEG analysis
## alligning Hisat2
the files are alligned to Tair10 Fasta files using Hisat2

## HTseq
counting of alligned reads is perfomed on the resulting BAM file in combination with a Tair10 GFF3 file using HTseq count

## DEG/Scripts/Analyse_DEG.R
DEG analysis is performed on the resulting counts file

## DEG/Scripts/Clustering_Heatmap.R
A heatmap is made on the DEG analysis

## DEG/Scripts/GEA.R
A GEA analysis is performed on the output of the DEG analysis

# Combining output
The output from the different analyses are compared using MEME. This has not been done yet, but can be done in follow-up research.
