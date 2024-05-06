# CHOP_Retrospective_GitHub
Contains the code used for: ¬Within-Host Influenza Viral Diversity in the Pediatric Population as a Function of Age, Vaccine and Health Status


## Overview

This document describes the analyses used to evaluate factors modulate intrahost viral diversity in the CHOA Retrospective flu study. The data was generated from next-generation sequencing of clinical isolates from the 2017-2018 flu season at CHOP. For post-processing of the CHOP sequence, please run the following code in this order:

1) Organize_PipelineOutput.R: This script takes in the output from FluPipeline, organizes it, generates complementary file structure in the R working directory, and copies the relevant files into R workspace. It also generates pileup files for the CHOP sequences. 

For this code to work, it is ESSENTIAL that the sequences following the following naming convention: 
Prefix-###_D#V#

Here, the prefix is CHOP, the numbers range from 001 - 197, D is the replicate #, and V is the version number. D represents the technical replicate, which denotes an individual RNA extraction from a clinical sample. V represents the version, which is used to identify repeated sequences run from the same RNA extaction. For example, CHOP-001_D1V1, indicates the sample was from CHOP-001, was from the initial RNA extraction, and the first time that extraction was sequenced. Conversely, CHOP-057-D3V1, represents a sample from CHOP-057 which was the 3rd time the RNA was extacted but the first time that extraction was sequenced. 

This script requires the following packages: 
tidyverse
Rsamtools


2) SampleQC.R: This script performs additional quality control steps for the CHOP sequences an cateogorizes all of the variants identified in the samples. The required input for this script includes: the references sequence (and a csv file with additional information on the CDS/protein coding regions for those references), a reference key that contains the PB2 accession number(s) for your references, the run dates for the sequences, and the additional quality control metrics you have chosen for minimim coverage, mininum Phred quality score, and mininum mapping score. 

This script requires the following packages: 
stringr
phylotools
tidyverse
stringi
seqinr
tidysq
patchwork 

3) Compile_data.R: This script combines the data from different sequencing runs, incorporates the clean clinical metadata, and identifies technical replicates and ensures there is concordance between the technical replicates (excluding those with poor concordance). This also generates Figure 1. 

This script requires the following packages: 
tidyverse
phylotools
gridExtra
seqinr
GenomicRanges
data.table

Once the sequence analysis is complete, the scripts for the data analysis can be completed. These scripts include: 
  1. Data_Analysis: Complete the data analysis for the CHOP sequences for: Table 1, Figures 3, and S2.  
  2. GenerateMutation_Table: Expands the information present in the table of mutations to create Table S2.
  3. Generate_Dandelion_Plots - organizes the data for Table 2 and Figures 4 - 5
  4. Generate_Coverage_Plots - generates coverage plots, organized by run, for Figures S2 - S8
  4. Calculate_NucletideDiversityStatistic - calculate the nucleotide diversity and creates the plot for Figure S3
  5. Partition_Compiled_info - generates fasta files containing sequences organized by strain and segment
  6. Generate_MixedInfection_plot - generated mixed infeciton plots for Figure S1
  7. Generate_SRA_datasheets - this code generates the citation document for Table S1, the acnknowledgement table for GISAID
  8. get_GISAID_subsets - selects sequence from GISAID sequences 
  9. FINAL_Generate_Augur - generates data structures for running nextstrain


# CHOP_Retrospective_FINAL
