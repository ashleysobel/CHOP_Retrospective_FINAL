# CHOP_Phylogenies 
This is a subfolder of the CHOP_Retrospective_IAV_Diversity repo that contains scripts used to generate the phylogenies. There are 2 primary R scripts and 2 text files. 

The R scripts are as follows: 

1) get_GISAID_subsets.R

This script is used to subsample sequence data downloaded from GISAID. The code imports relevant packages and datasets, defines functions, filters sequences based on specific time frames and randomly # selects sequences for further analysis. It handles two types of data - a dataset of sequences from the United States (both the metadata and associated fasta file) and sequences from Pennsylvania (both the metadata and associated fasta file). After loading and cleaning these datasets, the code generates a new dataset combining sequences from both. The code is processes both H3N2 and H1N1 sequence. 

This script uses the following packages: 
stringr
phylotools
dplyr
stringi
lubridate

2) FINAL_Generate_Augur.R

This script combines the subsampled GISAID data with the CHOP sequences for H3N2 and H1N1 and generates a fasta file and metadata file for each. It also adds the annotation "whose" that specifies where the sequence was obtained, i.e. United States vs. Pennsylvania for the GISAID sequences and CHOP for the CHOP sequences. 

This script uses the following packages: 
stringr
phylotools
dplyr
stringi
lubridate
data.table
seqinr


The output files: Final_CHOA_H3N2_HA_metadata.csv, CHOA_H3N2_HA_phyloseq.fasta, Final_CHOA_H1N1_HA_metadata.csv, and CHOA_H1N1_HA_phyloseq.fasta are copied into the local nextstrain directory to generate a phylogeny using Nextstrain. 

3-4) CHOA_NextStrain_Script_H1N1.txt and CHOA_NextStrain_Script_H3N2.txt

This file contains the commands for running Nextstrain for H1N1 and H3N2. 