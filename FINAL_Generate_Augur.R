# This script prepares the working environment for phylogenetic analysis by
# clearing old data,  setting file paths, loading essential packages, and
# generating citation files. It processes consensus sequence data and
# phylogenetic information for influenza strains H3N2 and H1N1, aligning
# sequences, filtering variants, and preparing data for visualization and
# further analysis.
# Author: Ashley Sobel Leonard
# Date: 5/1/2024

# Clear environment to get rid of old files and get paths  -------------------
rm(list = ls())
R_path <- getwd()

# Get paths to fasta consensusns sequences and where the phylogeny files are located
Fasta_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut","Consensus_seqs")
Phylogeny_path <- file.path(R_path,"CHOA_Phylogenies")

# Load packages ----- 
library(tidyverse)
library(phylotools) # use for read.fasta
library(seqinr)
library(data.table)


# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("tidyverse", "phylotools","seqinr","data.table")

# Initialize an empty character vector to hold citations
all_citations <- character(0)

# Loop through each package and get its citation in BibTeX format
for (pkg in package_list) {
  # Fetch the citation information for the package
  citation_info <- citation(pkg)
  
  # Convert the citation to BibTeX format
  bibtex_info <- utils::toBibtex(citation_info)
  
  # Extract publication year from the citation information
  pub_year <- sapply(citation_info, function(x) x$year[1])
  
  # Generate a unique key using the package name and publication year
  keys <- paste0(tolower(pkg), pub_year)
  
  # Replace the first occurrence of the comma with a comma followed by the key
  bibtex_info <- sub("@(\\w+)\\{,", paste0("@\\1{", keys, ","), bibtex_info)
  
  # Add the modified citation to the list
  all_citations <- c(all_citations, bibtex_info)
}

# Combine all citations into a single string
all_citations_text <- paste(all_citations, collapse = "\n\n")
# Write the combined citations to a .bib file
writeLines(all_citations_text, con = file.path("ScriptCitations","Generate_Augur.bib"))


# Define functions -----------
# Function to write a fasta file, typically used for exporting masked consensus sequences
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

GetFasta <- function(PipelineOutput_path) {
  path_fasta <- paste0(PipelineOutput_path,"/Consensus"); path_fasta
  fasta_names <- list.files(path_fasta, pattern = "\\.fasta"); fasta_names
  fasta_info <- list(path_fasta=path_fasta,fasta_names=fasta_names); fasta_info
}

# Add "whose" data to CHOP metadata file
adjust_metadata <- function(metadata) {
  # Adjust the metadata
  adjusted_metadata <- metadata %>%
    mutate(whose = case_when(
      # If 'name' contains 'CHOA', 'whose' should be 'CHOP'
      stringr::str_detect(name, "CHOA") ~ "CHOP",
      
      # If 'name' contains 'Pennsylvania', 'whose' should be 'PA'
      stringr::str_detect(name, "Pennsylvania") ~ "PA",
      
      # If the name contains neither 'CHOA' nor 'Pennsylvania', 'whose' should be 'USA'
      TRUE ~ "USA"
    ))
  
  # Return the adjusted metadata
  return(adjusted_metadata)
}

# Load data -----
CHOA_metadata <- fread(file.path(Phylogeny_path,'FINAL_CHOA_merged_data.csv'))

# Load "Bad_IDs_All" - the samples that were excluded due to poor replicability
Bad_IDs_All <- read.csv(file.path(Phylogeny_path,"CHOA_Bad_IDs_All.csv"))

# Remove the "Bad_IDs"
CHOA_metadata <- CHOA_metadata[!SubjectID %in% Bad_IDs_All]

# Extract sample names 
H3N2_name <- CHOA_metadata[Strain == "H3N2"]$SubjectID
H1N1_name <- CHOA_metadata[Strain == "H1N1"]$SubjectID

# Iterate over each of the H3 samples to load the sequence
H3n <- 1
H3_HA <- data.frame(seq.name=character(),seq.text=character())
for (H3n in 1:length(H3N2_name)){
  H3_name <- paste0(H3N2_name[H3n],"_consensus.fasta");H3_name
  H3_path <- file.path(Fasta_path,H3_name)
  H3_all_fasta <- phylotools::read.fasta(H3_path) # Phylotools package
  H3_HA_tmp <- H3_all_fasta[4,]
  H3_HA <- rbind(H3_HA,H3_HA_tmp)
}
rownames(H3_HA) <- NULL

# Write the alignment to a fasta file 
new_seq_name <- gsub(pattern = "_HA",replacement = "",x = H3_HA$seq.name)
new_fasta_df <- dplyr::data_frame(name=new_seq_name,seq=H3_HA$seq.text) #This is weird formatting, but it's the only way it works
writeFasta(data = new_fasta_df,filename = file.path(Phylogeny_path,"CHOA_H3N2_HA_seq.fasta"))

rm(new_seq_name)

# Repeat for H1 
H1n <- 1
H1_HA <- data.frame(seq.name=character(),seq.text=character())
for (H1n in 1:length(H1N1_name)){
  H1_name <- paste0(H1N1_name[H1n],"_consensus.fasta");H1_name
  H1_path <- file.path(Fasta_path,H1_name)
  H1_all_fasta <- phylotools::read.fasta(H1_path) # Phylotools package
  H1_HA_tmp <- H1_all_fasta[4,]
  H1_HA <- rbind(H1_HA,H1_HA_tmp)
}
H1_HA$seq.name
new_seq_name <- gsub(pattern = "_HA",replacement = "",x = H1_HA$seq.name)
rownames(H1_HA) <- NULL

new_fasta_df <- dplyr::data_frame(name=new_seq_name,seq=H1_HA$seq.text) #This is weird formatting, but it's the only way it works
writeFasta(data = new_fasta_df,filename = file.path(Phylogeny_path,"CHOA_H1N1_HA_seq.fasta"))

# Load files GISAID files
GISAID_subset_metadata <- read.csv(file.path(Phylogeny_path,"GISAID_H3N2_subset.csv"),row.names = NULL)
GISAID_subset_metadata <- GISAID_subset_metadata[,-1]; head(GISAID_subset_metadata)

# Now we will load the fasta files for each of the samples and get the metadata
CHOA_H3N2_HA_seq <- phylotools::read.fasta(file.path(Phylogeny_path,"CHOA_H3N2_HA_seq.fasta")) # Phylotools package

# Load CHOP metadata metadata 
CHOA_tree_metadata <- read.csv(file.path(Phylogeny_path,"CHOA_SeqDate.csv"))
CHOA_tree_metadata <- CHOA_tree_metadata[,-1]; head(CHOA_tree_metadata)
colnames(CHOA_tree_metadata) <- c("Sample","Collection")


# Generate Augur files for H3N2 -------------------
# Remove duplicate sequences
CHOA_H3N2_HA_seq <- CHOA_H3N2_HA_seq[!duplicated(CHOA_H3N2_HA_seq$seq.text),]
head(CHOA_H3N2_HA_seq$seq.name)

# Add the corresponding metadata (sampling date)
nmeta <- 1
CHOA_H3N2_HA_metadata <- data.frame()
for (nmeta in 1:nrow(CHOA_H3N2_HA_seq)){
  nmeta_tmp <- CHOA_H3N2_HA_seq$seq.name[nmeta]
  nmeta_tmp <- gsub("_HA","",nmeta_tmp); nmeta_tmp
  meta_data_tmp <- CHOA_tree_metadata[CHOA_tree_metadata$Sample == nmeta_tmp,]; meta_data_tmp
  meta_data_tmp <- meta_data_tmp[1,]; meta_data_tmp
  meta_data_tmp$Collection <- as.Date(meta_data_tmp$Collection, format = "%m/%d/%y")
  CHOA_H3N2_HA_metadata <- rbind(CHOA_H3N2_HA_metadata, meta_data_tmp)
}

# Combine the CHOP and GISAID metadata and save as a csv file 
CHOP_H3N2_metadata <- data.frame(name=CHOA_H3N2_HA_metadata$Sample,date=CHOA_H3N2_HA_metadata$Collection)
GISAID_H3N2_metadata <- data.frame(name=GISAID_subset_metadata$seq.name,date=GISAID_subset_metadata$collection_date)
CHOA_H3N2_HA_phylo_data <- rbind(CHOP_H3N2_metadata,GISAID_H3N2_metadata); head(CHOA_H3N2_HA_phylo_data)
CHOA_H3N2_HA_phylo_data <- adjust_metadata(CHOA_H3N2_HA_phylo_data)
write.csv(CHOA_H3N2_HA_phylo_data,file.path(Phylogeny_path,"Final_CHOA_H3N2_HA_metadata.csv"))

# Combine the CHOP and GISAID sequence data and save as a fasta files
GISAID_subset_fasta <- data.frame(seq.name= GISAID_subset_metadata$seq.name, seq.text = GISAID_subset_metadata$seq.text)
CHOA_H3N2_HA_phyloseq <- rbind(CHOA_H3N2_HA_seq,GISAID_subset_fasta); head(CHOA_H3N2_HA_phyloseq)
CHOA_H3N2_HA_phyloseq <- dplyr::data_frame(name = CHOA_H3N2_HA_phyloseq$seq.name,
                                           seq = CHOA_H3N2_HA_phyloseq$seq.text,
                                           other = NULL)
writeFasta(CHOA_H3N2_HA_phyloseq,file.path(Phylogeny_path,"CHOA_H3N2_HA_phyloseq.fasta"))


# Based on the initial tree, H3N2_CHOP_tree_1.nex (dowloaded time tree from
# nextstain), I am removing the following sequences for visualization purposes
# to allow the tree root to be rooted in a more recent time point. These sequences are: 
# A/Maryland/08/2012
# A/Maryland/05/2012
# A/Maryland/04/2012
# A/Maryland/09/2012
# A/Iowa/09/2011
# A/Iowa/08/2011
# A/West_Virginia/06/2011
# A/Utah/10/2012

# Step 1: Generate a list of the sequences to be removed named "CHOP_H3N2_RemovalSeqs"
# CHOP_H3N2_RemovalSeqs <- c("A/Maryland/08/2012", "A/Maryland/05/2012", "A/Maryland/04/2012",
#                           "A/Maryland/09/2012", "A/Iowa/09/2011", "A/Iowa/08/2011", 
#                            "A/West_Virginia/06/2011", "A/Utah/10/2012")
# 
# # Step 2: Remove the sequences in "CHOP_H3N2_RemovalSeqs" from CHOA_H3N2_HA_phyloseq
# CHOA_H3N2_HA_phyloseq_2 <- CHOA_H3N2_HA_phyloseq[!CHOA_H3N2_HA_phyloseq$name %in% CHOP_H3N2_RemovalSeqs, ]
# head(CHOA_H3N2_HA_phyloseq_2)
# 
# # Step 3: Remove the sequences in "CHOP_H3N2_RemovalSeqs" from CHOA_H3N2_HA_phylo_data
# CHOA_H3N2_HA_phylo_data_2 <- CHOA_H3N2_HA_phylo_data[!CHOA_H3N2_HA_phylo_data$name %in% CHOP_H3N2_RemovalSeqs, ]
# 
# # Step 4: Save the revised sequences as CHOA_H3N2_HA_phyloseq_2 using the function writeFasta
# writeFasta(CHOA_H3N2_HA_phyloseq_2, "CHOA_H3N2_HA_phyloseq_2.fasta")
# 
# # Step 5: Save the revised metadata as CHOA_H3N2_HA_phylo_data_2
# # Assuming you want this saved as a .csv file
# write.csv(CHOA_H3N2_HA_phylo_data_2, "CHOA_H3N2_HA_metadata_2.csv")


# Repeat for H1N1 -------------------
# Load files GISAID files
GISAID_subset_metadata <- read.csv(file.path(Phylogeny_path,"GISAID_H1N1_subset.csv"),row.names = NULL)
GISAID_subset_metadata <- GISAID_subset_metadata[,-1]; head(GISAID_subset_metadata)

# Now we will load the fasta files for each of the samples and get the metadata
CHOA_H1N1_HA_seq <- phylotools::read.fasta(file.path(Phylogeny_path,"CHOA_H1N1_HA_seq.fasta")) # Phylotools package

# Remove duplicate sequences
CHOA_H1N1_HA_seq <- CHOA_H1N1_HA_seq[!duplicated(CHOA_H1N1_HA_seq$seq.text),]
nmeta <- 1
CHOA_H1N1_HA_metadata <- data.frame()

# Generate metadata files for H3N2
for (nmeta in 1:nrow(CHOA_H1N1_HA_seq)){
  nmeta_tmp <- CHOA_H1N1_HA_seq$seq.name[nmeta]
  nmeta_tmp <- gsub("_HA","",nmeta_tmp); nmeta_tmp
  meta_data_tmp <- CHOA_tree_metadata[CHOA_tree_metadata$Sample == nmeta_tmp,]; meta_data_tmp
  meta_data_tmp <- meta_data_tmp[1,]; meta_data_tmp
  meta_data_tmp$Collection <- as.Date(meta_data_tmp$Collection, format = "%m/%d/%y")
  CHOA_H1N1_HA_metadata <- rbind(CHOA_H1N1_HA_metadata, meta_data_tmp)
}

# Combine the CHOP and GISAID metadata and save as a csv file 
CHOP_H1N1_metadata <- data.frame(name=CHOA_H1N1_HA_metadata$Sample,date=CHOA_H1N1_HA_metadata$Collection)
GISAID_H1N1_metadata <- data.frame(name=GISAID_subset_metadata$seq.name,date=GISAID_subset_metadata$collection_date)
CHOA_H1N1_HA_phylo_data <- rbind(CHOP_H1N1_metadata,GISAID_H1N1_metadata); head(CHOA_H1N1_HA_phylo_data)
CHOA_H1N1_HA_phylo_data <- adjust_metadata(CHOA_H1N1_HA_phylo_data)
write.csv(CHOA_H1N1_HA_phylo_data,file.path(Phylogeny_path,"Final_CHOA_H1N1_HA_metadata.csv"))

# Combine the CHOP and GISAID sequences and save as a fasta file
GISAID_subset_fasta <- data.frame(seq.name=GISAID_subset_metadata$seq.name, seq.text=GISAID_subset_metadata$seq.text)
CHOA_H1N1_HA_phyloseq <- rbind(CHOA_H1N1_HA_seq,GISAID_subset_fasta); head(CHOA_H1N1_HA_phyloseq)
CHOA_H1N1_HA_phyloseq <- dplyr::data_frame(name = CHOA_H1N1_HA_phyloseq$seq.name,
                                           seq = CHOA_H1N1_HA_phyloseq$seq.text,
                                           other = NULL)

writeFasta(CHOA_H1N1_HA_phyloseq,file.path(Phylogeny_path,"CHOA_H1N1_HA_phyloseq.fasta"))

