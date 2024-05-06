# Code description: The provided script is a comprehensive bioinformatics
# pipeline for dealing with sequence data related to the H3N2 strain of the
# Influenza virus. The code imports relevant packages and datasets, defines
# functions, filters sequences based on specific time frames and randomly
# selects sequences for further analysis. It handles two types of data - a
# global dataset (GISAID_H3N2_2017_FilteredMetadata.csv & associated fasta file)
# and a more localized dataset focusing on Pennsylvania
# (GISAID_H3N2_2017_PA_FilteredMetadata.csv & associated fasta file). After
# loading and cleaning these datasets, the code generates a new dataset
# combining sequences from both. It specifically focuses on the H3N2 strain of
# the Influenza virus and does not handle H1N1 data.
# Author: Ashley Sobel Leonard
# Date: 5/1/2024

# Load packages ----- 
library(stringr) #use for str_replace_all
library(phylotools) # use for read.fasta
library(dplyr) # use for select
library(stringi) # use for stri_extract_last
library(lubridate)

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("stringr", "phylotools","dplyr","stringi","lubridate")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","get_GISAID_subsets_citations.bib"))


# Prepare environment ----- 
rm(list = ls())
R_path <- getwd()
Phylogeny_path <- file.path(R_path,"CHOA_Phylogenies")

# Load reference data ----
# Load Flu season reference file
Flu_seasons <- read.csv(file.path(Phylogeny_path,"CHOP_Retro_Flu_seasons.csv"))

# Convert the season time points to dates
Flu_seasons$Start <- as_date(Flu_seasons$Start)
Flu_seasons$End <- as_date(Flu_seasons$End)
nSeasons <- nrow(Flu_seasons) - 1; nSeasons


# Define functions ---- 
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

# Flu_seasons_this <- Flu_seasons_H3N2_PA; GISAID_combined <- GISAID_H3N2_PA_combined; n_seq_general <- n_seq_PA; n_seq_SampleYear <- n_seq_PA_SampleYear;
Get_SequenceSubset <- function(Flu_seasons_this,GISAID_combined,n_seq_general,n_seq_SampleYear,nSeason_SampleYear){
  # Initialize data frames to hold randomly sampled sequences
  compiled_data <- data.frame()
  ns <- 1
  for (ns in 1:nrow(Flu_seasons_this)){
    Season_start <- Flu_seasons_this$Start[ns]
    Season_end <- Flu_seasons_this$End[ns]
    
    # Subset the GISAID sequence for the year 
    selected_data_all <- subset(GISAID_combined, GISAID_combined$collection_date > Season_start & GISAID_combined$collection_date < Season_end)
    all_seq <- nrow(selected_data_all); all_seq
    
    # Set the number of sequences to subsample for the year
    if (ns == nSeason_SampleYear){
      n_seq <- n_seq_SampleYear
    } else {
      n_seq <- n_seq_general
    }
    
    # Randomly subsample the sequences 
    if (n_seq <= all_seq){ 
      selected_data <- selected_data_all[sample(1:nrow(selected_data_all),n_seq,replace = FALSE),]; selected_data
    } else {
      selected_data <- selected_data_all
    } 
    
    # Add the selected sequences and metadata to the running list
    compiled_data <- rbind(compiled_data,selected_data); head(compiled_data)
    rm(Season_end, Season_start, selected_data_all, selected_data)
  }
  
  # Reset row names for compiled sequences and metadata
  rownames(compiled_data) = seq(length=nrow(compiled_data))
  
  # Check for duplicates
  duplicates_seq <- duplicated(compiled_data$seq.text)
  
  # Remove duplicates
  compiled_data <- compiled_data[!duplicates_seq,]
  return(compiled_data)
}


metadata_file <- "GISAID_H3N2_US_2017_FilteredMetadata.csv"
fasta_file <- "GISAID_H3N2_USA_2017_filtered.fasta"

load_and_combine_data <- function(metadata_file, fasta_file) {
  metadata <- read.csv(file.path(Phylogeny_path,metadata_file), row.names = NULL)
  metadata <- metadata[,-1]
  metadata$collection_date <- as_date(metadata$collection_date)
  fasta_data <- phylotools::read.fasta(file.path(Phylogeny_path,fasta_file))
  
  if (sum(fasta_data$seq.name != metadata$seq.name) > 0) {
    stop("The sequences from the fasta files and metadata are not equivalent/in the same")
  } else {
    combined_data <- cbind(metadata, seq.text = fasta_data$seq.text)
  }
  
  return(combined_data)
}

# nSeason_Start <- 2
# nSeason_End <- 8
# combined_data <- GISAID_H3N2_USA_combined
# n_seq_general <- 15
# n_seq_SampleYear <- 30 
# nSeason_SampleYear <- 5

generate_subset <- function(Flu_seasons, nSeason_Start, nSeason_End, combined_data, n_seq_general, n_seq_SampleYear, nSeason_SampleYear) {
  Flu_seasons_subset <- Flu_seasons[nSeason_Start:nSeason_End,]
  row.names(Flu_seasons_subset) <- NULL
  
  subset_data <- Get_SequenceSubset(Flu_seasons_subset, combined_data, n_seq_general, n_seq_SampleYear, nSeason_SampleYear)
  return(subset_data)
}

remove_duplicates <- function(data) {
  duplicates <- duplicated(data$seq.text)
  unique_data <- data[!duplicates,]
  return(unique_data)
}

# # Run code for H3N2 ------
# # Load and combine H3N2 USA data
# GISAID_H3N2_USA_combined <- load_and_combine_data("GISAID_H3N2_US_2017_FilteredMetadata.csv", "GISAID_H3N2_USA_2017_filtered.fasta")
# 
# # Generate H3N2 USA subset
# GISAID_H3N2_USA_subset <- generate_subset(Flu_seasons, 2, 8, GISAID_H3N2_USA_combined, 15, 30, 5)
# 
# # Load and combine H3N2 PA data
# GISAID_H3N2_PA_combined <- load_and_combine_data("GISAID_H3N2_PA_2017_FilteredMetadata.csv", "GISAID_H3N2_PA_2017_filtered.fasta")
# 
# # Generate H3N2 PA subset
# GISAID_H3N2_PA_subset <- generate_subset(Flu_seasons, 8, 8, GISAID_H3N2_PA_combined, 25, 50, 1)
# 
# # Combine and remove duplicates for H3N2
# GISAID_H3N2_subset <- rbind(GISAID_H3N2_USA_subset, GISAID_H3N2_PA_subset)
# GISAID_H3N2_subset <- remove_duplicates(GISAID_H3N2_subset)
# 
# # Save the subsetted GISAID sequences
# write.csv(GISAID_H3N2_subset,file.path(Phylogeny_path,"GISAID_H3N2_subset.csv"))

# Run code for H1N1 -----
# Load and combine H1N1 USA data
GISAID_H1N1_USA_combined <- load_and_combine_data("GISAID_H1N1_US_2017_FilteredMetadata.csv", "GISAID_H1N1_USA_2017_filtered.fasta")

# Generate H1N1 USA subset
GISAID_H1N1_USA_subset <- generate_subset(Flu_seasons, 4, 8, GISAID_H1N1_USA_combined, 15, 30, 5)

# Load and combine H1N1 PA data
GISAID_H1N1_PA_combined <- load_and_combine_data("GISAID_H1N1_PA_2017_FilteredMetadata.csv", "GISAID_H1N1_PA_2017_filtered.fasta")

# Generate H1N1 PA subset
GISAID_H1N1_PA_subset <- generate_subset(Flu_seasons, 8, 8, GISAID_H1N1_PA_combined, 25, 50, 1)

# Combine and remove duplicates for H1N1
GISAID_H1N1_subset <- rbind(GISAID_H1N1_USA_subset, GISAID_H1N1_PA_subset)
GISAID_H1N1_subset <- remove_duplicates(GISAID_H1N1_subset)

# Save the subsetted GISAID sequences
write.csv(GISAID_H1N1_subset,file.path(Phylogeny_path,"GISAID_H1N1_subset.csv"))
