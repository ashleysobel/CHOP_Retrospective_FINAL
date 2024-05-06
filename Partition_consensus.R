# This script partitions consensus sequences and associated mutation data by
# strain and gene segment, saving the results in structured directories and files. It
# also manages file paths, reads sample and mutation data, and filters based on
# specific criteria.

# Prepare environment -----------------------------------------------------
# Clear the workspace to ensure a clean environment for running the script
rm(list = ls())

# Get the current working directory to use as the base for all relative paths
R_path <- getwd()

# Define paths to where data structures are stored/should be written
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
full_fasta_path <- file.path(Compiled_OutPut_path,"Consensus_seqs")


# Load Packages -----------------------------------------------------
library(tidyverse)
library(phylotools)

# Generate .bib citation file for statistical packages used in analysis
# List packages for citation and write citations to a .bib file
package_list <- c("tidyverse","phylotools")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","Partition_Conensus.bib"))



# Set options for run -----
RunDate <-
  c(33122, 50422, 61722, 91222, 22323, 31023, 50123) # These are the run dates of the sample tracker to include

# Study <- "CHOA" # Identify the study of interest, here it is CHOA
# Strain <- "H3N2" # Identify the strain of interest, here is H3N2. This should be repeated for H1N1.

# Define functions ----- 
# Function to check if the fasta file name contains any selected sample names
contains_sample <- function(fasta_name, samples) {
  any(sapply(samples, function(sample) grepl(sample, fasta_name)))
}


# Create a function that partitions the fasta files and mutations by strain and study
Partition_Compiled_Info <- function(Compiled_OutPut_path,Selected_sample_list,Study,Strain,MutationOutput_Compiled,full_fasta_path){
  
  # Create the directory where you will save your fasta output
  partitioned_fasta_path <- file.path(Compiled_OutPut_path,paste0("Consensus-",Study,"_",Strain))
  dir.create(partitioned_fasta_path)
  
  Selected_samples_this <- Selected_sample_list[Selected_sample_list$Strain == Strain & Selected_sample_list$Study == Study ,]
  Selected_samples_this
  
  # List all fasta files in the directory
  fasta_files <- list.files(full_fasta_path, pattern = "\\.fasta$")
  fasta_files
  
  filtered_fasta_files <- fasta_files[sapply(fasta_files, contains_sample, Selected_samples_this$Sample)]
  
  
  # Define the gene segments
  gene_segments <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")
  
  # Iterate over gene segments
  gs <- 1
  for (gs in 1:length(gene_segments)) {
    # Initialize a data.frame to store the sequences for the current gene segment
    
    gene <- gene_segments[gs]
    gene_data <- data.frame(seq.name = character(), seq.text = character(), stringsAsFactors = FALSE)
    
    # Iterate over fasta files
    i <- 1
    for (i in 1:length(filtered_fasta_files)) {
      # Read the fasta file
      fasta_path <- file.path(full_fasta_path, filtered_fasta_files[i])
      seqs <- phylotools::read.fasta(fasta_path)
      
      # Iterate through the sequences and find the matching gene segment
      j <- 1
      for (j in 1:length(seqs$seq.name)) {
        seq_id <- seqs$seq.name[j]
        seq_text <- seqs$seq.text[j]
        seq_id_parts <- strsplit(seq_id, "_")[[1]]
        gene_segment <- seq_id_parts[length(seq_id_parts)]
        
        # Check if the current sequence is for the target gene segment
        if (gene_segment == gene) {
          # Add the sequence to the gene_data data.frame
          gene_data <- rbind(gene_data, data.frame(seq.name = seq_id, seq.text = seq_text, stringsAsFactors = FALSE))
        }
      }
    }
    # Save the current gene segment data in a fasta file
    output_filename <- paste0(Study,"_",Strain,"-",gene, ".fasta")
    output_file_path <- file.path(partitioned_fasta_path,output_filename)
    dat2fasta(gene_data,output_file_path)
  }
  
  Mutations_this <- MutationOutput_Compiled %>%
    filter(subject %in% Selected_samples_this$Sample)
  Mutations_this$X <- NULL
  
  Mutation_name <- paste0("Mutations_",Study,"_",Strain,".csv")
  write.csv(Mutations_this,file.path(Compiled_OutPut_path,Mutation_name))
}

# Load reference data -----

# Load samples selected for further analysis
CHOA_Selected_Samples <- read.csv(file.path(Compiled_OutPut_path,"CHOA_Selected_Samples_Revised.csv"))
CHOA_Selected_Samples$X  <- NULL # Remove unwanted index column

# Load compiled mutation output data
MutationOutput_Compiled <- read.csv(file.path(Compiled_OutPut_path,"Final_MutationOutput_minvar3.csv"))
MutationOutput_Compiled$X <- NULL # Remove unwanted index column

# Generate Data Structures for sequence organization ----- 
# Prepare a list of selected samples with their associated strain for filtering
Selected_sample_list <- data.frame(Sample=CHOA_Selected_Samples$Sample,Strain=CHOA_Selected_Samples$Strain)

# This section is to ensure only the sequences for the current study (CHOA) are
# included. For reference, "CHOA" is the retropsective flu study using sequences
# from 2017-2018, "CHOB" is the pilot year of the prospective flu study
# (2021-2022), "CHOC" is the first true year of the prospective flu study
# (2022-2023)
Selected_sample_list <- Selected_sample_list %>%
  mutate(Study = case_when(
    grepl("CHOA", Sample) ~ "CHOA",
    grepl("CHOB", Sample) ~ "CHOB",
    grepl("CHOC", Sample) ~ "CHOC",
    TRUE ~ NA_character_
  ))

# Run code for consensus generation ---
Partition_Compiled_Info(Compiled_OutPut_path,Selected_sample_list,"CHOA","H3N2",MutationOutput_Compiled,full_fasta_path)
Partition_Compiled_Info(Compiled_OutPut_path,Selected_sample_list,"CHOA","H1N1",MutationOutput_Compiled,full_fasta_path)

