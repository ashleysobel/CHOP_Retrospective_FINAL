# This script clears the R environment, sets file paths, loads necessary
# libraries, and processes mutation data to analyze and update entries based on
# reference sequences, integrating functions to identify matching mutations
# across samples, and generating a final updated mutation table. Specifically,
# it includes all mutations present at 1%, adds the codon position, incorporates
# the number of samples the mutation was found in and which samples those are
# Author: Ashley Sobel Leonard Date: 5/1/2024


# Clear the environment to get rid of old files and get path -------------------
rm(list = ls())
R_path <- getwd()


# Define paths to where data structures are stored/should be written
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
Mutation_analysis_path <- file.path(Compiled_OutPut_path,"Mutation_analysis")
dir.create(Mutation_analysis_path)


# Load Packages -----------------------------------------------------
library(tidyverse)

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("tidyverse")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","GenerateMutation_Table.citations.bib"))

# Define functions -------

# Function to find matching mutations in a dataframe based on certain criteria
find_matching_mutations <- function(this_mut, df) {
  # Convert relevant columns to the appropriate data types for accurate comparison
  current_strain <- as.character(this_mut["Strain"])
  current_segment_num <- as.numeric(this_mut["segment_num"])
  current_position <- as.numeric(this_mut["position"])
  current_position
  
  # Use dplyr to filter the dataframe
  # Filter rows where Strain, segment_num, and position match the current mutation
  matching_rows <- df %>%
    filter(Strain == current_strain, 
           segment_num == current_segment_num, 
           position == current_position)

  # Create a sorted unique list of subject IDs from the matching rows
  # This removes duplicates and sorts the IDs for consistency
  unique_subjects <- sort(unique(matching_rows$Sample_ID_mut))

  # Return a list containing:
  # N_samples: the count of unique samples where the mutation was found
  # SamplesIdentified: a string concatenating the sample IDs where the mutation was found, separated by commas
  list(N_samples = length(unique_subjects), 
       SamplesIdentified = paste(unique_subjects, collapse = ", "))
}


# Load data --------
Mutation_table <- read.csv(file.path(Compiled_OutPut_path,"FINAL_MutationOutput_minvar3.csv"))
Mutation_table$X <- NULL

# Replace all entries for location of "NA" with "NA." to distinguish from empty
Mutation_table$location[is.na(Mutation_table$location)] <- "NA."

# Enter the csv file containing the references you want to use. For the CHOP
# 2017-2018 retrospective study, we want to use H3N2 2022 and H1N1 2019.
H3N2_Ref <- read.csv("Reference_files/H3N2_Ref_A_2016.csv")
H1N1_Ref <- read.csv("Reference_files/H1N1_Ref_A_2015.csv")

# Update mutation table ------
nMut <- nrow(Mutation_table)
nm <- 1
Updated_Mutation_table <- data.frame()
for (nm in 1:nMut){
  # Extract mutation for this 
  Mutation_tmp <- Mutation_table[nm,]; Mutation_tmp
  nt_position <- Mutation_tmp$position
  seg <- Mutation_tmp$segment_num
  
  # Assign reference based on strain 
  if (Mutation_tmp$Strain == "H3N2"){
    Ref <- H3N2_Ref
  } else if (Mutation_tmp$Strain == "H1N1"){
    Ref <- H1N1_Ref
  } else {
    stop(print("Error: the reference is incorrect"))
  }
  
  # Select reference segment
  Ref_seg <- Ref[seg,]
  
  # Now we determine the codon position of the mutation. We start by checking if
  # the mutation occurs in a non-coding region. If so, the coding position is
  # irrelevant
  if (Mutation_tmp$location == "NonCoding"){
    Mutation_tmp$CodonPosition <- "NonCoding"
  } else { 
    # Mutation occurs in the coding region, now we determine the coding position
    # by comparing the nucleotide position relative to the start of the ORF in
    # the reference segment
    Mutation_tmp$CodonPosition <- (Mutation_tmp$position - Ref_seg$Start) %% 3 + 1
  }
  
  # We repeat this for the second coding region (if present)
  if (Mutation_tmp$location2 == "NaN"){
    # There isn't a 2nd ORF
    Mutation_tmp$CodonPosition2 <- NaN
  } else if (Mutation_tmp$location2 == "NonCoding"){
    # The mutation is in a non-coding region of the 2nd ORF
    Mutation_tmp$CodonPosition2 <- "NonCoding"
  } else {
    Mutation_tmp$CodonPosition2 <- (Mutation_tmp$position - Ref_seg$Start) %% 3 + 1
  }
  
  # Next we generate a column to reflect if the mutation was called by bbmap in
  # the replicate. This information is aleady contained in the Presence column,
  # but we add an additional column to avoid confusion
  Mutation_tmp$called_by_bbmap <- Mutation_tmp$Presence
  
  # We can also remove the Presence_corrected to avoid confusion (this is where the frequency was adjusted to that determined from the pileup)
  Mutation_tmp$Presence_corrected <- NULL
    
  # Apply the find_matching_mutations function to the current mutation
  matching_muts <- find_matching_mutations(Mutation_tmp, Mutation_table)
  
  # Unlist the results from the function and assign them to Mutation_tmp
  Mutation_tmp$N_samples <- matching_muts$N_samples
  Mutation_tmp$SamplesIdentified <- matching_muts$SamplesIdentified

    # Append the processed Mutation_tmp to the Updated_Mutation_table
  Updated_Mutation_table <- rbind(Updated_Mutation_table, Mutation_tmp)
}

write.csv(x = Updated_Mutation_table,file = file.path(Mutation_analysis_path,"20624-Updated_MutationTable.csv"))

