# This script identifies and visualizes variants at mismatch sites in influenza
# sequences from mixed infections, focusing on their frequency and distribution
# across specific gene segments. It extracts, aligns, and compares sequences
# from individual patients, computes variant frequencies, and produces detailed
# graphical outputs to illustrate the locations and prevalence of these
# mismatches. This analysis distinguishes itself by integrating sequence
# alignment, variant filtering, and visual representation to provide insights
# into viral diversity within hosts.

# Prepare environment ----------------------------------------------------- 
# Clear the workspace to ensure a clean environment for running the script
rm(list = ls())

# Get the current working directory to use as the base for all relative paths
R_path <- getwd()

# Define paths to where data structures are stored/should be written
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
Masked_Fasta_path <- file.path(R_path,"PostProcessing_OutPut","MaskedConsensus_CHOP")
# Mutation_analysis_path <- file.path(Compiled_OutPut_path,"Mutation_analysis")

# Generate a new directory to hold the output of the mixed infection plot 
MixedInfection_OutPut_path <- file.path(Compiled_OutPut_path,"MixedInfection")
dir.create(MixedInfection_OutPut_path)

# Load packages -----------------
library(tidyverse)
library(tidysq)
library(Biostrings)
library(gridExtra)

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("tidyverse","tidysq","Biostrings", "gridExtra")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","MixedInfections_citations.bib"))


# Define functions --------

# Converts alignment data into a structured data frame for analysis
convert_alignment_to_df <- function(pattern, subject, mismatch_positions) {
  # Assuming pattern is your sequence and subject is the reference
  
  seq_length <- max(nchar(pattern), nchar(subject))
  positions <- 1:seq_length
  
  # Split sequences into individual characters
  pattern_chars <- strsplit(pattern, "")[[1]]
  subject_chars <- strsplit(subject, "")[[1]]
  
  # Initialize a data frame
  alignment_df <- data.frame(
    Position = positions,
    Base = pattern_chars,
    ReferenceBase = subject_chars,
    Match = !positions %in% mismatch_positions,
    stringsAsFactors = FALSE
  )
  
  return(alignment_df)
}

# Filters variants based on frequency and coverage thresholds
Filter_variants <- function(CHOA_HiVar_Expanded_intrahost, min_var_final, min_cov_final,ID_this,subject1, subject2) {
  
  CHOA_Muts <- CHOA_HiVar_Expanded_intrahost[CHOA_HiVar_Expanded_intrahost$Sample_ID_mut == ID_this,]
  
  CHOA_Muts_1_filtered <- CHOA_Muts %>%
    filter(CHOA_Muts$subject == subject1, CHOA_Muts$var_freq > min_var_final, CHOA_Muts$depth > min_cov_final)
  dim(CHOA_Muts_1_filtered)
  
  # Filter CHOA_Muts_101R based on conditions
  CHOA_Muts_2_filtered <-  CHOA_Muts %>%
    filter(CHOA_Muts$subject == subject2, CHOA_Muts$var_freq > min_var_final, CHOA_Muts$depth > min_cov_final)
  head(CHOA_Muts_1_filtered)
  dim(CHOA_Muts_2_filtered)
  
  # Extract entries present in both dataframes
  common_entries <- inner_join(CHOA_Muts_1_filtered, CHOA_Muts_2_filtered, by = c("segment_num", "position", "var_nt"))
  
  # Calculate mean of var_freq.x and var_freq.y and create mean_freq column
  common_entries <- common_entries %>%
    mutate(mean_freq = (var_freq.x + var_freq.y) / 2)
  
  # Return common entries
  return(common_entries)
}

# Splits variant data into those occurring at mismatch sites versus others
split_variants <- function(merged_data) {
  library(dplyr)
  
  # Filter for variants at mismatches (VariantFrequency != NA AND Match = FALSE)
  Variants_at_mismatches <- merged_data |>
    filter(!is.na(VariantFrequency) & Match == FALSE)
  
  # Filter for variants not at mismatches (VariantFrequency != NA AND Match = TRUE)
  Variants_not_mismatches <- merged_data |>
    filter(!is.na(VariantFrequency) & Match == TRUE)
  
  # Return a list containing both dataframes
  return(list(Variants_at_mismatches = Variants_at_mismatches,
              Variants_not_mismatches = Variants_not_mismatches))
}

# Generates a plot visualizing the frequency and distribution of variants along a gene segment
create_plot <- function(data, subject, gene_segment, min_nt, max_nt) {
  # Calculate the maximum position rounded to the nearest hundred
  max_pos_rounded <- floor(max_nt / 100) * 100

  if (min_nt == 1){
    # Generate the vector for x-ticks from 100 to max_value, and concatenate with 1
    xticks <- c(1, seq(100, max_pos_rounded, by = 100))
  } else {
    xticks <- seq(min_nt, max_pos_rounded, by = 100)
  }
  
  # Plot the location of intrahost variants relative to mismatch sites
  p <- ggplot(data) +
    geom_bar(data = subset(data, !is.na(VariantFrequency)), 
             aes(x = Position, y = VariantFrequency, width = tile_width), 
             stat = "identity", fill = "blue") +
    geom_tile(aes(x = Position, y = -tile_height/2, height = tile_height, width = tile_width, fill = Match)) +
    geom_hline(yintercept = 0.03, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("TRUE" = "lightgrey", "FALSE" = "red")) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(angle = 90, vjust = 0.5, size = 12, face = "plain"),
      axis.text.y = element_text(size = 12, face = "plain"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, face = "plain"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    labs(
      title = paste(subject, gene_segment, "Segment Variant Frequencies"),
      x = "Nucleotide Position",
      y = "Variant Frequency"
    ) +
    coord_cartesian(clip = 'off', ylim = c(-0.01, 0.1), xlim = c(min_nt, max_nt)) +
    scale_x_continuous(expand = c(0, 0), limits = c(min_nt, max_nt), breaks = xticks) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.05))
  return(p)
}

# Function to merge alignment and variant data 
Merge_AlignmentsVariants <- function(alignment_data_101_HA, CHOA101_HA_var) {
  merged_data_101_HA <- alignment_data_101_HA %>%
    left_join(CHOA101_HA_var, by = c("Position" = "position"))
  
  merged_data_101_HA <- merged_data_101_HA %>%
    mutate(VariantFrequency = mean_freq) %>%
    select(Position, Base, ReferenceBase, Match, VariantFrequency, effect.x, aa_pos.x, aa_ref.x,aa_sub.x)
  colnames(merged_data_101_HA) <- c("Position","Base","ReferenceBase","Match","VariantFrequency","Effect","aa_pos","aa_ref","aa_sub")
  
  return(merged_data_101_HA)
}

# Function to identify mismatches
identify_mismatches <- function(seq1, seq2) {
  mismatches <- seq1 != seq2
  mismatches[seq1 == "-"] <- FALSE # Ignore gaps
  mismatches[seq2 == "-"] <- FALSE # Ignore gaps
  mismatches[seq1 == "N"] <- FALSE # Ignore masked nts
  mismatches[seq2 == "N"] <- FALSE # Ignore masked nts
  return(which(mismatches))
}


# Set Run Options ------
# Set selected replicates for CHOA-101 and CHOA-117
CHOA101 <- "CHOA-101_D2V1"
CHOA101R <- "CHOA-101_D3V1"
CHOA117 <- "CHOA-117_D1V1"
CHOA117R <- "CHOA-117_D2V1"
sample_list <- c(CHOA101,CHOA101R,CHOA117,CHOA117R)

# Set chosen parameters for variant filtering 
min_var_final <- 0.01
min_cov_final <- 100


# Load data -------------
# Load compiled mutations
MutationOutput_Compiled <- read.csv(file.path(Compiled_OutPut_path,"MutationOutput_Compiled.csv"))


# Fix inter/intrahost variants
# -------------------------------------------------------- As in previous code,
# there are some instances where an intrahost variant is called as a consensus
# level variant and vice versa due to different reference sequencess (i.e.
# stated reference for consensus level variation vs consensus sequence for
# intrahost level) 

# Initialize an empty DataFrame to collect mutation data from specified samples
CHOA_HiVar_Expanded_MutationOutput <- data.frame()

# Iterate over each sample to extract and compile mutation data
nv <- 1
for (nv in 1:length(sample_list)) {
  # Print the number of remaining iterations to the console
  print(length(sample_list) - nv)
  
  # Get the current sample from the list
  sample_tmp <- sample_list[nv]
  
  # Check and proceed if the current sample is valid
  if (is.nan(sample_tmp) == FALSE) {
    # Extract the mutations present in that sample from the compiled mutation output data frame
    sample_muts <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject == sample_tmp, ]
    
    # Add the extracted mutations to the mutation output data frame for selected samples
    CHOA_HiVar_Expanded_MutationOutput <- rbind(CHOA_HiVar_Expanded_MutationOutput, sample_muts)
    
    # Remove the temporary data frame to free up memory
    rm(sample_muts)
    
    # Extract the subjects from the mutation output data frame and store them in a new variable
    subject_mut <- CHOA_HiVar_Expanded_MutationOutput$subject
    
    # Extract various parts of the subject field for further analysis
    
    # Extract the sample ID (first 8 characters of the subject field)
    Sample_ID_mut <- stringr::str_sub(CHOA_HiVar_Expanded_MutationOutput$subject, 1, 8)
    
    # Extract the replicate number (characters 10 and 11 of the subject field)
    Replicates_mut <- stringr::str_sub(CHOA_HiVar_Expanded_MutationOutput$subject, 10, 11)
    
    # Extract the version number (characters 12 and 13 of the subject field)
    Versions_mut <- stringr::str_sub(CHOA_HiVar_Expanded_MutationOutput$subject, 12, 13)
    
    # Combine the extracted fields with the original mutation output data frame
    CHOA_Expanded_mutations <-
      cbind(Sample_ID_mut,
            Replicates_mut,
            Versions_mut,
            CHOA_HiVar_Expanded_MutationOutput)
    
  }
}

# Extract mutation data specifically flagged as 'Intra' for intra-host analysis
CHOA_HiVar_Expanded_intrahost <-
  CHOA_Expanded_mutations[CHOA_Expanded_mutations$variant_level == "Intra", ]

# Remove the row names from the intra-host data frame
rownames(CHOA_HiVar_Expanded_intrahost) <- NULL
write.csv(CHOA_HiVar_Expanded_intrahost,file = file.path(MixedInfection_OutPut_path,"CHOA_HiVar_Expanded_IntrahostCompiled.csv"))

# Extract and filter mutation data for selected samples and gene segments
CHOA_101_filtered_muts <- Filter_variants(CHOA_HiVar_Expanded_intrahost,min_var_final,min_cov_final,ID_this = "CHOA-101",subject1 = CHOA101,subject2 = CHOA101R)
CHOA_117_filtered_muts <- Filter_variants(CHOA_HiVar_Expanded_intrahost,min_var_final,min_cov_final,ID_this = "CHOA-117",subject1 = CHOA117,subject2 = CHOA117R)

# Load sequence data for analysis and visualization
CHOA101_HA_var <- CHOA_101_filtered_muts[CHOA_101_filtered_muts$segment_num == 4,]
CHOA117_HA_var <- CHOA_117_filtered_muts[CHOA_117_filtered_muts$segment_num == 4,]

# Load fasta files for CHOA-101 and CHOA-117
CHOA_101_id <- paste0(CHOA101,"_Masked.fasta")
CHOA_117_id <- paste0(CHOA117,"_Masked.fasta")
  
# Read in CHOA-101 consensus sequence and extact the HA gene segment
CHOA_101_consensus <- list.files(path=Masked_Fasta_path,pattern = CHOA_101_id,recursive = TRUE,full.names = TRUE)
CHOA_101_fasta <- tidysq::read_fasta(CHOA_101_consensus,alphabet="dna")
CHOA_101_HA <- CHOA_101_fasta$sq[4]

# Repeat for CHOA-117
CHOA_117_consensus <- list.files(path=Masked_Fasta_path,pattern = CHOA_117_id,recursive = TRUE,full.names = TRUE)
CHOA_117_fasta <- tidysq::read_fasta(CHOA_117_consensus,alphabet="dna")
CHOA_117_HA <- CHOA_117_fasta$sq[4]

# Load fasta file for A/Kansas/14/2017 (3C.3a1 reference)
HA_3c3a1 <- tidysq::read_fasta(file.path(R_path,"Reference_files","A_Kansas_14_2017_HA.fasta"),alphabet = "dna")
HA_3c3a1$name <- "A/Kansas/14/2017"

# Assuming HA_3c3a1 and CHOA_101_HA are DNAStringSet objects or similar
HA_3c3a1_char <- as.character(HA_3c3a1$sq)
CHOA_101_HA_char <- as.character(CHOA_101_HA)
CHOA_117_HA_char <- as.character(CHOA_117_HA)

# Perform pairwise alignment to identify and map mismatches with the 3c3a1 reference strain
CHOA101_HA_alignment <- pairwiseAlignment(CHOA_101_HA_char, HA_3c3a1_char)
CHOA117_HA_alignment <- pairwiseAlignment(CHOA_117_HA_char, HA_3c3a1_char)


# Extract aligned sequences
aligned_seq1_101_HA <- as.character(CHOA101_HA_alignment@pattern)
aligned_seq2_101_HA <- as.character(CHOA101_HA_alignment@subject)
aligned_seq1_117_HA <- as.character(CHOA117_HA_alignment@pattern)
aligned_seq2_117_HA <- as.character(CHOA117_HA_alignment@subject)

# Identify mismatches in the aligned sequences to filter relevant variants
mismatch_positions_101_HA <- identify_mismatches(strsplit(aligned_seq1_101_HA, "")[[1]], strsplit(aligned_seq2_101_HA, "")[[1]])
mismatch_positions_117_HA <- identify_mismatches(strsplit(aligned_seq1_117_HA, "")[[1]], strsplit(aligned_seq2_117_HA, "")[[1]])

# Convert aligned sequences into a data frame for further analysis
alignment_data_101_HA <- convert_alignment_to_df(aligned_seq1_101_HA, aligned_seq2_101_HA, mismatch_positions_101_HA)
alignment_data_117_HA <- convert_alignment_to_df(aligned_seq1_117_HA, aligned_seq2_117_HA, mismatch_positions_117_HA)

# Merge alignment data with variant information to prepare for plotting
merged_data_101_HA <- Merge_AlignmentsVariants(alignment_data_101_HA,CHOA101_HA_var)
merged_data_117_HA <- Merge_AlignmentsVariants(alignment_data_117_HA,CHOA117_HA_var)

# Determine the proportion of variants that occur at mismatch sites --------

# Analyze the proportion of variants that occur at mismatch sites and report them for CHOA-101 and CHOA-117
HA_101_result <- split_variants(merged_data_101_HA)
Variants_at_mismatches_HA_101 <- HA_101_result$Variants_at_mismatches
Variants_not_mismatches_HA_101 <- HA_101_result$Variants_not_mismatches
print(paste("The proportion of CHOA-101 HA variants occurring at mismatch sites =",
            nrow(Variants_at_mismatches_HA_101) / 
              sum(nrow(Variants_at_mismatches_HA_101), nrow(Variants_not_mismatches_HA_101))))

HA_117_result <- split_variants(merged_data_117_HA)
Variants_at_mismatches_HA_117 <- HA_117_result$Variants_at_mismatches
Variants_not_mismatches_HA_117 <- HA_117_result$Variants_not_mismatches
nrow(Variants_at_mismatches_HA_117)/sum(nrow(Variants_at_mismatches_HA_117),nrow(Variants_not_mismatches_HA_117))
print(paste("The proportion of CHOA-117 HA variants occurring at mismatch sites =",
            nrow(Variants_at_mismatches_HA_117) / 
              sum(nrow(Variants_at_mismatches_HA_117), nrow(Variants_not_mismatches_HA_117))))


# Generate mixed infection plots ------ ---
# Set dimensions for the square-like appearance of the tiles
tile_height <- 0.0025  # Smaller height for a square appearance
tile_width <- 1     # Width close to 1 for square-like appearance

# Create plots illustrating the distribution and frequency of variants at mismatch sites
CHOA101_HA_plot <- create_plot(data = merged_data_101_HA, subject = "CHOA-101", gene_segment = "HA", min = 1, max = max(merged_data_101_HA$Position))
print(CHOA101_HA_plot)

CHOA117_HA_plot <- create_plot(data = merged_data_117_HA, subject = "CHOA-117", gene_segment = "HA",min = 300, max = 700)
print(CHOA117_HA_plot)

# Arrange and save the plots in EPS format for high-quality publication outputs
arranged_plots_HA <- grid.arrange(CHOA101_HA_plot, CHOA117_HA_plot, ncol = 1)
file_name_HA <- file.path(MixedInfection_OutPut_path, "MixedInfection_HA_plot.eps")
ggsave(file = file_name_HA, arranged_plots_HA, device = "eps",width = 10, height = 6)
