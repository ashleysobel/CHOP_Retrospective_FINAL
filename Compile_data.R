# Title: Pipeline for Processing NGS Data for Influenza Analysis
# Description: This script processes NGS data to identify influenza variants,
# perform quality control, and generate necessary output directories and citation files.
# The script ensures reproducibility and efficient data management by handling file operations,
# setting QC parameters, and characterizing mutations based on reference sequences.
# Inputs: Path to the NGS data, run dates for sample processing.
# Outputs: Processed data files organized in specified directories, citation file for referenced R packages.
# Author: Ashley Sobel Leonard
# Date: 5/1/2024


# Prepare environment -----------------------------------------------------
# Clear the workspace to ensure a clean environment for running the script
rm(list = ls())
R_path <- getwd() # Get the current working directory to use as the base for all relative paths


# Set options for run ------------
RunDate <-
  c(33122, 50422, 61722, 91222, 22323, 31023, 50123) # These are the run dates of the sample tracker to include

# Apply filters to select mutations based on frequency and coverage criteria (we
# set this later in the code, but you could do it here instead, make sure to
# delete/comment out the later instances though)
# min_var_final <- 0.03
# min_var_replicate_final <- 0.03
# min_cov_final <- 100

# Load reference data ----- 
# Specify segment names for the influenza analysis, corresponding to different parts of the virus genome
segment_name <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")

# Load reference files for H3N2 and H1N1 influenza strains 
H3N2_Ref <- read.csv("Reference_files/H3N2_Ref_A_2016.csv")
H1N1_Ref <- read.csv("Reference_files/H1N1_Ref_A_2015.csv")

# Create directories for storing output - ensure directories exist to avoid errors in file writing
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")
dir.create(Compiled_OutPut_path)
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
dir.create(Figure_OutPut_path)

# Load Packages -----------------------------------------------------

# Load required packages for data manipulation, sequence analysis, and plotting
library(tidyverse)
library(phylotools) # 'phylotools' for reading FASTA files,
library(gridExtra) #  arranging grid graphical objects
library(seqinr) # 'seqinr' for biological sequences retrieval and analysis
library(ShortRead) # Also loads biostrings
library(data.table) # fast aggregation of large data
library(cowplot)

# Generate .bib citation file for statistical packages used in analysis
# List packages for citation and write citations to a .bib file
package_list <- c("tidyverse", "phylotools","gridExtra","seqinr","ShortRead","data.table","cowplot")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","Compile_data_citations.bib"))

# Define Functions -----------------------------------------------------

#  Writes sequence data to a FASTA format file.
writeFasta <- function(data, filename) {
  fastaLines = c() # Initialize an empty vector to hold FASTA formatted strings
  for (rowNum in 1:nrow(data)) { # Loop through each row in the data frame
    # Create FASTA header and sequence strings, then add to fastaLines
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum, "name"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, "seq"]))
  }
  fileConn <- file(filename) # Open a file connection for writing
  writeLines(fastaLines, fileConn) # Write the FASTA formatted data to file
  close(fileConn) # Close the file connection
}

# Compares two strings and returns differences, excluding specified characters.
GetStringDiff <- 
  function(a,
           b,
           exclude = c("-", "?"),
           ignore.case = TRUE,
           show.excluded = FALSE) {
    if (nchar(a) != nchar(b))
      stop("Lengths of input strings differ. Please check your input.") # Ensure strings are the same length
    if (ignore.case)
    {
      a <- toupper(a) # Convert strings to uppercase if ignoring case
      b <- toupper(b)
    }
    split_seqs <- strsplit(c(a, b), split = "") # Split strings into individual characters
    only.diff <- (split_seqs[[1]] != split_seqs[[2]]) # Identify positions where characters differ
    only.diff[(split_seqs[[1]] %in% exclude) | 
                (split_seqs[[2]] %in% exclude)] <- NA # Exclude specified characters
    diff.info <- data.frame(which(is.na(only.diff) | only.diff),
                            split_seqs[[1]][only.diff], split_seqs[[2]][only.diff]) # Create a data frame of differences
    names(diff.info) <- c("position", "poly.seq.a", "poly.seq.b") # Set column names for clarity
    if (!show.excluded) 
      diff.info <- na.omit(diff.info) # Optionally remove rows with NAs if exclusions are not to be shown
    diff.info # Return the difference information
  }

# Checks for the presence of a mutation in the replicate samples and assesses read depth and frequency from pileup data.
Check_ReplicateMut <- function(CHOA_targets, mut_tmp, CHOA_Selected_Samples, pileup_path) {
  # Filter out mutations from the same subject to identify mutations in the replicate subject
  paired_mutation <-
    CHOA_targets[CHOA_targets$subject != mut_tmp$subject, ]
  
  # Identify the replicate subject based on the unique ID but different sample identifier
  replicate_subject <-
    CHOA_Selected_Samples[CHOA_Selected_Samples$ID == mut_tmp$Sample_ID_mut &
                            CHOA_Selected_Samples$Sample != mut_tmp$subject, ]$Sample
  
  # Narrow down mutations to those that match the replicate subject's identifier
  paired_mutation <-
    paired_mutation[paired_mutation$subject == replicate_subject, ]
  
  # Further filter to match both the genomic segment and exact position of the mutation
  paired_mutation <-
    paired_mutation[paired_mutation$segment_num == mut_tmp$segment_num, ]
  
  paired_mutation <-
    paired_mutation[paired_mutation$position == mut_tmp$position, ]
  
  # Check if there is more than one mutation that matches the variant nucleotide
  if (nrow(paired_mutation) > 1) {
    # Filter the paired mutation data frame to only include rows where the variant nucleotide matches the current mutation's variant nucleotide
    paired_mutation <-
      paired_mutation[paired_mutation$var_nt == mut_tmp$var_nt, ]
  }
  
  # If exactly one matching mutation is found, confirm its presence as a replicate
  if (nrow(paired_mutation) == 1) {
    # mutation is present as replicate, huzzah!, the mutation is present as a replicate
    ReplicateSubject <- replicate_subject
    Presence <- "Y"
    RepFreq <- paired_mutation$var_freq
    RepDepth <- paired_mutation$depth
  } else {
    # Check if the replicate exists but the mutation is not detected
    IDtmp <- mut_tmp$Sample_ID_mut
    test <-
      CHOA_Selected_Samples[CHOA_Selected_Samples$ID == IDtmp, ]
    if (is.na(test$ID[2] == TRUE)) {
      # There is no replicate mutation because there's no replicate, populate the subsequent entries to reflect that 
      ReplicateSubject <- "No Rep"
      Presence <- "No Rep"
      RepFreq <- NaN
      RepDepth <- NaN
    } else {
      # There is a replicate but the mutation doesn't seem to be called in the
      # mutaiton file. We see if there is a mutation present based on nucleotide
      # frequency, but not called, by assessing the pileup file
      Presence <- "N"
      replicate <-
        test[test$Sample != mut_tmp$subject, ]
      
      # identify and load the corresponding pileup file
      ReplicateSubject <- replicate_subject
      target_file <-
        list.files(path = pileup_path, pattern = replicate$Sample)
      pileup_tmp <- read.csv(file.path(pileup_path, target_file))
      pileup_tmp <- data.table::data.table(pileup_tmp)
      
      # Extract the information from the pileup file about the current mutation location 
      Rep_reads <-
        pileup_tmp[SegNum == mut_tmp$segment_num & pos == mut_tmp$position]
      RepDepth <- sum(Rep_reads$count)
      Var_reads <-
        Rep_reads[nucleotide == mut_tmp$var_nt]
      if (nrow(Var_reads) == 0) {
        RepFreq <- 0
      } else {
        RepFreq <- Var_reads$count / RepDepth
      }
    }
  }
  # Save the information about the replicate (whether the mutation was called,
  # present based on frequency but not called, the frequency and the depth)
  replicate_mut_info_this <-
    data.frame(ReplicateSubject=ReplicateSubject, Presence=Presence, RepFreq=RepFreq, RepDepth=RepDepth)
  return(replicate_mut_info_this)
}

# Generates a consensus sequence from multiple DNA sequences by identifying the most common nucleotide at each position.
get_majority_consensus <- function(seq_text) {
  # Create a matrix with the counts of each nucleotide for each position
  counts <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(seq_text), baseOnly = TRUE)
  # Extract the majority nucleotides using the max count for each column
  majority_nucleotides <- apply(counts, 2, function(x) names(which.max(x)))
  majority_nucleotides <- gsub("other", "N", majority_nucleotides)
  # Combine the majority nucleotides into a single consensus sequence
  consensus_seq <- paste(majority_nucleotides, collapse = "")
  return(consensus_seq)
}

# Constructs consensus sequences for each influenza segment from high-quality samples and tracks mismatches against the majority sequence.
Generate_consensus_sequences <- function(uID, SampleTracker_Compiled, R_path, segment_name) {
  # Generate dataframes to hold information
  mismatches <-
    data.frame(
      Sample = character(),
      RunDate = character(),
      SegNum = double(),
      Position = double(),
      NT_1 = character(),
      NT_2 = character()
    )
  PassingSample_Matches <- data.frame()
  
  for (ns in 1:length(uID)) {
    IDtmp <- uID[ns]
    # Retrieves all valid sample versions based on the unique sample identifier.
    passing_samples <-
      SampleTracker_Compiled[grep(IDtmp, SampleTracker_Compiled$Sample), ]
    CHOA_PB2 <- data.frame(seq.name = character(), seq.text = character())
    CHOA_PB1 <- data.frame(seq.name = character(), seq.text = character())
    CHOA_PA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_HA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NP <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NA <- data.frame(seq.name = character(), seq.text = character())
    CHOA_MP <- data.frame(seq.name = character(), seq.text = character())
    CHOA_NS <- data.frame(seq.name = character(), seq.text = character())
    
    # Processes each sample to generate consensus sequences by evaluating all masked FASTA sequences.
    ps <- 1
    for (ps in 1:nrow(passing_samples)) {
      # Get name of masked consensus sequence
      passing_tmp <- passing_samples[ps, ]
      masked_path <-
        paste0(R_path,
               "/PostProcessing_Output/MaskedConsensus_CHOP/",
               passing_tmp$SeqDate,
               "_MaskedConsensus")
      target_fasta <-
        list.files(path = masked_path, pattern = passing_tmp$Sample)
      tmp_fasta <-
        phylotools::read.fasta(file.path(masked_path, target_fasta))
      
      # Aggregates sequence data for each segment to later determine the majority consensus sequence.
      CHOA_PB2 <- rbind(CHOA_PB2, tmp_fasta[1, ])
      CHOA_PB1 <- rbind(CHOA_PB1, tmp_fasta[2, ])
      CHOA_PA <- rbind(CHOA_PA, tmp_fasta[3, ])
      CHOA_HA <- rbind(CHOA_HA, tmp_fasta[4, ])
      CHOA_NP <- rbind(CHOA_NP, tmp_fasta[5, ])
      CHOA_NA <- rbind(CHOA_NA, tmp_fasta[6, ])
      CHOA_MP <- rbind(CHOA_MP, tmp_fasta[7, ])
      CHOA_NS <- rbind(CHOA_NS, tmp_fasta[8, ])
    }
    # Compute the consensus sequence for each influenza segment.
    PB2_consensus <-
      get_majority_consensus(CHOA_PB2$seq.text)
    PB1_consensus <-
      get_majority_consensus(CHOA_PB1$seq.text)
    PA_consensus <-
      get_majority_consensus(CHOA_PA$seq.text)
    HA_consensus <-
      get_majority_consensus(CHOA_HA$seq.text)
    NP_consensus <-
      get_majority_consensus(CHOA_NP$seq.text)
    NA_consensus <-
      get_majority_consensus(CHOA_NA$seq.text)
    MP_consensus <-
      get_majority_consensus(CHOA_MP$seq.text)
    NS_consensus <-
      get_majority_consensus(CHOA_NS$seq.text)
    
    # Prepares new FASTA files for each consensus sequence.
    new_name <- paste0(IDtmp, "_consensus.fasta")
    new_frame <- data.frame(seq.name = character(), seq.text = character())
    seq_names <-
      c(
        paste0(IDtmp, "_PB2"),
        paste0(IDtmp, "_PB1"),
        paste0(IDtmp, "_PA"),
        paste0(IDtmp, "_HA"),
        paste0(IDtmp, "_NP"),
        paste0(IDtmp, "_NA"),
        paste0(IDtmp, "_MP"),
        paste0(IDtmp, "_NS")
      )
    seq_names
    seqs <-
      rbind(
        PB2_consensus,
        PB1_consensus,
        PA_consensus,
        HA_consensus,
        NP_consensus,
        NA_consensus,
        MP_consensus,
        NS_consensus
      )
    rownames(seqs) <- NULL
    new_frame <- data.frame(seq.name <- seq_names, seq.text <- seqs)
    colnames(new_frame) <- c("seq.name", "seq.text")
    
    new_fasta_df <-
      dplyr::data_frame(name = new_frame$seq.name, seq = new_frame$seq.text) #This is weird formatting, but it's the only way it works
    consensus_seq_path <- file.path(R_path,"PostProcessing_OutPut/CompiledOutPut/Consensus_seqs",new_name)
    writeFasta(data = new_fasta_df, filename = consensus_seq_path) # Writes the newly generated consensus sequences to FASTA files.
    
    # Initialize a matrix to store information about segment mismatches.
    segment_match <-
      matrix(data = NaN,
             nrow = nrow(passing_samples),
             ncol = 10)
    colnames(segment_match) <- c(segment_name, 'All', "nMismatch")
    segment_match <- as.data.frame(segment_match)
    segment_match$nMismatch <- as.numeric(segment_match$nMismatch)
    
    # Compares each segment of the sample to the newly created consensus sequences
    pps <- 1
    for (pps in 1:nrow(passing_samples)) {
      # Get name of masked consensus sequence
      passing_tmp <- passing_samples[pps, ]
      passing_tmp
      Sample <- passing_tmp$Sample
      masked_path <-
        paste0(R_path,
               "/PostProcessing_Output/MaskedConsensus_CHOP/",
               passing_tmp$SeqDate,
               "_MaskedConsensus")
      target_fasta <-
        list.files(path = masked_path, pattern = passing_tmp$Sample)
      tmp_fasta <-
        phylotools::read.fasta(file.path(masked_path, target_fasta))
      # Iterate over each gene segment to identify places where the nucleotides differ
      for (seg in 1:8) {
        if (nchar(tmp_fasta$seq.text[seg]) == nchar(new_fasta_df$seq[seg])) {
          setDiff <-
            GetStringDiff(a = tmp_fasta$seq.text[seg],
                          b = new_fasta_df$seq[seg],
                          exclude = "N")
          if (nrow(setDiff) == 0) {
            segment_match[pps, seg] <- "Y"
          } else {
            segment_match[pps, seg] <- "N"
            Position <- setDiff[1]
            
            NT_1 <- setDiff[2]
            NT_2 <- setDiff[3]
            # Add new mismatches to tabke of mismatches
            mismatches_this <-
              cbind(Sample, seg, Position, NT_1, NT_2)
            mismatches_this
            mismatches <- rbind(mismatches, mismatches_this) # Records mismatches for further analysis. Samples that are highly mismatched may reflect one poor sequence that can cause exclusion of other, good sequences. 
          }
        } else{
          segment_match[pps, seg] <- "N"
        }
      }
      
      mismatches_this <- NULL
      # We add another column that determines if any of the segments are mismatched
      if ('N' %in% segment_match[pps, ]) {
        segment_match[pps, 9] = "N"
      } else {
        segment_match[pps, 9] = 'Y'
      }
      # We record the number of mismatches
      nMismatch <- sum(mismatches$Sample == Sample)
      segment_match[pps, 10] <- nMismatch
      nMismatch <- NULL
    }
    # We aggregate the information about whether the samples match
    PassingSample_Matches_tmp <-
      cbind(passing_samples$Sample, segment_match)
    PassingSample_Matches <-
      rbind(PassingSample_Matches, PassingSample_Matches_tmp)
    PassingSample_Matches_tmp
    
  }
  return(PassingSample_Matches)
}

# Evaluates the reproducibility of mutations between original samples and their replicates based on a specified variant frequency threshold.
Check_Original_Replicates <- function(CHOA_Muts,osID,rsID,min_var_rep){
  # Split the data into original and replicate mutations.
  Original_Muts <- CHOA_Muts[CHOA_Muts$subject %in% osID, ]
  Replicate_Muts <- CHOA_Muts[CHOA_Muts$subject %in% rsID, ]
  # Initializes a new column 'Presence_corrected' to 'N' in the original mutations dataframe for tracking mutation reproducibility.
  Original_Muts[, Presence_corrected := "N"]
  # Updates 'Presence_corrected' to 'Y' for mutations in original samples that are confirmed by replicate samples meeting the frequency threshold.
  Original_Muts$Presence_corrected <- ifelse(Original_Muts$RepFreq >= min_var_rep, "Y", "N")
  return(Original_Muts)
}

# Identifies samples with poor mutation reproducibility between original and replicate samples based on a frequency threshold.
Get_Bad_IDs <- function(CHOA_Muts, osID, rsID,min_var_rep,ID_version){
  # Initializes a dataframe to track sample reproducibility details
  SampleReproducibility_data <- data.frame(
    Basic_ID = character(),
    SampleID = character(),
    ReplicateID = character(),
    niSNV = integer(),
    nSharediSNV = integer(),
    pShared = double(),
    stringsAsFactors = FALSE
  )
  # Assess reproducibility of mutations between original and replicate samples
  Original_HiFreq_Muts <- Check_Original_Replicates(CHOA_Muts,osID,rsID,min_var_rep=min_var_rep)
  
  # Iterate through each original sample ID to evaluate shared mutations
  ix <- 1
  for (ix in 1:length(osID)) {
    SampleID <- osID[ix]
    ReplicateID <- rsID[ix]
    Basic_ID <- sub("_.*", "", SampleID) # Extract the base ID without additional identifiers
    
    # Gather mutations marked as 'Y' in Presence_corrected for each sample
    mutations <- Original_HiFreq_Muts[subject == SampleID, ]
    
    # Count total and shared mutations to compute the share proportion
    niSNV <- nrow(mutations)
    nSharediSNV <- sum(mutations$Presence_corrected == "Y")
    
    # Calculate the proportion of shared mutations
    pShared <- nSharediSNV / niSNV
    
    # Compile reproducibility data for each sample
    SampleReproducibility_data <- rbind(SampleReproducibility_data, 
                                        data.frame(Basic_ID, SampleID, ReplicateID, niSNV, nSharediSNV, pShared))
  }
  # Plot reproducibility data highlighting the proportion of mutations not shared
  graph <-
    ggplot(SampleReproducibility_data, aes(x = Basic_ID, y = 1 - pShared)) +
    geom_col(fill = "darkgrey") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    xlab("ID") +
    ylab("p(Not Shared)") +
    ggtitle("Sample Reproducibility Data")
  
  ggsave(filename = paste0(Figure_OutPut_path,"/",ID_version,"reproducibility_plot.pdf"), plot = graph)
  
  # Determine samples with less than 50% mutation sharing and flag as poor reproducibility
  bad_OG_samples <- SampleReproducibility_data %>%
    filter(!is.na(pShared) & pShared < 0.5 & niSNV > 1)
  
  # Extract the basic IDs of these samples
  Bad_IDs <- bad_OG_samples$Basic_ID
  
  return(Bad_IDs)
}

# Function to visualize and analyze the relationship between initial and replicate iSNV frequencies using a scatter plot and linear regression.

Plot_muts_pub_arrangmement <- function(iMuts, Analysis_date) {
  # Prepare a subset of the data focusing on necessary variables for plotting.
  plot_muts <- iMuts[, .(Sample_ID_mut, var_freq, RepFreq)]
  
  # Fit a linear model to predict replicate variant frequency based on initial variant frequency.
  fit <- lm(RepFreq ~ var_freq, data = plot_muts)
  
  # Calculate the R-squared value from the linear model fit to represent the goodness of fit.
  r_squared <- summary(fit)$r.squared
  
  # Generate a scatter plot to visually compare variant frequencies between original and replicate samples.
  var_comp <- ggplot2::ggplot(plot_muts, aes(x = var_freq, y = RepFreq)) +
    geom_point(alpha = 1, color = "black") + # Use black points to represent data points.
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +   # Add a linear regression line without confidence interval.
    labs(x = "iSNV frequency", y = "Replicate iSNV frequency") +  # Label axes.
    scale_x_continuous(limits = c(0.02, 0.505), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) +  # Set X-axis limits and breaks.
    scale_y_continuous(limits = c(0.02, 0.505), breaks = seq(0, 0.5, 0.1), expand = c(0, 0)) + # Set Y-axis limits and breaks.
    theme_classic() +  # Use the classic theme for a clean look.
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),  # Customize the panel border.
          panel.grid.major = element_blank(),  # Remove major grid lines.
          panel.grid.minor = element_blank(), # Remove minor grid lines.
          axis.text = element_text(size = 10),  # Set size of axis text.
          axis.title = element_text(size = 10))  # Set size of axis titles.
  
  # Display the created scatter plot.
  print(var_comp)
  
  # Define the filename for the saved plot
  scatter_name <- paste0(Analysis_date,"-Replicate_Frequency_Comparison_VPC.pdf")
  scatter_eps <- paste0(Analysis_date,"-Replicate_Frequency_Comparison_VPC.eps")
  
  # Save the scatter plot as a PDF file.
  ggplot2::ggsave(file.path(Figure_OutPut_path,scatter_name),plot = var_comp, width = 6, height = 5)
  
  # Save the scatter plot as an EPS file.
  ggplot2::ggsave(file.path(Figure_OutPut_path,scatter_eps), plot = var_comp, width = 6, height = 5, device = "eps")
  
  # Return the R-squared value to indicate the fit of the linear model.
  return(var_comp)
}


# Compile SampleTrackers -----------------------------------------------------
# Initialize an empty data frame to compile sample trackers from various run dates
SampleTracker_Compiled <- data.frame()
sd <- 1

# Loop through each run date specified in the RunDate vector
for (sd in 1:length(RunDate)) {
  # Construct the file path for the sample tracker CSV based on the run date
  tracker_name <-
    paste("PostProcessing_Output/",RunDate[sd], "_SampleTracker_CHOA.csv", sep = "")
  
  # Load the sample tracker file
  tracker_tmp <- read.csv(tracker_name)
  
  tracker_tmp <-
    tracker_tmp[-c(1)] # Get rid of the first index column
  # Append the loaded data to the main compiled sample tracker
  SampleTracker_Compiled <-
    rbind(SampleTracker_Compiled, tracker_tmp)
  # Reset row names to ensure consistency and avoid potential issues with default row naming
  rownames(SampleTracker_Compiled) <- NULL
}

# Define column names explicitly for clarity and ensure they are correctly assigned
colnames(SampleTracker_Compiled) <-
  c("Sample",
    "Strain",
    "SeqDate",
    "PB2_accession",
    "propN",
    "QC_Cov")

# Save the complete list of compiled sample trackers to a CSV file for all data
if (length(unique(SampleTracker_Compiled$Sample)) != nrow(SampleTracker_Compiled)) {
  stop("There are duplicate samples IDs that need to be fixed before we proceed.")
}

# Save the compiled sample trackers
SampleTracker_Compiled_path <- file.path(Compiled_OutPut_path,"SampleTracker_Compiled_ALL.csv")
write.csv(SampleTracker_Compiled, SampleTracker_Compiled_path)

# Filter the compiled sample trackers to include only those samples passing the QC for coverage
SampleTracker_Compiled <-
  SampleTracker_Compiled[SampleTracker_Compiled$propN < 5, ]

# Extract passing sample identifiers
passingSamples <- SampleTracker_Compiled$Sample

# Save the filtered list of sample trackers to a separate CSV file
SampleTracker_Compiled_path <- file.path(Compiled_OutPut_path,"SampleTracker_Compiled.csv")
write.csv(SampleTracker_Compiled,SampleTracker_Compiled_path)


# Compile Mutation Output -----------------------------------------------------
# Initialize an empty data frame to compile mutation outputs
MutationOutput_Compiled <- data.frame()
md <- 1

# Loop through each date provided in RunDate to process mutation output files
for (md in 1:length(RunDate)) {
  # Construct the file name for the current mutation output based on the run date
  mutation_name <-
    paste0("PostProcessing_Output/MutationOutput_", RunDate[md], ".csv")
  # Read the mutation output CSV file
  mutation_tmp <- read.csv(mutation_name)
  # Append the data from the current file to the compiled mutation output dataframe
  MutationOutput_Compiled <-
    rbind(MutationOutput_Compiled, mutation_tmp)
}
# Remove any redundant column named 'X' that might be added during file reading
MutationOutput_Compiled$X <- NULL

# Filter the compiled mutation outputs to include only those from samples with low propN values
MutationOutput_Compiled <-
  MutationOutput_Compiled[MutationOutput_Compiled$subject %in% passingSamples, ]

# Create Fasta Copies -----------------------------------------------------

# Create a main directory for masked consensus fasta files if it doesn't already exist
dir.create("PostProcessing_OutPut/MaskedConsensus_CHOP")
rd <- 1
# Loop through each run date to manage Fasta files
for (rd in 1:length(RunDate)) {
  run <- RunDate[rd]
  # Define the source directory where original Fasta files are stored
  src_dir <-
    paste0(R_path, "/Pipeline_Output_Files/",RunDate[rd],"_PipelineOutput/MaskedConsensus")
  
  # Define the destination directory where Fasta files will be copied
  dest_dir <-
    paste0(R_path, "/PostProcessing_OutPut/MaskedConsensus_CHOP/", run, "_MaskedConsensus")
  
  # Ensure the destination directory exists; create if not
  dir.create(dest_dir)
  
  # Retrieve a list of all files in the source directory
  files <- list.files(src_dir)
  # Copy all files from the source to the destination directory, allowing overwrites
  file.copy(file.path(src_dir, files),
            file.path(dest_dir, files),
            overwrite = TRUE)
}


# Fix Mutations -----------------------------------------------------
# Set default state for 'Modified' to 0 indicating no changes have been made yet
MutationOutput_Compiled$Modified <- 0
# Remove row names from the compiled mutation output for cleaner data handling
rownames(MutationOutput_Compiled) <- NULL
# Initialize data frames to track rows to delete and add during mutation fixing
to_delete <- data.frame()
to_add <- data.frame()

# Select mutations that are classified as intrahost and have a variant frequency greater than 50%
target_muts <-
  MutationOutput_Compiled[MutationOutput_Compiled$var_freq > 0.5 &
                            MutationOutput_Compiled$variant_level == "Intra", ]

# Loop through the list of mutations and update the mutation and consensus sequence
MutationOutput_Unmodified <- MutationOutput_Compiled

# Preserve the original state of mutation output before making any modifications and save as a csv
MutationOutput_Unmodified_path <- paste0(Compiled_OutPut_path,"/MutationOutput_Unmodified.csv")
write.csv(MutationOutput_Unmodified, MutationOutput_Unmodified_path)
mut <- 1
# Start processing each mutation in the target list
for (mut in 1:nrow(target_muts)) {
  print(nrow(target_muts) - mut)  # Debug: Show remaining mutations to process
  mut_tmp <- target_muts[mut, ]
  row_number <- as.numeric(rownames(mut_tmp))  # Get the row number for potential updates
  
  # Extract the current segment_num, position, and var_nt
  segment_num <- as.integer(mut_tmp$segment_num)
  position <- as.integer(mut_tmp$position)
  var_nt <- mut_tmp$var_nt
  ref_nt <- mut_tmp$ref_nt
  
  # Load sample information from the tracker
  sample_tmp <-
    SampleTracker_Compiled[SampleTracker_Compiled$Sample %in% mut_tmp$subject, ]
  
  # Load the appropriate reference sequence based on the strain
  if (sample_tmp$Strain == "H3N2") {
    Ref <- H3N2_Ref
  } else if (sample_tmp$Strain == "H1N1") {
    Ref <- H1N1_Ref
  } else {
    stop("Non-existing strain!!!") # Stop execution if strain is not recognized
  }
  
  # Path to the fasta files for sequence adjustments
  masked_path <-
    paste0(R_path,
           "/PostProcessing_Output/MaskedConsensus_CHOP/",
           sample_tmp$SeqDate,
           "_MaskedConsensus")
  
  target_fasta <-
    intersect(
      list.files(path = masked_path, pattern = mut_tmp$subject),
      list.files(masked_path, pattern = ".fasta")
    )
  
  # Load the fasta data for the targeted mutation
  tmp_fasta <-
    phylotools::read.fasta(file.path(masked_path, target_fasta))
  
  # Retrieve the sequence for the specific segment and determine the nucleotide at the mutation position
  seq <- tmp_fasta$seq.text[segment_num]
  
  # Check if the nucleotide at the mutation position matches the reference
  target_nucleotide <-
    substring(tmp_fasta$seq.text[segment_num], position, position)
  target_nucleotide
  
  # Make sure the nucleotide isn't "N", if it is the the nucleotide is masked anyways
  if (target_nucleotide != "N") {
    # Update the sequence to include the variant nucleotide
    if (target_nucleotide == ref_nt) {
      # Replace the nucleotide at the specified position
      seq_modified <- stringr::str_sub(seq, 1, position - 1) %>%
        stringr::str_c(var_nt) %>%
        stringr::str_c(stringr::str_sub(seq, position + 1))
      
      # Update the fasta sequence in the dataset
      tmp_fasta$seq.text[segment_num] <- seq_modified
      
      # Save the modified fasta sequence back to the file
      updated_fasta <-
        dplyr::data_frame(name = tmp_fasta$seq.name, seq = tmp_fasta$seq.text) #This is weird formatting, but it's the only way it works
      writeFasta(data = updated_fasta, file.path(masked_path, target_fasta)) 
      rm(updated_fasta, seq_modified) # Clean up variables
    } else{
      stop("The reference nt's don't match!")
    }
  }
  
  # Evaluate the mutation type and perform updates or deletions as needed
  # First, check if there's a consensus mutation at the population level that matches the target mutation characteristics
  ref_seq <- Ref$NT[segment_num]
  ref_seq <- gsub("U", "T", ref_seq)
  consensus_ref_nt <-
    substring(ref_seq, position, position)
  consensus_ref_nt
  
  # See if there's an existing population level mutation for the current mutation
  pop_tmp <-
    MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                              MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                              MutationOutput_Compiled$position %in% mut_tmp$position &
                              MutationOutput_Compiled$variant_level == "Pop", ]
  if (nrow(pop_tmp) == 0) {
    # If no corresponding population-level mutation is found, reclassify or update the mutation as necessary
    if (consensus_ref_nt != mut_tmp$var_nt) {
      # If the consensus reference does not match the variant nucleotide, update the mutation level to population
      mut_updated <- mut_tmp
      mut_updated$variant_level <- "Pop"
      mut_updated$Modified <- 1 # Mark the mutation entry as modified
      MutationOutput_Compiled[row_number, ] <- mut_updated
      rm(mut_updated)
    } else if (consensus_ref_nt == mut_tmp$var_nt) {
      # If the mutation matches the consensus but was incorrectly classified, schedule it for deletion
      row_number <- as.numeric(rownames(mut_tmp))
      to_delete <- rbind(to_delete, row_number)
    }
  }
  
  # Also check if the mutation can be reclassified as an intrahost variation that hasn't been documented yet
  swap_vf <- 1 - mut_tmp$var_freq # Calculate the alternate variant frequency
  
  if (swap_vf >= 0.01) {   # Consider for reclassification only if >1%
    # Check existing entries for the same mutation classified as intrahost
    intra_tmp <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                                MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                                MutationOutput_Compiled$position %in% mut_tmp$position &
                                MutationOutput_Compiled$var_nt == mut_tmp$var_nt &
                                MutationOutput_Compiled$variant_level == "Intra", ]
    # Check for other intrahost records to evaluate depth information
    depth_check <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject %in% mut_tmp$subject &
                                MutationOutput_Compiled$segment_num %in% mut_tmp$segment_num &
                                MutationOutput_Compiled$position %in% mut_tmp$position &
                                MutationOutput_Compiled$variant_level == "Intra", ]
    if (nrow(intra_tmp) >= 1 && nrow(pop_tmp) == 1) {
      # If the mutation is already listed and has a corresponding population record
      update_tmp <- mut_tmp
      # Swap reference and variant nucleotides to update the mutation
      update_tmp$ref_nt <- mut_tmp$var_nt
      update_tmp$var_nt <- mut_tmp$ref_nt
      
      # Update the variant frequency to reflect intrahost variation
      update_tmp$var_freq <- 1 - mut_tmp$var_freq
      
      # Update amino acid references to reflect the swap
      update_tmp$aa_ref <- mut_tmp$aa_sub
      update_tmp$aa_sub <- mut_tmp$aa_ref
      update_tmp$aa_ref2 <- mut_tmp$aa_sub2
      update_tmp$aa_sub2 <- mut_tmp$aa_ref2
      update_tmp$variant_level <- "Intra"
      # Update variant depth based on available data
      if (nrow(depth_check) > 1) {
        update_tmp$var_depth <- NaN
      } else {
        update_tmp$var_depth <- mut_tmp$depth - mut_tmp$var_depth
      }
      update_tmp$Modified <- 1
      # Replace the original mutation entry with the updated one
      MutationOutput_Compiled[row_number, ] <- update_tmp
      rm(update_tmp)
    }
    else if (nrow(intra_tmp) == 1 && nrow(pop_tmp) == 0) {
      stop("check me out!") # Trigger error if an inconsistency is found
    } else if (nrow(intra_tmp) == 0) {
      # Side note - this will be negative because we've already changed the initial mutation from intra --> pop
      #This qualifies as an intrahost variant as well, so we need to generate one
      update_tmp <- mut_tmp
      # Swap and update information as previously described
      update_tmp$ref_nt <- mut_tmp$var_nt
      update_tmp$var_nt <- mut_tmp$ref_nt
      update_tmp$var_freq <- 1 - mut_tmp$var_freq
      update_tmp$aa_ref <- mut_tmp$aa_sub
      update_tmp$aa_sub <- mut_tmp$aa_ref
      update_tmp$aa_ref2 <- mut_tmp$aa_sub2
      update_tmp$aa_sub2 <- mut_tmp$aa_ref2
      update_tmp$variant_level <- "Intra"
      if (nrow(depth_check) > 1) {
        update_tmp$var_depth <- NaN
      } else {
        update_tmp$var_depth <- mut_tmp$depth - mut_tmp$var_depth
      }
      update_tmp$Modified <- 1
      update_tmp
      to_add <- rbind(to_add, update_tmp)   # Prepare to add the new mutation
      rm(update_tmp)
    } else {
      stop("Check me out!")  # Alert for unexpected cases
    }
  }
  rm(mut_tmp)  # Clean up the temporary mutation data
}

if (nrow(to_delete) > 0) {
  # Assign column name 'Row' to the dataframe 'to_delete' to facilitate row identification
  colnames(to_delete) <- "Row"
  
  # Sort the rows marked for deletion in descending order to prevent row shifting during deletion
  to_delete_sorted <- sort(to_delete$Row, decreasing = TRUE)
  
  # Iteratively remove rows from MutationOutput_Compiled based on the sorted deletion indices
  for (row in to_delete_sorted) {
    MutationOutput_Compiled <- MutationOutput_Compiled[-row,]
  }
}


# Append rows stored in 'to_add' to the MutationOutput_Compiled dataframe
MutationOutput_Compiled <- rbind(MutationOutput_Compiled, to_add)

# Sort MutationOutput_Compiled first by 'position', then by 'segment_num', and finally by 'subject'
# to organize the data for easier readability and processing
MutationOutput_Compiled <-
  MutationOutput_Compiled[order(MutationOutput_Compiled$position), ]
MutationOutput_Compiled <-
  MutationOutput_Compiled[order(MutationOutput_Compiled$segment_num), ]
MutationOutput_Compiled <-
  MutationOutput_Compiled[order(MutationOutput_Compiled$subject), ]

# Define the path for the final compiled mutation output
Mutation_Output_Compiled_path <- paste0(Compiled_OutPut_path,"/MutationOutput_Compiled.csv")

# Write the compiled mutation output to a CSV file at the specified path
write.csv(MutationOutput_Compiled, Mutation_Output_Compiled_path)

# Generate Consensus Sequences -----------------------------------------------------
# Extract the first 8 ("CHOA-XXX) characters from each sample identifier to separate unique IDs
Sample_ID <- stringr::str_sub(SampleTracker_Compiled$Sample, 1, 8)

# Deduplicate the Sample_IDs to ensure each unique identifier is only processed once
uID <- unique(Sample_ID)

# Create a directory to store the generated consensus sequences for organization
dir.create(file.path(R_path,"PostProcessing_OutPut/CompiledOutPut/Consensus_seqs"))

# Inspection of the generated consensus sequences for CHOP-009_D3V1,
# CHOA-009_D4V1, and CHOA-009_D5V1 indicated that all subjects contained a
# deletion in their consensus sequence not present on the other runs where the
# sample was sequenced (6/17/2022 and 5/01/2023). As it would be highly
# implausible for a deletion to be present at the beginning of the PA gene
# segment and because this deletion was not detected in the other 2 technical
# replicates, we have excluded the 3 sequences from subsequent analysis. We
# remove them here as the resulting frameshift form the deletion leads too many
# mismatches between the sample consensus sequences for PA and prevents
# reasonable evaluation of the remaining replicates.

# These samples are excluded from further analysis due to the impact such errors can have on data integrity.
SampleTracker_Compiled <- SampleTracker_Compiled %>%
  filter(!Sample %in% c("CHOA-009_D3V1", "CHOA-009_D4V1", "CHOA-009_D5V1"))

# Generate consensus sequences for all remaining samples using the defined unique IDs
PassingSample_Matches <- Generate_consensus_sequences(uID, SampleTracker_Compiled, R_path, segment_name)

# Identify samples with significant mismatches in their consensus sequences,
# which may indicate issues in sequence alignment or quality. 
Mismatched_Samples <-
  PassingSample_Matches[PassingSample_Matches$nMismatch > 100, ]

# Store the details of mismatched samples in a CSV file for potential reevaluation or correction.
Mismatched_Samples_path <- paste0(Compiled_OutPut_path,"/CHOA_Mismatched_Samples.csv")
write.csv(Mismatched_Samples,Mismatched_Samples_path)


# Exclude samples with excessive mismatches from the main compiled datasets to prevent skewed analysis results,
# planning to reassess these samples for potential recovery or reprocessing.
SampleTracker_Compiled <-
  SampleTracker_Compiled[!SampleTracker_Compiled$Sample %in% Mismatched_Samples$`passing_samples$Sample`, ]

MutationOutput_Compiled <-
  MutationOutput_Compiled[!MutationOutput_Compiled$subject %in% Mismatched_Samples$`passing_samples$Sample`, ]


# Repeat consensus generation --------------------------------------------------

# Extract sample IDs, replicate numbers, and version numbers from the compiled sample tracker
Sample_ID_revised <- stringr::str_sub(SampleTracker_Compiled$Sample, 1, 8)
uID_revised <- unique(Sample_ID_revised)

# Regenerate consensus sequences, now excluding samples with significant mismatches identified previously
PassingSample_Matches_Revised <- Generate_consensus_sequences(uID_revised, SampleTracker_Compiled, R_path, segment_name)

# Organize the compiled sample tracker and revised matches by sample ID for consistency and easier management
Passing_CompiledSamples <-
  SampleTracker_Compiled[order(SampleTracker_Compiled$Sample), ]
PassingSample_Matches_Revised <-
  PassingSample_Matches_Revised[order(PassingSample_Matches_Revised$`passing_samples$Sample`), ]
Passing_CompiledSamples
PassingSample_Matches_Revised

# Break down sample identifiers into separate components for detailed tracking and analysis
SampleID_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 1, 8)
Replicates_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 10, 11)
Versions_only <-
  stringr::str_sub(Passing_CompiledSamples$Sample, 12, 13)

# Remove row names for cleaner dataframe management
rownames(Passing_CompiledSamples) <-
  NULL
rownames(PassingSample_Matches) <- NULL

# Compile a detailed Quality Control (QC) dataframe from the sorted sample data
PassingSample_QC <-
  data.frame(
    Sample = Passing_CompiledSamples$Sample,
    ID = SampleID_only,
    Replicate = Replicates_only,
    Version = Versions_only,
    Strain = Passing_CompiledSamples$Strain,
    SeqDate = Passing_CompiledSamples$SeqDate,
    PB2_accession = Passing_CompiledSamples$PB2_accession,
    propN = Passing_CompiledSamples$propN,
    QC_Cov = Passing_CompiledSamples$QC_Cov
  )

# Verify that the samples in QC matches correspond to those in revised matches
if (!all(PassingSample_Matches_Revised$`passing_samples$Sample` == PassingSample_QC$Sample)) {
  stop("The rows are not equivalent")
}

# Combine QC data with matching sample data for a comprehensive overview
tmp_combo <- cbind(PassingSample_QC, PassingSample_Matches_Revised)

# Save the comprehensive QC data for future reference and analyses
CHOA_SampleQC <- tmp_combo
write.csv(CHOA_SampleQC, paste0(Compiled_OutPut_path,"/CHOA_SampleQC.csv"))

# Extract Samples that pass all QC metrics unambiguously
passing_CHOASamples <- CHOA_SampleQC[CHOA_SampleQC$nMismatch == 0, ]

# Identify the IDs for te unique samples 
uID_pass <- unique(passing_CHOASamples$ID)

# We also save those samples with mismatches as they may still meet quality control standards with manual inspection
mismatch_CHOASamples <- CHOA_SampleQC[CHOA_SampleQC$nMismatch > 0, ]
write.csv(mismatch_CHOASamples, paste0(Compiled_OutPut_path,"/mismatch_CHOASamples.csv"))
mismatch_CHOASamples <- read.csv(paste0(Compiled_OutPut_path,"/mismatch_CHOASamples.csv"))


# Select Replicates -----------------------------------------------------

# Initialize data frames for storing sample information
sample_mat <- data.frame()
selected_samples <- data.frame()
not_replicate_samples <- data.frame()
rs <- 1

# Iterate over each unique ID passed
for (rs in 1:length(uID_pass)) {
  uID_tmp <- uID_pass[rs]
  
  # Filter samples by ID and sort by 'propN' in ascending order for lower error rates
  uPassing <-
    passing_CHOASamples[passing_CHOASamples$ID == uID_tmp, ]
  uPassing <-
    uPassing[order(uPassing$propN, decreasing = FALSE), ]
  # Check for any duplicates within the filtered samples
  if (any(duplicated(uPassing$Sample))) {
    stop("Error: Duplicated samples found in 'Sample' column of 'uPassing'")
  }
  # Get distinct replicates based on ID and Replicate number, retaining all associated data
  distinctPassing <-
    uPassing %>% distinct(uPassing$ID, uPassing$Replicate, .keep_all = TRUE) # Get the rows with the unique replicates
  # Remove unnecessary columns added for sorting or other operations
  distinctPassing <-
    distinctPassing[-c(21, 22)] # Get rid of the two sorting columns
  
  # Decision block to handle whether sufficient replicates are available
  if (nrow(uPassing) == 1) {
    # Skip further processing if no replicates are present
    next
  }
  if (nrow(distinctPassing) >= 2) {
    # If at least two replicates are available, select them
    selected_samples_tmp <- distinctPassing[1:2, ]
    selected_samples_tmp
    print(selected_samples_tmp)
    selected_samples <- rbind(selected_samples, selected_samples_tmp)
  } else {
    # If not enough replicates, select the two samples with the lowest 'propN'
    uPassing <-
      uPassing[order(uPassing$propN, decreasing = FALSE), ]
    not_replicate_samples_tmp <-
      uPassing[1:2, ]
    not_replicate_samples <-
      rbind(not_replicate_samples, not_replicate_samples_tmp)
  }
}

# Output the selected samples and not selected for further analysis
CHOA_Selected_Samples <- selected_samples
write.csv(CHOA_Selected_Samples, paste0(Compiled_OutPut_path,"/CHOA_Selected_Samples.csv"))

CHOA_NotReplicated_Samples <- not_replicate_samples
write.csv(CHOA_NotReplicated_Samples,paste0(Compiled_OutPut_path,"/CHOA_NotReplicated_Samples.csv"))


# Identify Samples with mismatches that can still be considered ----------

# Determine samples that were not included in the selected sets, possibly due to mismatches or other QC issues
NotSelected_Samples <-
  dplyr::setdiff(CHOA_SampleQC$ID, CHOA_Selected_Samples$ID)

# A total of 48 samples were sequenced but not included in the final sample set
# due to the absence of a paired sample replicate. This does not include samples
# with propN > 5%.

# Extract unmatched samples
Unmatched_Samples <-
  CHOA_SampleQC[CHOA_SampleQC$ID %in% NotSelected_Samples, ]

write.csv(Unmatched_Samples, paste0(Compiled_OutPut_path,"/CHOA_UnMatched_Samples.csv"))

# We inspect the unmatched samples manually to see if the samples may still be
# included. The reason for automated exclusion is either: 1) no sequence
# replicate (which can't be fixed) or the presence of mismatches between the
# consensus sequences. We can consider consensus mismatches if they occur within
# 20bp of the segment ends, represent a mismatch between a called and masked
# nucleotide, are in a region of poor sequence quality, or occur positions with
# high variation. The final set of selected samples can be found in the CSV
# file: CHOA_Selected_Samples_revised.csv. See below for a summary of the
# samples we have chosen to include following manual inspection and the
# reasoning.


# CHOA-051_D1V1 and CHOA-051_D3V1. There was a single mismatch between the
# consensus sequences for these samples at NS 306. Inspection of the
# CHOA-051_D3V1 isolate showed a high variant (49%) iSNV at NS 306 with the
# variant allele matching the consensus for D1V1, so we will use these samples.

# CHOA-128_D1V1 and CHOA-128_D2V1. These can be kept as there was a single
# mimatch at PB2 813. Both samples had high frequency variants at that site.


# Fix inter/intrahost variants --------------------------------------------------------

# Read in the csv file that contains selected samples
CHOA_Selected_Samples <-
  read.csv(paste0(Compiled_OutPut_path,"/CHOA_Selected_Samples_Revised.csv"))

# Create a list of samples from the data frame
sample_list <- CHOA_Selected_Samples$Sample

# Initialize an empty data frame to store the mutation output
CHOA_MutationOutput <- data.frame()

# Loop over the sample list to process each sample individually
for (nv in 1:length(sample_list)) {
  # Print the number of remaining iterations to the console
  length(sample_list) - nv
  
  # Get the current sample from the list
  sample_tmp <- sample_list[nv]
  
  # Check if the current sample is valid and not a NaN (Not a Number)
  if (is.nan(sample_tmp) == FALSE) { # Ensures that the sample name is not NaN, which is critical for data integrity.
    # Extract the mutations present in that sample from the compiled mutation output data frame
    sample_muts <-
      MutationOutput_Compiled[MutationOutput_Compiled$subject == sample_tmp, ]
    
    # Extract the mutations present in that sample from the compiled mutation output data frame
    CHOA_MutationOutput <- rbind(CHOA_MutationOutput, sample_muts)
    
    # Remove the temporary data frame to free up memory
    rm(sample_muts)
  }
}

# After processing all samples, extract and process identifiers for further analysis
subject_mut <- CHOA_MutationOutput$subject

# Extract different components from the subject field for detailed analysis
Sample_ID_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 1, 8)
Replicates_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 10, 11)
Versions_mut <- stringr::str_sub(CHOA_MutationOutput$subject, 12, 13)

# Combine the extracted fields with the original mutation output data frame to form a comprehensive data structure
CHOA_mutations <-
  cbind(Sample_ID_mut,
        Replicates_mut,
        Versions_mut,
        CHOA_MutationOutput)

# Fix variants that have NA as a frequency. These were hard to calculate for positions with multiple alleles present

# Remove the row names from the mutation data frame
rownames(CHOA_mutations) <- NULL

# Find the rows with NA in the var_depth column
to_fix <- CHOA_mutations[is.na(CHOA_mutations$var_depth), ]

# Get the row indices of the rows to fix
indices <- as.numeric(rownames(to_fix))

# Calculate var_depth by multiplying var_freq by depth and round it off
to_fix$var_depth <- round(to_fix$var_freq * to_fix$depth)

# Condition to check if there are any rows to fix
if (isEmpty(to_fix) != TRUE){  # Checks if the to_fix dataframe is not empty, meaning there are rows that need adjustments.
  # Initialize a loop counter
  tf <- 1
  # Loop through each row that needs to be fixed
  for (tf in 1:nrow(to_fix)) {
    # Update the original CHOA_mutations data frame with fixed values
    ind_target <- indices[tf]
    
    # Directly replaces the problematic entries in the main data frame with the corrected ones from to_fix.
    CHOA_mutations[ind_target, ] <- to_fix[tf, ]
  }
}

# Filter the mutation data frame for entries labeled as 'Intra'
CHOA_intrahost <-
  CHOA_mutations[CHOA_mutations$variant_level == "Intra", ]

# Print the number of rows in the intra-host data frame to the console
nrow(CHOA_intrahost)

# Remove the row names from the intra-host data frame
rownames(CHOA_intrahost) <- NULL

# Write the intra-host mutation data to a CSV file
write.csv(CHOA_intrahost,paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))

# Find samples to exclude ------------------------------------------------------

# Load all mutations
CHOA_intrahost <- read.csv(paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))

# Filter samples into originals vs replicates 
Original_Samples <- CHOA_Selected_Samples %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup()
# Generate a list of samples tagged as "original"
osID <- Original_Samples$Sample

# Filter the second occurrence of each ID and store in Replicate_Samples
Replicate_Samples <- CHOA_Selected_Samples %>%
  group_by(ID) %>%
  slice_tail(n = 1) %>%
  ungroup()
# Generate a list of samples tagged as "replicate"
rsID <- Replicate_Samples$Sample


# Create a data frame to store information about each replicate mutation
replicate_mut_info <-
  data.frame(
    ReplicateSubject = character(),
    Presence = character(),
    RepFreq = double(),
    RepDepth = double()
  )

# Set the maximum frequency threshold for considering mutations
maxFreq <- 0.03

# Set the frequency threshold for replicate analysis
replicateFreq <- 0.03

# Set the minimum coverage threshold for mutation analysis
min_cov <- 100

# Filter the CHOA_intrahost data frame for mutations above the specified frequency threshold
CHOA_targets <-
  CHOA_intrahost[CHOA_intrahost$var_freq >= maxFreq, ]

# Define the path for pileup files used in mutation analysis
pileup_path <-
  paste0(R_path, "/PostProcessing_OutPut/Pileup_renamed/")

# Extract unique IDs from the selected samples for tracking or further analysis
Selected_IDs <- unique(CHOA_Selected_Samples$ID)

# Initialize an empty vector for IDs that may need to be excluded from analysis
Exclude_IDs <- c() # We haven't identifed samples to exclude yet

nm <- 1
# Loop through each target mutation to check for replicates
for (nm in 1:nrow(CHOA_targets)) {
  print(nrow(CHOA_targets) - nm) # Display the number of mutations left to process, useful for tracking progress during long operations.
  
  # Extract the current mutation for analysis
  mut_tmp <- CHOA_targets[nm, ]
  replicate_mut_info_this <- Check_ReplicateMut(CHOA_targets, mut_tmp, CHOA_Selected_Samples, pileup_path)
  
  # Append the information about replicate mutations to a collective dataframe
  replicate_mut_info <- rbind(replicate_mut_info,replicate_mut_info_this)
}

# Ensure the format of replicate information is numeric for frequency and depth
replicate_mut_info$RepFreq <-
  as.numeric(replicate_mut_info$RepFreq)
replicate_mut_info$RepDepth <-
  as.numeric(replicate_mut_info$RepDepth)

# Combine the target mutations with their corresponding replicate information into one table
CHOA_HiFreq_Muts <- as.data.table(cbind(CHOA_targets, replicate_mut_info))
CHOA_HiFreq_Muts$X <- NULL # Remove index column 

# Filter mutations to include only those with sufficient depth
CHOA_HiFreq_Muts <- CHOA_HiFreq_Muts[depth >= min_cov & RepDepth >= min_cov]

# Identify samples with problematic mutation frequencies using a custom function
Bad_IDs_1 <- Get_Bad_IDs(CHOA_HiFreq_Muts,osID, rsID,min_var_rep=0.03,ID_version = "OG")
Bad_IDs_2 <- Get_Bad_IDs(CHOA_HiFreq_Muts,rsID, osID,min_var_rep=0.03,ID_version = "Rep")

# Combine and deduplicate IDs from both rounds of analysis
Bad_IDs_all <- unique(c(Bad_IDs_1, Bad_IDs_2))

# Sort the IDs for better readability and further analysis
Bad_IDs_all <- Bad_IDs_all[order(Bad_IDs_all)]

# Output the bad IDs to a CSV for record-keeping or further examination
write.csv(Bad_IDs_all, paste0(Compiled_OutPut_path,"/CHOA_Bad_IDs_All.csv"))

# Save the finalized sample set by removing Bad_IDs_All from CHOA_Selected_Samples
CHOA_SampleSet_FINAL <- CHOA_Selected_Samples[!CHOA_Selected_Samples$ID %in% Bad_IDs_all,]
write.csv(CHOA_SampleSet_FINAL, paste0(Compiled_OutPut_path,"/CHOA_SampleSet_FINAL.csv"))


# Finalize mutations at chosen frequency parameters ------------------

# Load previously identified bad IDs and remove the default read.csv column naming
Bad_IDs_all <- read.csv(paste0(Compiled_OutPut_path,"/CHOA_Bad_IDs_All.csv"))
Bad_IDs_all$X <- NULL # Remove extraneous column created during CSV read

# Load intra-host mutation data and remove any extraneous column
CHOA_intrahost <- fread(paste0(Compiled_OutPut_path,"/CHOA_IntrahostCompiled.csv"))
CHOA_intrahost$V1 <- NULL # Remove extraneous column created during CSV read

# Exclude mutations related to the bad IDs from further analysis
Muts <- CHOA_intrahost[!Sample_ID_mut %in% Bad_IDs_all$x]

# Generate empty dataframe to contain information about replicates
replicate_mut_info <- data.frame()
# Loop through mutations to check for replicates using a custom function
nm <- 1
for (nm in 1:nrow(Muts)) {
  print(nrow(Muts) - nm) # Display remaining mutations to process for user awareness
  
  # Get the current mutation from the data frame
  mut_tmp <- Muts[nm, ]
  replicate_mut_info_this <- Check_ReplicateMut(Muts, mut_tmp, CHOA_Selected_Samples, pileup_path)
  
  # Store collected replicate information
  replicate_mut_info <- rbind(replicate_mut_info,replicate_mut_info_this)
}

# Ensure numerical representation for replication frequencies and depths for accurate computations
replicate_mut_info$RepFreq <-
  as.numeric(replicate_mut_info$RepFreq)
replicate_mut_info$RepDepth <-
  as.numeric(replicate_mut_info$RepDepth)

# Combine mutations with their replication details
CHOA_All_Muts <- as.data.table(cbind(Muts, replicate_mut_info))
CHOA_All_Muts$X <- NULL # Remove index column 

# Define minimum criteria for variant frequency and coverage
min_var_final <- 0.03
min_var_replicate_final <- 0.03
min_cov_final <- 100

# Apply filters to select mutations based on frequency and coverage criteria
Muts_this <- CHOA_All_Muts[var_freq >= min_var_final & RepFreq >= min_var_replicate_final]
CHOA_Muts_1 <- Check_Original_Replicates(Muts_this,osID,rsID,min_var_rep=min_var_replicate_final)
CHOA_Muts_1 <- CHOA_Muts_1[depth >= min_cov & RepDepth >= min_cov]

# Calculate average variant frequency and coverage for more robust analysis
CHOA_Muts_1[, avg_var_freq := (var_freq + RepFreq) / 2]
CHOA_Muts_1[, avg_cov := (depth + RepDepth)/2 ]

# Repeat the process for the reverse pairing of original and replicate samples
CHOA_Muts_2 <- Check_Original_Replicates(Muts_this,rsID,osID,min_var_rep=min_var_replicate_final)
CHOA_Muts_2 <- CHOA_Muts_2[depth >= min_cov & RepDepth >= min_cov]
CHOA_Muts_2[, avg_var_freq := (var_freq + RepFreq) / 2]
CHOA_Muts_2[, avg_cov := (depth + RepDepth)/2 ]

# Combine and deduplicate mutations from both original and replicate analyses
CHOA_Muts_FINAL <- rbind(CHOA_Muts_1,CHOA_Muts_2)
CHOA_Muts_FINAL_FILTERED <- distinct(CHOA_Muts_FINAL, Sample_ID_mut, position, segment_num, var_nt, .keep_all = TRUE)

# Assign mutation effect categories based on presence of any non-synonymous mutations
CHOA_Muts_FINAL <- CHOA_Muts_FINAL_FILTERED

# Assign mutation effect categories based on presence of any non-synonymous mutations
CHOA_Muts_FINAL$effect_any <- ifelse(CHOA_Muts_FINAL$effect == "Nonsynonymous" | CHOA_Muts_FINAL$effect2 == "Nonsynonymous", "Nonsynonymous", "Synonymous")

# Output final processed mutation data
write.csv(CHOA_Muts_FINAL,paste0(Compiled_OutPut_path,"/FINAL_MutationOutput_minvar3.csv"))

# Calculate summary statistics for iSNV frequencies
mean_iSNV_frequency <- mean(CHOA_Muts_FINAL$avg_var_freq)
std_iSNV_frequency <- sd(CHOA_Muts_FINAL$avg_var_freq)

# Generate a plot visualizing replicate concordance
replicate_concordance_plot <- Plot_muts_pub_arrangmement(iMuts = CHOA_Muts_FINAL, Analysis_date = "81123")

# Save the finalized sample set by removing Bad_IDs_All from CHOA_Selected_Samples
CHOA_SampleSet_FINAL <- CHOA_Selected_Samples[!CHOA_Selected_Samples$ID %in% Bad_IDs_all,]
write.csv(CHOA_SampleSet_FINAL, paste0(Compiled_OutPut_path,"/CHOA_Muts_FINAL.csv"))

# Save parameters chosen run parameters -----
# Save the final parameters used in the analysis into a CSV file for record-keeping.
CHOA_RunParam <- data.frame(
  min_var = min_var_final,
  min_var_replicate = min_var_replicate_final,
  min_cov = min_cov_final
)
write.csv(CHOA_RunParam, file = paste0(Compiled_OutPut_path,"/CHOA_RunParam.csv"), row.names = FALSE)

# Incorporate metadata --------------
# Load metadata from a clean CSV file and remove any extraneous index column that may have been imported.
CHOA_metadata <- read.csv("CLEAN_CHOA_metadata.csv")
CHOA_metadata$X <- NULL

# Merge the loaded metadata with sample data based on subject IDs to ensure each sample is linked with the correct metadata.
merged_data<- merge(CHOA_metadata, Original_Samples, by.x = "SubjectID", by.y = "ID")

# Add columns to the merged data to clearly label the original and replicate samples for each subject.
merged_data$Subject_Original <- Original_Samples$Sample
merged_data$Subject_Replicate <- Replicate_Samples$Sample


# Generate mutation frequency plot -------

# Create and display a histogram of average variant frequencies for identified mutations.
iSNV_freq_hist <- ggplot2::ggplot(CHOA_Muts_FINAL, aes(x = avg_var_freq)) + 
  geom_histogram(breaks = seq(0, 0.5, by = 0.01), fill = "white", color = "black") +  # Histogram bins defined by 0.01 increments
  labs(x = "iSNV frequency", y = "Count") +
  scale_x_continuous(limits = c(0, 0.51), breaks = seq(0, 0.5, by = 0.1), expand = c(0, 0)) +  # X-axis properties
  scale_y_continuous(limits = c(0, 62), breaks = seq(0, 60, by = 20), expand = c(0, 0)) +  # Y-axis properties
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Display the plot
print(iSNV_freq_hist)

# Save the plot in eps format
ggplot2::ggsave(filename = paste0(Figure_OutPut_path,"/iSNV_frequency_distribution.eps"), plot = iSNV_freq_hist, device = "eps", width = 6, height = 5, units = "in")


# Filter and save non-synonymous mutations for HA and NA segments into separate CSV files.
NonSynonymous_Muts_HA <- CHOA_Muts_FINAL[effect_any == "Nonsynonymous" & segment_num == 4]
write.csv(NonSynonymous_Muts_HA,paste0(Compiled_OutPut_path,"/NonSynonymous_Muts_HA.csv"))

NonSynonymous_Muts_NA <- CHOA_Muts_FINAL[effect_any == "Nonsynonymous" & segment_num == 6]
write.csv(NonSynonymous_Muts_HA,paste0(Compiled_OutPut_path,"/NonSynonymous_Muts_NA.csv"))

# Determine the number of iSNVs per sample ----------------------
# Count various types of mutations for each sample and update the merged data frame with these counts.
merged_data$all_iSNV_present <- 0
merged_data$nonSyn_iSNV_present <- 0
merged_data$Syn_iSNV_present <- 0
merged_data$nonSyn_HA_present <- 0
merged_data$nonSyn_NA_present <- 0

si <- 1
for (si in 1:length(merged_data$SubjectID)) {
  # count number of matching entries in Var_Present
  sample_id <- merged_data$SubjectID[si]
  sampleSNV_present <- CHOA_Muts_FINAL[CHOA_Muts_FINAL$Sample_ID_mut == sample_id,]
  sampleSNV_nonSynonymous <- CHOA_Muts_FINAL[CHOA_Muts_FINAL$Sample_ID_mut == sample_id & CHOA_Muts_FINAL$effect_any == "Nonsynonymous",]
  
  merged_data$all_iSNV_present[si] <- nrow(sampleSNV_present)
  merged_data$nonSyn_iSNV_present[si] <- nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous",])
  merged_data$Syn_iSNV_present[si] <- nrow(sampleSNV_present[sampleSNV_present$effect_any == "Synonymous",])
  merged_data$nonSyn_HA_present[si] <-  nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous" & sampleSNV_present$segment_num == 4,])
  merged_data$nonSyn_NA_present[si] <-  nrow(sampleSNV_present[sampleSNV_present$effect_any == "Nonsynonymous" & sampleSNV_present$segment_num == 6,])
  merged_data$mean_varFreq_all[si] <- mean(sampleSNV_present$avg_var_freq)
  merged_data$mean_varFreq_nonSyn[si] <- mean(sampleSNV_nonSynonymous$avg_var_freq)
}


# Generate association plot for all iSNVs and CT -----
# Generate and display a scatter plot to visualize the relationship between CT values and the number of iSNVs per sample.
CT_iSNV_plot <- ggplot2::ggplot(merged_data, aes(x = CT, y = all_iSNV_present)) + 
  geom_point(color = "black", fill = "black", shape = 21, size = 1.25) +  # Define the appearance of points
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add a linear model line without a confidence envelope
  scale_x_continuous(limits = c(15, 30), breaks = seq(15, 30, by = 5), expand = c(0, 0)) + # X-axis properties
  scale_y_continuous(limits = c(-2, 82), breaks = seq(0, 80, by = 10), expand = c(0, 0)) + # Y-axis properties
  labs(x = "CT", y = "iSNV per Sample") +  # Axis labels
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        text = element_text(size=12, family="Helvetica"),
        axis.title = element_text(size=12))
# Display the plot
print(CT_iSNV_plot)

# Fit a linear model to analyze the relationship between CT values and the count of iSNVs in each sample.
lm_model <- lm(all_iSNV_present ~ CT, data = merged_data)

# Calculate and print the R-squared value to assess the model's explanatory power.
r_squared <- summary(lm_model)$r.squared
print(paste("R-squared value:", r_squared))

# Save the CT vs. iSNV scatter plot in EPS format, providing high-quality vector graphics for publications.
ggplot2::ggsave(filename = paste0(Figure_OutPut_path,"/iSNV_CT_association.eps"), plot = CT_iSNV_plot, device = "eps", width = 6, height = 5, units = "in")

# Arrange multiple plots into a single figure using a grid layout for cohesive presentation.
final_QC_plot <- cowplot::plot_grid(replicate_concordance_plot, iSNV_freq_hist, CT_iSNV_plot, ncol = 1, nrow = 3, align = "v", rel_widths = c(1, 1))
print(final_QC_plot)

# Save the composite QC plot in EPS format for high-resolution output. This is Figure 1 in the paper
ggplot2::ggsave(paste0(Figure_OutPut_path,"/final_QC_plot.eps"), plot = final_QC_plot, width = 3.5, height = 7, device = "eps")


# Generate final merged data set ------------

# Save the merged data frame containing comprehensive QC data to a CSV file. This will be used for generating the phylogeny
write.csv(merged_data,paste0(Compiled_OutPut_path,"/FINAL_CHOA_merged_data.csv"))

# Reload the saved data to ensure accuracy and completeness.
merged_data <- read.csv(file.path("PostProcessing_OutPut/CompiledOutPut/FINAL_CHOA_merged_data.csv"))

# Create a data frame summarizing key variables for each sample, such as total iSNVs, mean variant frequency, and counts of non-synonymous and synonymous mutations.
sample_DataTable <- data.frame(SubjectID = merged_data$SubjectID,all_iSNV = merged_data$all_iSNV_present, 
                               mean_VarFreq = merged_data$mean_varFreq_all,
                               nonSynonymous_iSNV = merged_data$nonSyn_iSNV_present,
                               Synonymous_iSNV = merged_data$Syn_iSNV_present,
                               mean_NonSynon_VarFreq = merged_data$mean_varFreq_nonSyn,
                               nonSynoymous_HA = merged_data$nonSyn_HA_present,
                               nonSynonymous_NA = merged_data$nonSyn_NA_present)

# Order the data by the count of iSNVs in decreasing order for better visibility of high-variant samples.
sample_DataTable <- sample_DataTable[order(sample_DataTable$all_iSNV, decreasing = TRUE),]
row.names(sample_DataTable) <- NULL

# Save the ordered table to a CSV file.
Name_datatable <- paste0(Compiled_OutPut_path,"/","FinalMuts-CHOA_SampleDataTable.csv")
write.csv(sample_DataTable,Name_datatable)


