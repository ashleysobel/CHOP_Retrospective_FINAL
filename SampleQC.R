# Title: Pipeline for Processing NGS Data for Influenza Analysis
# Description: This script processes NGS data to identify influenza variants,
# perform quality control, and generate necessary output directories and citation files.
# The script ensures reproducibility and efficient data management by handling file operations,
# setting QC parameters, and characterizing mutations based on reference sequences.
# Inputs: Path to the NGS data, run dates for sample processing.
# Outputs: Processed data files organized in specified directories, citation file for referenced R packages.
# Author: Ashley Sobel Leonard
# Date: 5/1/2024


# Clear environment to get rid of old files and get path -------------------
rm(list = ls())
R_path <- getwd()

# Set QC parameters -----------------------
RunDate <- c(22323) # Specify the RunDate(s) of interest
min_Qual <- 30 #Minimum PHRED score
min_Map <- 40  #Minimum mapping quality
min_coverage <- 100 #minimum coverage
options(warn=1) # Set warning level to show first warning

# Generate directory for holding the processed pipeline output ------
PostProcessingOutput_path <- file.path(R_path,"PostProcessing_OutPut")

dir.create(PostProcessingOutput_path)

# Load packages ---------- 
library(stringr) # Used for string manipulation functions like str_replace_all
library(phylotools)   # Used for reading fasta files (read.fasta)
library(tidyverse)  # General data manipulation and plotting functions
library(stringi) #  # String manipulation tools such as stri_extract_last
library(seqinr)  # Biological sequence analysis, including read and write fasta such as write.fasta and translate
library(tidysq) # Used for handling biological sequences, such as as.sq
library(patchwork)  # Used for combining ggplot2 plots

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("stringr", "phylotools", "tidyverse","stringi","seqinr","tidysq","patchwork")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","SampleQC_citations.bib"))

# Load reference data for the script ----------  
# Load Reference Key (matches PB2 accession # with strain). Make sure this has been updated if your run contains different reference files
Reference_key <- read.csv("Reference_files/Reference_key.csv")

# Set segment names
segment_name <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
# Enter the csv file containing the references you want to use. For the CHOP
# 2017-2018 retrospective study, we want to use H3N2 2022 and H1N1 2019.
H3N2_Ref <- read.csv("Reference_files/H3N2_Ref_A_2016.csv")
H1N1_Ref <- read.csv("Reference_files/H1N1_Ref_A_2015.csv")
H3N2_Ref_consensus <- phylotools::read.fasta("Reference_files/A_Washington_17_2016_H3N2.fasta")
H1N1_Ref_consensus <- phylotools::read.fasta("Reference_files/A_Michigan_45_2015_H1N1_ref.fasta")



# Define functions ------------------------------------------

# Define a function to retrieve fasta file information from a given directory path
# This function collects all fasta filenames within a directory and returns these details in a list
GetFasta <- function(PipelineOutput_path) {
  path_fasta <- paste0(PipelineOutput_path,"/Consensus"); path_fasta
  fasta_names <- list.files(path_fasta, pattern = "\\.fasta"); fasta_names
  fasta_info <- list(path_fasta=path_fasta,fasta_names=fasta_names); fasta_info
}

# Define a function to retrieve coverage data from a specified directory
# It collects names of coverage files, used to further analyze data quality
GetCov <- function(R_path, Pipeline_date) {
  path_basecov <- paste0(PipelineOutput_path,"/BaseCov"); path_basecov
  base_cov <- list.files(path_basecov, pattern = "basecov.txt"); base_cov
  cov_info <- list(path_basecov=path_basecov,basecov_names=base_cov)
}

# Function to create a path for coverage plots by sample name
# It constructs file paths for saving coverage plot PDFs by sample
GetCovPlot_path <- function(R_path,Pipeline_date,sample_name) {
  CovPlot_path <- paste0(CovPlot_path,"/",sample_name,"_Coverage.pdf")
} 

# Function to identify differences between two nucleotide or amino acid sequences
# This function is crucial for mutation analysis in genomic sequencess
GetStringDiff<- function(a, b, exclude = c("-", "?"), ignore.case = TRUE, show.excluded = FALSE){
  if(nchar(a)!=nchar(b)) stop("Lengths of input strings differ. Please check your input.")
  if(ignore.case)
  {
    a <- toupper(a)
    b <- toupper(b)
  }
  split_seqs <- strsplit(c(a, b), split = "")
  only.diff <- (split_seqs[[1]] != split_seqs[[2]])
  only.diff[
    (split_seqs[[1]] %in% exclude) |
      (split_seqs[[2]] %in% exclude)
  ] <- NA
  diff.info<-data.frame(which(is.na(only.diff)|only.diff),
                        split_seqs[[1]][only.diff],split_seqs[[2]][only.diff])
  names(diff.info)<-c("position","poly.seq.a","poly.seq.b")
  if(!show.excluded) diff.info<-na.omit(diff.info)
  diff.info
}


# Function to plot read coverage for a specific segment
# It visualizes read depth across a genome segment to assess sequencing coverage adequacy
Plot_Segment_Coverage <- function(sorted_base_cov, segment_index, segment_name) {
  # Subset data to include only the current segment based on index
  segment <- subset(sorted_base_cov, sorted_base_cov$Segment == unique(sorted_base_cov$Segment)[segment_index])
  segment_Position <- as.numeric(segment$Position)
  segment_Coverage <- as.numeric(segment$Coverage)
  
  # Create a data frame specifically for plotting
  segment_df <- data.frame(cbind(segment_Position,segment_Coverage))
  
  # Generate the coverage plot using ggplot2
  segment_plot <- ggplot(segment_df, aes(x=segment_Position, y=segment_Coverage)) + 
    geom_area() +
#    scale_y_log10() + # add this line to set y-axis to log scale
    ggtitle(segment_name[segment_index]) + 
    theme(plot.title=element_text(hjust=0.5)) +
    labs(x = "Position",y = "Coverage") 
  
  #  Returning the ggplot object for further manipulation or saving outside this function.
  return(segment_plot)
}

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


# Function to characterize mutations, categorizing them based on their impact on the protein sequence
# This function analyzes genetic mutations found in sequencing data,
# determining the type and effects of mutations at the nucleotide and amino acid levels.
# Inputs:
# - stmp: Data frame containing mutation data for a given sample
# - subject: Identifier for the subject from whom the sample was taken
# - Strain: Viral strain of the influenza virus
# - Ref: Data frame containing reference genetic sequences
# Outputs: data frame containing the refernece genetic sequence 

Characterize_mutations <- function(stmp, subject, Strain, Ref) {
  sample_mat <- list()  # Initialize a list to hold mutation data for each mutation
  mutation_output_tmp <- data.frame()  # Temporary data frame to accumulate results for each sample.
  
  # Filter the dataset for substitution type mutations only, which are the primary focus here
  subs <- stmp[stmp$TYPE == "SUB", ] # Filter mutations for substitutions only 
  nSub <- nrow(subs) # Number of substitutions
  
  # Process each substitution if there are any
  if (nSub != 0){
    mt <- 1
    for (mt in 1:nSub){
      sub_tmp <- subs[mt,]; sub_tmp
      # Extract important identifiers like accession number and segment number from mutation data.
      accession <- strsplit(sub_tmp$CHROM,'-')[[1]][1]  # Extract accession number
      segment_num <- as.numeric(strsplit(sub_tmp$CHROM,'_')[[1]][3])  # Extract segment number
      position <- as.numeric(sub_tmp$POS) # Mutation position in the sequence
      ref_nt <- sub_tmp$REF  # Reference nucleotide
      var_nt <- sub_tmp$ALT  # Variant nucleotide
      var_freq <- as.numeric(sub_tmp$AF) # Variant frequency
      var_qual <- as.numeric(sub_tmp$QUAL)  # Quality score of the variant
      
      # Handle different data fields for mapping quality from variant calling outputs.
      if (is.null(sub_tmp$MQ)){
        sub_tmp$MQ <- sub_tmp$MQM  # Update missing MQ field using MQM if MQ is null.
      }
      var_map <- as.numeric(sub_tmp$MQ) # Mapping quality of the variant
      var_type <- sub_tmp$TYPE  # Type of the variant, expected to be "SUB" here.
      if (var_type != "SUB"){
        stop("There is an error!!!!!") # We're only considering substitutions for now, so this would be bad
      }
      depth <- as.numeric(sub_tmp$DP)  # Total depth at the mutation position.
      var_depth <- as.numeric(sub_tmp$AD) # Depth of the variant nucleotide.
      var_info <- sub_tmp$DP4 # Additional detailed information on variant depth
      
      # Compile all mutation information into a data frame.
      mut <- data.frame(subject,Strain,segment_num,position,ref_nt,var_nt,var_freq,var_qual,var_map,var_type,depth,var_depth,var_info)
      
      # Analyze the effect of the mutation based on the reference sequence.
      # This includes determining whether the mutation occurs in a coding region and its impact.
      Ref_seg <- Ref[Ref$Segment == mut$segment_num, ]
      if (mut$position < Ref_seg$Start) {
        # Handle mutations outside the coding regions.
        location <- "NonCoding"
        effect <- NaN
        aa_pos <- NaN
        aa_sub  <- NaN
        aa_ref <- NaN
        
      } else if (mut$position > Ref_seg$Stop) {
        # Variant is outside the coding region here as well
        location <- "NonCoding"
        effect <- NaN
        aa_pos <- NaN
        aa_sub  <- NaN
        aa_ref <- NaN
      } else {
        # Process mutations within coding regions to determine their impact on the protein sequence.
        location <- Ref_seg$Protein # Specify protein affected
        seq_ref <- Ref_seg$NT # Reference nucleotide sequence
        
        # Manipulate sequences to reflect mutations and translate them to assess the impact on the amino acid sequence.
        ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")  # Convert sequence to sq object
        ref_CDS <- bite(ref_nt_seg, indices = Ref_seg$Start:Ref_seg$Stop)  # Extract coding DNA sequence
        ref_AA <- c2s(seqinr::translate(s2c(as.character(ref_CDS))))   # Translate to amino acids
        
        seq_sub <- seq_ref # Copy reference sequence to introduce mutation.
        substr(seq_sub,mut$position,mut$position) <- mut$var_nt # Introduce mutation
        sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc") # Convert mutated sequence
        sub_CDS <- bite(sub_nt, indices = Ref_seg$Start:Ref_seg$Stop) # Extract mutated CDS
        sub_AA <- c2s(seqinr::translate(s2c(as.character(sub_CDS)))) # Translate mutated CDS
        
        aa_pos <- floor((mut$position - Ref_seg$Start)/3) + 1  # Calculate amino acid position
        aa_ref <- substr(x = ref_AA, start = aa_pos, stop = aa_pos) # Reference amino acid
        aa_sub <- substr(x = sub_AA, start = aa_pos, stop = aa_pos) # Mutated  amino acid
        
        # Determine the nature of the mutation (synonymous or nonsynonymous).
        if (aa_ref == aa_sub){
          effect <- "Synonymous"
        } else if (aa_ref != aa_sub){
          effect <- "Nonsynonymous"
        }
      }  
      # Compile mutation effects into a data frame
      mutation_effect1 <- data.frame(location,effect,aa_pos,aa_ref,aa_sub); mutation_effect1
      
      # Check for a second protein encoded in the same gene segment (if applicable)
      if (Ref_seg$Protein2 ==""){
        # Determine the intervals of the second protein within the gene segment.  If there is not another protein, than the related outputs are set to NaN as placeholders
        location2 <- NaN; aa_pos2 <- NaN; effect2 <- NaN; aa_ref2 <- NaN; aa_sub2 <- NaN;
      } else {
        # IF there is a second portien, the process is similar to above 
        # Retrieve the start and stop positions for potentially multiple coding regions of the second protein.
        intervals <- c(Ref_seg$Start2_1,Ref_seg$Stop2_1,Ref_seg$Start2_2,Ref_seg$Stop2_2); intervals
        
        # Check if the mutation falls within the coding region of the second protein
        if (mut$position >= intervals[1] && mut$position <= intervals[2]){
          # Handling mutations within the first coding segment of the second protein.
          location2 <- Ref_seg$Protein2  # Set the protein affected by the mutation.
          seq_ref <- Ref_seg$NT # Retrieve the reference nucleotide sequence.
          seq_sub <- seq_ref  # Duplicate the reference sequence to create a mutated sequence.
          position_check_seq <- seq_sub  # Copy for position checking.
           
          # Introducing and marking the mutation position for visual tracking in sequence analysis.
          substr(seq_sub,mut$position,mut$position) <- mut$var_nt  # Introduce the mutation.
          substr(position_check_seq,mut$position,mut$position) <- 'Q' #  Locating the mutation in the combined sequence to correctly calculate its impact.
          
          if (is.na(intervals[3]) == FALSE){
            # Handle the protein encoded in multiple disjoint segments
            seq_part1 <- substr(position_check_seq,intervals[1],intervals[2]) # Extract the portion for the first part of the protein sequence
            seq_part2 <- substr(position_check_seq,intervals[3],intervals[4]) # Extract the portion for the 2nd part of the protein sequence
            seq_whole <- paste(seq_part1,seq_part2,sep=""); seq_whole # combine the two sections 
            position_new <- unlist(gregexpr('Q', seq_whole)); position_new # Identify the position of the placeholder mutatioon
            
            ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc") # Extract the reference nucleotide sequence
            ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
            ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2)))) # Convert nucleotides to amino acids
            sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc") # identify the subject's nucleotide sequence 
            sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract then CDS from the subject's nucleotide sequence
          } else {
            # Simplified handling for a single segment coding sequence within the second protein.
            seq_part1 <- substr(position_check_seq,intervals[1],intervals[2]) 
            seq_whole <- seq_part1; seq_whole
            position_new <- unlist(gregexpr('Q', seq_whole)); position_new 
            
            ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
            ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2])) # Extract the CDS using the indices
            ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2))))
            sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
            sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2]))
          }
          sub_AA2 <- c2s(seqinr::translate(s2c(as.character(sub_CDS2)))) # Here we translate the subjects CDS to then amino acid sequence
          aa_pos2 <- floor((position_new - 1)/3) + 1 # Identify the position of the amino acid
          aa_ref2 <- substr(x = ref_AA2, start = aa_pos2, stop = aa_pos2) # Reference amino acid
          aa_sub2 <- substr(x = sub_AA2, start = aa_pos2, stop = aa_pos2) # Subject amino acid
          # Determine if the mutation is synonymous or nonsynonymous by comparing the refererence and subject's amino acid at the location
          if (ref_AA2 == sub_AA2){
            effect2 <- "Synonymous"
          } else if (ref_AA2 != sub_AA2){
            effect2 <- "Nonsynonymous"
          }
        } else if (is.na(intervals[3]) == FALSE  && mut$position >= intervals[3] && mut$position <= intervals[4]){
          # Variant is within a coding region for the 2nd protein. We repeat the same process as above
          location2 <- Ref_seg$Protein2
          seq_ref <- Ref_seg$NT # Get the sequence of the reference
          seq_sub <- seq_ref # Create the sample sequence
          position_check_seq <- seq_sub
          substr(seq_sub,mut$position,mut$position) <- mut$var_nt
          substr(position_check_seq,mut$position,mut$position) <- 'Q'
          seq_part1 <- substr(position_check_seq,intervals[1],intervals[2])
          seq_part2 <- substr(position_check_seq,intervals[3],intervals[4])
          seq_whole <- paste(seq_part1,seq_part2,sep=""); seq_whole
          position_new <- unlist(gregexpr('Q', seq_whole)); position_new
              # Translate the nucleotide sequences to amino acids.

          ref_nt_seg <- as.sq(seq_ref, alphabet = "dna_bsc")
          ref_CDS2 <- bite(ref_nt_seg, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
          ref_AA2 <- c2s(seqinr::translate(s2c(as.character(ref_CDS2))))
          sub_nt <- as.sq(seq_sub, alphabet = "dna_bsc")
          sub_CDS2 <- bite(sub_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4]))
          sub_AA2 <- c2s(seqinr::translate(s2c(as.character(sub_CDS2))))
          aa_pos2 <- floor((position_new - 1)/3) + 1 # Get aa residue (Thanks to Andrew for code!)
          aa_ref2 <- substr(x = ref_AA2, start = aa_pos2, stop = aa_pos2) # Reference amino acid
          aa_sub2 <- substr(x = sub_AA2, start = aa_pos2, stop = aa_pos2) # Subject amino acid
          
          #  Assessing the synonymous or nonsynonymous nature of the mutation based on amino acid changes.
          if (ref_AA2 == sub_AA2){
            effect2 <- "Synonymous"
          } else if (ref_AA2 != sub_AA2){
            effect2 <- "Nonsynonymous"
          }
        }
        else {
          location2 <- "NonCoding"
          effect2 <- NaN
          aa_pos2 <- NaN
          aa_sub2 <- NaN
          aa_ref2 <- NaN
        }
      }
      # Now we add the effect of the mutation for the 2nd coding sequence
      mutation_effect2 <- data.frame(location2,effect2,aa_pos2,aa_ref2,aa_sub2); mutation_effect2
      
      # We convert the mutation information generated above into a data frame format 
      mutation_mat <- data.frame(mut,mutation_effect1,mutation_effect2)
      # Append the mutation data frame to the list
      sample_mat[[mt]] <- mutation_mat
    }
  }
  # Combine all the data frames in the list
  mutation_output_tmp <- do.call(rbind, sample_mat)
  return(mutation_output_tmp)
}

# Function 'create_dir': Ensures a directory exists at the specified path or creates it if it does not exist.
create_dir <- function(parent_path, dir_name) { # Helper function to create directories and return their paths.
  full_path <- file.path(parent_path, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path)
  }
  return(full_path)
}

# Function to generate a reference for the sequence of interest
Generate_Reference <- function(Ref, consensus_sequence, segment_name_num) {
  Ref$NT <- NaN; Ref$AA <- NaN; Ref$AA2 <- NaN #Clear original sequences, need to replace w/consensus
  
  # Load the consensus sequence and use it to update the reference
  i <- 1
  for (i in 1:8){
    consensus_tmp <- consensus_sequence[i,]
    consensus_tmp
    # Convert sequence to character vector
    consensus_char <- as.character(consensus_tmp$seq.text)
    consensus_char
    segment_tmp <- as.double(segment_name_num$SegNum[i]); segment_tmp
    Ref$NT[segment_tmp] <- consensus_char #Update NT sequence
    
    Ref_seg <- Ref[segment_tmp,]; Ref_seg
    ref_nt <- as.sq(consensus_char, alphabet = "dna_bsc")
    
    # Update first protein sequence 
    ref_CDS <- bite(ref_nt, indices = Ref_seg$Start:Ref_seg$Stop); ref_CDS # Extract the CDS using the indices
    ref_AA <- c2s(seqinr::translate(s2c(as.character(ref_CDS))))
    Ref$AA[segment_tmp] <- ref_AA
    
    #If there's a second protein, update the second AA sequence
    if (Ref_seg$Protein2 != ""){
      # Check if the mutation falls within the coding region of the second protein
      intervals <- c(Ref_seg$Start2_1,Ref_seg$Stop2_1,Ref_seg$Start2_2,Ref_seg$Stop2_2); intervals  
      if (is.na(intervals[3])){
        # There is a single interval
        ref_CDS_2 <- bite(ref_nt, indices = c(intervals[1]:intervals[2]))
      } else{
        # There are multiple starts
        ref_CDS_2 <- bite(ref_nt, indices = c(intervals[1]:intervals[2],intervals[3]:intervals[4])) # Extract the CDS using the indices
      }
      Ref$AA2[segment_tmp] <- c2s(seqinr::translate(s2c(as.character(ref_CDS_2))))
    } 
  }
  return(Ref)
}

# Run QC ------------------------------------------
# Initialize RD to 1 before entering the loop (we will iterate over RD, a placeholder variable for the index, to assist with debugging)
RD <- 1
# Loop through each date in the RunDate vector
for (RD in 1:length(RunDate)){
  # Assign the current date from RunDate to Pipeline_date for this iteration
    Pipeline_date <- RunDate[RD]
    
    # Construct a path for pipeline output and create necessary directories for output
    PipelineOutput_path <- paste0(R_path,"/Pipeline_Output_Files/",Pipeline_date,"_PipelineOutput"); PipelineOutput_path
    CovPlot_path <- create_dir(PipelineOutput_path, "CovPlot")
    MaskedConsensus_path <- create_dir(PipelineOutput_path, "MaskedConsensus")
    
    # Initialize vectors to store various types of data collected during processing
    prop_N <- numeric()
    sample_list <- character()
    pass_Coverage <- character()
    sequencing_date <- character()
    sequence_strain <- character()
    PB2_accession_list <- character()
    
    # Load fasta and coverage data from the specified pipeline output path
    fasta_info <- GetFasta(PipelineOutput_path)
    cov_info <- GetCov(PipelineOutput_path)
    
    # Extract sample names that appear in both fasta and coverage data
    cov_conserved_names <- sub("_basecov\\.txt$", "", cov_info$basecov_names)
    fasta_conserved_names <- sub("_consensus_sequence\\.fasta$", "", fasta_info$fasta_names)
    not_common_entries <- !(fasta_conserved_names %in% cov_conserved_names) # Finding the entries that exist in both vectors
    fasta_info$fasta_names <- fasta_info$fasta_names[!not_common_entries] # Remove entries from fasta_info that are not in common

    # Loop through each fasta file to process and rename
    n <- 1 # Again, we index over a placeholder variable
    for (n in 1:length(fasta_info$fasta_names)){
      parts <- strsplit(fasta_info$fasta_names[n],"_")[[1]]
      sample_list[n] <- paste0(parts[[1]],"_",parts[[2]])
      sequencing_date[n] <- as.character(as.integer(Pipeline_date))
    }
    
    # Continue processing each sample in sample_list
    ind <- 1 # Another placeholder variable 
    for (ind in 1:length(sample_list)){
      # Process fasta files for each sample
      fasta_name_path <- paste0(fasta_info$path_fasta,"/",fasta_info$fasta_names[ind])
      tmp_fasta <- phylotools::read.fasta(fasta_name_path) # Phylotools package
      
      # Order fasta sequences by the segment number extracted from the sequence name
      qzar_last <- as.numeric(stringr::str_sub(tmp_fasta$seq.name, start = -1, end = -1))
      tmp_fasta_new <- cbind(tmp_fasta,qzar_last)
      tmp_qzar <- tmp_fasta_new[order(tmp_fasta_new$qzar_last),]
      tmp_fasta <- data.frame(seq.name=tmp_qzar$seq.name,seq.text=tmp_qzar$seq.text)
      # Check if the first sequence corresponds to an unintended reference and remove it
      if (tmp_fasta$seq.name[1] == "CY087752_1_4"){
        # This sequence mapped to H5N1, so we need to get rid (this is a control because we should not be seeing bird flu. This code initially written before 2024)
        to_delete <- sample_list[ind]
        file.remove(list.files(pattern=to_delete, recursive=TRUE))
        next
        
      }
      
      # Rename fasta files with new identifiers and update names within the file
      sample_name <- sample_list[ind]
      fasta_name_string <- paste0(MaskedConsensus_path,"/",sample_name,"_Masked.fasta")
      
      # Rename the individual gene segments
      old_name <- tmp_fasta$seq.name
      new_name <- paste(sample_name,segment_name,sep = "_") # Generate new names
      ref2 = data.frame(old_name,new_name)
      new_fasta <- paste(sample_name,"_",Pipeline_date,".fasta",sep = "")
      
      # Process coverage data to mask low coverage areas in the consensus sequence
      base_cov_path <- paste(cov_info$path_basecov,"/",cov_info$basecov_names[ind],sep="")
      base_cov <- read.delim(base_cov_path)
      colnames(base_cov) <- c("Segment","Position","Coverage")
      sorted_base_cov <- base_cov[order(stringi::stri_extract_last(base_cov$Segment,regex="\\w")),]
      
      # Convert the coverage and position columns to numeric so they can be manipulated
      base_cov$Coverage <- as.numeric(base_cov$Coverage)
      base_cov$Position <- as.numeric(base_cov$Position) 
      
      # Adjust position index (0-based to 1-based indexing)
      base_cov$Position <- base_cov$Position + 1
      
      # Save the updated Base Coverage file with corrected position as a CSV file
      new_base_cov_path <- gsub('basecov.txt','updated_basecov.csv',base_cov_path)
      write.csv(base_cov,new_base_cov_path)
      
      # Identify and mask positions with insufficient coverage
      coverage <- as.double(base_cov$Coverage)
      coverage_fail <- which(coverage < min_coverage) # Identify the positions with insufficient coverage, termed min_coverage
      # Update the consensus sequence so the sites where the coverage is < min_cov, are masked with "N"
      segx <- 1 # Another dummy variable for indexing
      for (segx in 1:8){ 
        seg_accession <- tmp_fasta$seq.name[segx] # Determine the accession number for segment and extract the base coverage for that segment
        seg_base_cov <- base_cov[base_cov$Segment == seg_accession,]
        seg_cov <- seg_base_cov$Coverage # Get the coverage for a segment
        
        bad_sites <- which(seg_cov < min_coverage) # Find positions where coverage is below minimum
        if (length(bad_sites) != 0){
          fasta_seg <- tmp_fasta[tmp_fasta$seq.name == seg_accession,] # Get fasta sequence for current segment
          seg_seq <- fasta_seg$seq.text
          
          for (i in bad_sites){
            # Mask positions with 'N' where coverage is insufficient
            stringr::str_sub(seg_seq,start = i, end = i) <-"N"    
          }
          tmp_fasta$seq.text[segx] <- seg_seq   # Update the sequence in the fasta object
        } 
      }
      # Save the masked fasta sequence 
      tmp_fasta$seq.name <- ref2$new_name 
      new_fasta_df <- dplyr::data_frame(name=tmp_fasta$seq.name,seq=tmp_fasta$seq.text) #This is weird formatting, but it's the only way it works
      writeFasta(data = new_fasta_df,filename = fasta_name_string)
      
      # Calculate the proportion of bases with insufficient coverage
      N_prop <- length(coverage_fail)/length(coverage) 
      prop_N[ind] <- N_prop*100
      if (prop_N[ind] > 5){
        pass_Coverage[ind] <- "No" # Mark as failing coverage quality control if more than 5% of bases are masked
      } else {
        pass_Coverage[ind] <- "Yes"  # Otherwise, mark as passing
      }
      accession <- unique(sorted_base_cov[,1])
      
      # Extract PB2 accession number from sequence metadata (this is the accession number I've chosen to index by to make sure I have contained the appropriate refernce file 
      PB2_accession <- unlist(strsplit(accession[1],'_'))[1]
      
      # Check if the PB2 accession is listed and assign it along with the corresponding strain
      accession_index <- Reference_key[Reference_key$Accession_PB2 == PB2_accession,]
      if (grepl(PB2_accession,accession_index$Accession_PB2,fixed=TRUE) != TRUE){
        stop("Need to update accession list")
      }
      
      if (is_empty(accession_index$Accession_PB2) == FALSE){
        PB2_accession_list[ind] <- accession_index$Accession_PB2
        sequence_strain[ind] <- accession_index$Strain
      } else {
        PB2_accession_list[ind] <- NaN
        sequence_strain[ind] <- NaN
      }
      
      # Generate segment-wise coverage plots      
      for (ns in 1:length(segment_name)){
        tmp_cov <- sorted_base_cov[sorted_base_cov[,1]==accession[1],]
        sorted_base_cov <- data.frame(lapply(sorted_base_cov, function(x){
          gsub(accession[ns],segment_name[ns],x)
        }))
      }
      
      # Generate segment-wise coverage plots
      # Each segment is processed individually to create a coverage plot
      PB2_plot <- Plot_Segment_Coverage(sorted_base_cov,1,segment_name)
      PB1_plot <- Plot_Segment_Coverage(sorted_base_cov,2,segment_name)
      PA_plot <- Plot_Segment_Coverage(sorted_base_cov,3,segment_name)
      HA_plot <- Plot_Segment_Coverage(sorted_base_cov,4,segment_name)
      NP_plot <- Plot_Segment_Coverage(sorted_base_cov,5,segment_name)
      NA_plot <- Plot_Segment_Coverage(sorted_base_cov,6,segment_name)
      MP_plot <- Plot_Segment_Coverage(sorted_base_cov,7,segment_name)
      NS_plot <- Plot_Segment_Coverage(sorted_base_cov,8,segment_name)
      
      # Saving the arranged plots to a file, typically a PDF, which provides high-quality output suitable for presentations and publications.
      sample_output <- GetCovPlot_path(R_path, Pipeline_date, sample_name)
      
      # Set the dimensions and file path for the PDF output
      pdf(file = sample_output,   # The directory you want to save the file in
          width = 8, # The width of the plot in inches
          height = 10) # The height of the plot in inches
    
      # Set the name for the coverage map  
      samplelabel <- paste("Coverage plots for", sample_name, sep = " ")

      # Arrange the plots using patchwork
      arranged_plot <- (PB2_plot | PB1_plot ) /
        (PA_plot | HA_plot ) /
        ( NP_plot | NA_plot) /
        (MP_plot | NS_plot) +
        patchwork::plot_annotation(title = samplelabel, theme = theme(plot.title = element_text(hjust = 0.5)))
      
      # Print the arranged plot to the device
      print(arranged_plot)
      
      dev.off() # Close the PDF device to finalize the file
    } 
    
    
    # Generate sample tracker containing all samples
    SampleTracker_all <- data.frame(sample_list,sequence_strain,sequencing_date,PB2_accession_list,prop_N,pass_Coverage)
    colnames(SampleTracker_all) <- c("Sample","Strain","Sequencing Date","Ref_Accession_PB2","%N","Pass QC_Cov")
    rownames(SampleTracker_all) <- NULL
    
    # Save sample tracker as csv file
    tracker_name <- paste0("PostProcessing_OutPut/",as.character(as.integer(Pipeline_date)),"_SampleTracker_all.csv")
    write.csv(x = SampleTracker_all,file = tracker_name)
    
    # Generate a sample tracker for just the CHOA samples
    SampleTracker_CHOA <- SampleTracker_all[str_detect(SampleTracker_all$Sample, "CHOA"),]
    SampleTracker_CHOA <- SampleTracker_CHOA[is.na(SampleTracker_CHOA$Strain) == FALSE,]
    CHOA_tracker_name <- paste0("PostProcessing_OutPut/",as.character(as.integer(Pipeline_date)),"_SampleTracker_CHOA.csv")
    write.csv(x = SampleTracker_CHOA,file = CHOA_tracker_name)
    
    
    # Characterize sample variants -----------------------------
    H3N2_Ref$NT <- gsub("U", "T", H3N2_Ref$NT) # Convert to DNA
    H3N2_Ref$Segment <- as.numeric(H3N2_Ref$Segment)
    H1N1_Ref$NT <- gsub("U", "T", H1N1_Ref$NT) # Convert to DNA
    H1N1_Ref$Segment <- as.numeric(H1N1_Ref$Segment)
    
    # assign path to variants identified in Pipeline
    consensus_variant_path <- file.path(PipelineOutput_path,"Variants/Consensus")
    intrahost_variant_path <- file.path(PipelineOutput_path,"Variants/Intrahost")
    
    # Generate a list of samples so that we can classify the variants for each of the samples.
    samples <- SampleTracker_CHOA$Sample
    N_sample <- length(samples)
    
    mutation_output <- data.frame() # Empty dataframe for holding mutation info 
    segment_name_num <- data.frame(SegNum <- as.double(c(1:8)),Segment <- c("PB2","PB1","PA","HA","NP","NA","MP","NS"))
    colnames(segment_name_num) <- c("SegNum","Segment")
    
    # We consider both intrahost variants (variants present at less than 50%
    # frequency in the subjects, these are called relative to the corresponding
    # subject's consensus sequence) and consensus level variants (variants that
    # are present >50% frequency defined relative to the reference sequence)
    dir.create(file.path(PipelineOutput_path,"Reference"))
    dir.create(file.path(PipelineOutput_path,"Reference","Intrahost"))
    dir.create(file.path(PipelineOutput_path,"Reference","Consensus"))
    
    ### Run variant classification loops
    new_reference_intra_path <- paste0(PipelineOutput_path,"/Reference/Intrahost"); new_reference_intra_path
    new_reference_pop_path <- paste0(PipelineOutput_path,"/Reference/Consensus"); new_reference_pop_path
    
    # We iterate over each sample to identify the consensus level and intrahost level variants.
    smp <- 1 # another dummy variable for indexing 
    for (smp in 1:N_sample){
      # Get information for this sample
      sample_tmp <- samples[smp]
      subject_parts <- unlist(strsplit(sample_tmp, "_"))
      subject <- paste(subject_parts[1],subject_parts[2],sep="_") # Get subject name
      sample_info <- SampleTracker_CHOA[SampleTracker_CHOA$Sample == subject,] # Get the subject information from the sample tracker
      Strain <- sample_info$Strain # Get subject's strain
      
      # Ensure there are no duplicate entries, this would be an error
      if(nrow(sample_info) > 1){
        print("There is a duplicate sample, this is an error")
        {stop(TRUE)}
      }
      Strain <- sample_info$Strain # Get subject's strain
      
      # Load the masked consensus file for the current subject
      name_path <- file.path(MaskedConsensus_path,paste0(sample_info$Sample,"_Masked.fasta")); name_path
      
      # Generate reference file for the current sample (either H3N2 or H1N1)
      if (is.na(Strain)){
        print("There is no matching strain.")
        {stop(TRUE)}
      }
      if (Strain == "H3N2"){
        # The subject's strain is H3N2, so we load the information for the H3N2 reference sequence
        Ref <- H3N2_Ref
        consensus_sequence_pop <- H3N2_Ref_consensus
        consensus_sequence_tmp <- tidysq::read_fasta(name_path,alphabet="dna")  
        # We set the reference for identifying intrahost variants as that subject's consensus sequence 
        consensus_sequence_intrahost <- tibble(seq.name = consensus_sequence_tmp$name,seq.text = consensus_sequence_tmp$sq)
      } else if (Strain == "H1N1"){
        # The subject's strain is H1N1, so we load the information for the H1N1 reference sequence
        Ref <- H1N1_Ref
        consensus_sequence_pop <- H1N1_Ref_consensus
        consensus_sequence_tmp <- tidysq::read_fasta(name_path,alphabet="dna") 
        # We set the reference for identifying intrahost variants as that subject's consensus sequence 
        consensus_sequence_intrahost <- tibble(seq.name = consensus_sequence_tmp$name,seq.text = consensus_sequence_tmp$sq)
      }
      

      # Gerenate a reference sequence for the subject for both the consensus level (Ref_pop) and intrahost level (Ref_intrahost)
      Ref_pop <- Generate_Reference(Ref, consensus_sequence_pop, segment_name_num)
      Ref_intrahost <- Generate_Reference(Ref, consensus_sequence_intrahost, segment_name_num)
      
      # Save sample's reference file as a csv
      Sample_Ref_name_pop <- paste(Strain,'_Ref','_',subject,"_pop",sep=""); Sample_Ref_name_pop
      Sample_Ref_name_intra <- paste(Strain,'_Ref','_',subject,"_intra",sep=""); Sample_Ref_name_intra
      
      # Save the reference sequences 
      write.csv(x=Ref_pop,file=paste0(new_reference_pop_path,'/',Sample_Ref_name_pop,'.csv'))
      write.csv(x=Ref_intrahost,file=paste0(new_reference_intra_path,'/',Sample_Ref_name_intra,'.csv'))
      
      
      # Characterize mutations in the target sample from consensus and intrahost variant data
      # List files matching the sample identifier within specified variant directoriese
      stmp_path_consensus <- list.files(path=consensus_variant_path,pattern = samples[smp])
      stmp_path_intrahost <- list.files(path=intrahost_variant_path,pattern = samples[smp])
      
      # Read consensus variant data, set 'ALT' column as character to handle genetic codes properly
      stmp_consensus <- read.csv(file.path(consensus_variant_path,stmp_path_consensus))
      original_colclasses <- lapply(stmp_consensus,class) # Store original data types of columns
      original_colclasses$ALT <- "character"  # Ensure 'ALT' column is treated as character
      stmp_consensus <- read.csv(file.path(consensus_variant_path,stmp_path_consensus),colClasses = original_colclasses)
      
      # Read intrahost variant data with similar handling for 'ALT' column
      stmp_intra <- read.csv(file.path(intrahost_variant_path,stmp_path_intrahost))
      original_colclasses <- lapply(stmp_intra,class)
      original_colclasses$ALT <- "character"   # Ensure 'ALT' column is treated as character
      stmp_intra <- read.csv(file.path(intrahost_variant_path,stmp_path_intrahost),colClasses = original_colclasses)
      
      # Characterize variants using custom function and assign variant levels for downstream analysis
      mutation_output_pop <- Characterize_mutations(stmp_consensus,subject, Strain, Ref_pop)
      mutation_output_pop$variant_level <- "Pop"
      
      mutation_output_intrahost <- Characterize_mutations(stmp_intra,subject, Strain, Ref_intrahost)
      
      if (is_empty(mutation_output_intrahost) == FALSE){
        mutation_output_intrahost$variant_level <- "Intra"
        mutation_output_tmp <- rbind(mutation_output_pop,mutation_output_intrahost)
      } else {
        mutation_output_tmp <- mutation_output_pop
      }
      # Order mutations by segment number and position for standardized reporting
      mutation_output_tmp <- mutation_output_tmp[
        with(mutation_output_tmp, order(segment_num,position)),
      ]
      # Combine mutation data from current processing batch to the main output
      mutation_output <- rbind(mutation_output, mutation_output_tmp)
    }
    
    # Save the final mutation output to a CSV file in the designated post-processing output directory
    mutation_ouptut_name <-
      paste("MutationOutput", '_', Pipeline_date, '.csv', sep = "")
    mutation_output_path <- file.path(PostProcessingOutput_path,mutation_ouptut_name)
    write.csv(mutation_output, mutation_output_path)

}

