# This code generates the information needed for uploading files to the SRA toolkit.

# Prepare environment -----------------------------------------------------
rm(list = ls())
R_path <- getwd()

SRA_info_path <- file.path(R_path,"SRA_info")
dir.create(SRA_info_path)

# Load Packages -----------------------------------------------------
library(dplyr)
library(tidyverse)
library(stringr)
library(lubridate)

# Citations information for script --------------------

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("dplyr","tidyverse")

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
writeLines(all_citations_text, con = file.path("ScriptCitations","SRA_info.bib"))


# Define Functions -----------------------------------------------------

# Define a function called GetFastq_name that takes sample_tmp as an input and returns the date the sample was run and the full sample name
GetFastq_name <- function(sample_tmp) {
  # Construct the path to the sample folder using R_path and the subdirectory "Pipeline_Output_Files".
  sample_folder <- file.path(R_path, "Pipeline_Output_Files")
  
  # List all files in the sample folder, including files in subdirectories (recursive = TRUE).
  files <- list.files(sample_folder, recursive = TRUE)
  
  # Filter files that contain sample_tmp followed by any characters and "consensus_sequence.fasta" in their names.
  matching_files <- files[grep(paste0(sample_tmp, ".*", "consensus_sequence.fasta"), files)]
  matching_files
  # Print the list of matching files (for debugging purposes).
  matching_files
  
  # Check the number of matching files:
  if (length(matching_files) == 0) {
    # If there are no matching files, stop the function and display an error message.
    stop("The matching file cannot be found")
  } else if (length(matching_files) > 1) {
    # If there are multiple matching files, print the list and stop the function with an error message.
    print(matching_files)
    stop("There are multiple potential matches")
  } else if (length(matching_files == 1)) {
    # If there is exactly one matching file:
    
    # Extract the seqDate from the matching file name using regular expression.
    seqDate <- sub("^(\\d+).*", "\\1", matching_files)
    
    # Extract the sample_tmp from the matching file name using regular expression.
    full_name_tmp <- sub(".*/([^/]+)_consensus_sequence.fasta", "\\1", matching_files)
    # Create a vector containing seqDate and the input sample_tmp.
    fastq_components <- c(seqDate, full_name_tmp)
  }
  
  # Return the list containing seqDate and sample_tmp as the result.
  return(fastq_components)
}


# Function to convert date strings to Date objects using lubridate
convert_to_date <- function(date_string) {
  # Determine the length of the date string
  date_length <- nchar(date_string)
  
  if (date_length == 5) {
    # For date strings like '12123' (January 21, 2023)
    year <- substr(date_string, 4, 5)
    month <- substr(date_string, 1, 1)
    day <- substr(date_string, 2, 3)
  } else if (date_length == 4) {
    # For date strings like '123' (January 2, 2023)
    year <- substr(date_string, 3, 4)
    month <- substr(date_string, 1, 1)
    day <- substr(date_string, 2, 2)
  } else {
    # Handle other cases or throw an error
    stop("Invalid date format")
  }
  
  # Reformat the string to 'YYYY-MM-DD'
  formatted_date <- paste0("20", year, "-", 
                           str_pad(month, width = 2, pad = "0"), "-", 
                           str_pad(day, width = 2, pad = "0"))
  
  # Convert to Date
  as.Date(formatted_date, format = "%Y-%m-%d")
}



# Load data ------------------------------------------------------------

# Load csv file containing information about the final sample set used in the
# CHOA analysis
CHOA_merged_data <- read.csv(file.path(R_path,"CHOA_Phylogenies/FINAL_CHOA_merged_data.csv"))

# Load the csv files containing the sampling date information
# viruses
seq_date_info <- read.csv(file.path(R_path,"CHOA_Phylogenies/CHOA_SeqDate.csv"))

All_IAV_metadata <- data.frame("SubjectID" = seq_date_info$Unique.Study.Subject.ID,"collection_date"=seq_date_info$Sample.date ) 

# Combine the CHOA_merged_data and the All_IAV_metadata into a single dataframe by joining on the SubjectID column entry
combined_data <- merge(All_IAV_metadata, CHOA_merged_data, by = "SubjectID")

# Extract the sample names by combining the entries in the Subject_Original and Subject_Replicate columns
Original_CHOA_Samples <- combined_data$Subject_Original
Replicate_CHOA_Samples <- combined_data$Subject_Replicate
View(Original_CHOA_Samples)
# Record information regarding teh reference fasta files used for assembly
Ref_H1N1pdm <- "A/Michigan/45/2015(H1N1), Accession MK622934-MK622941 "
Ref_H3N2 <- "A/Washington/17/2016(H3N2), Accession KX414254-KX414261"

# # Generate SRA sample attributes table -----------------------------------------------------

# Create the tibble for containing the SRA sample attributes table. The format is based on: Package Pathogen: clinical or host-associated; version 1.0
CHOA_Original_attributes <- tibble(
  sample_name = character(),
  sample_title = character(),
  bioproject_accession = character(),
  organism = character(),
  strain = character(),
  isolate = character(),
  collected_by = character(),
  collection_date = character(),
  geo_loc_name = character(),
  host = character(),
  host_disease = character(),
  isolation_source = character(),
  lat_lon = character(),
  culture_collection = character(),
  genotype = character(),
  host_age = character(),
  host_description = character(),
  host_disease_outcome = character(),
  host_disease_stage = character(),
  host_health_state = character(),
  host_sex = character(),
  host_subject_id = character(),
  host_tissue_sampled = character(),
  passage_history = character(),
  pathotype = character(),
  serotype = character(),
  serovar = character(),
  specimen_voucher = character(),
  subgroup = character(),
  subtype = character(),
  description = character(),
  replicate = double()
)

# Generate the same data holder for the replicate samples
CHOA_Replicate_attributes <- CHOA_Original_attributes

# Loop over the original and replicate sample entries to populate the CHOA_attributes tibble
os<- 1
for (os in 1:length(Original_CHOA_Samples)){ # where cs is index for the CHOA_sample
  sample_tmp <- Original_CHOA_Samples[os]; sample_tmp

  # Get the subject information for the current sample by extracting the first component of the name
  sample_name_tmp <- unlist(strsplit(x = sample_tmp,"_"))[1]
  sample_name_tmp
  
  # Extract the relevant entry from combined_data dataframe
  this_sample <- combined_data[combined_data$SubjectID == sample_name_tmp,]
  CHOA_Original_attributes[os,]$replicate <- 1
  CHOA_Original_attributes[os,]$sample_name <- sample_tmp
  CHOA_Original_attributes[os,]$sample_title <- sample_name_tmp
  CHOA_Original_attributes[os,]$bioproject_accession <- ""
  CHOA_Original_attributes[os,]$organism <- "Influenza A virus"
  CHOA_Original_attributes[os,]$strain <- this_sample$Strain
  CHOA_Original_attributes[os,]$subtype <- this_sample$Strain
  # Indicate that the isolate is the "original" by adding the _1 suffix for the isolate info
  CHOA_Original_attributes[os,]$isolate <- paste0(this_sample$SubjectID,"_1")
  CHOA_Original_attributes[os,]$collected_by <- "Children's Hopsital of Philadelphia"
  CHOA_Original_attributes[os,]$collection_date <- this_sample$collection_date
  CHOA_Original_attributes[os,]$geo_loc_name <- "USA: Pennsylvania: Philadelphia"
  CHOA_Original_attributes[os,]$lat_lon <- '39.9487 N 75.1939 W'
  CHOA_Original_attributes[os,]$host <- "Homo sapiens"
  CHOA_Original_attributes[os,]$host_disease <- "Infuenza"
  CHOA_Original_attributes[os,]$isolation_source <- "Nasopharyngeal swab"
  CHOA_Original_attributes[os,]$host_age <- paste0(as.character(round(this_sample$Age_mo/12, 2), " years"))
  CHOA_Original_attributes[os,]$host_health_state <- this_sample$PMCA
  CHOA_Original_attributes[os,]$host_sex <- this_sample$Gender
  CHOA_Original_attributes[os,]$passage_history <- "Direct from clinical specimen"

}
View(CHOA_Original_attributes)

# Populate the CHOA_Replicate_attributes table
# Loop over the original and replicate sample entries to populate the CHOA_attributes tibble
rs <- 1
for (rs in 1:length(Replicate_CHOA_Samples)){ # where cs is index for the CHOA_sample
  sample_tmp <- Replicate_CHOA_Samples[rs]; sample_tmp
  
  # Get the subject information for the current sample by extracting the first component of the name
  sample_name_tmp <- unlist(strsplit(x = sample_tmp,"_"))[1]
  sample_name_tmp
  
  # Extract the relevant entry from combined_data dataframe
  this_sample <- combined_data[combined_data$SubjectID == sample_name_tmp,]
  CHOA_Replicate_attributes[rs,]$replicate <- 2
  CHOA_Replicate_attributes[rs,]$sample_name <- sample_tmp
  CHOA_Replicate_attributes[rs,]$sample_title <- sample_name_tmp
  CHOA_Replicate_attributes[rs,]$bioproject_accession <- ""
  CHOA_Replicate_attributes[rs,]$organism <- "Influenza A virus"
  CHOA_Replicate_attributes[rs,]$strain <- this_sample$Strain
  CHOA_Replicate_attributes[rs,]$subtype <- this_sample$Strain
  CHOA_Replicate_attributes[rs,]$isolate <- paste0(this_sample$SubjectID,"_2")
  CHOA_Replicate_attributes[rs,]$collected_by <- "Children's Hopsital of Philadelphia"
  CHOA_Replicate_attributes[rs,]$collection_date <- this_sample$collection_date
  CHOA_Replicate_attributes[rs,]$geo_loc_name <- "USA: Pennsylvania: Philadelphia"
  CHOA_Replicate_attributes[rs,]$lat_lon <- '39.9487 N 75.1939 W'
  CHOA_Replicate_attributes[rs,]$host <- "Homo sapiens"
  CHOA_Replicate_attributes[rs,]$host_disease <- "Infuenza"
  CHOA_Replicate_attributes[rs,]$isolation_source <- "Nasopharyngeal swab"
  CHOA_Replicate_attributes[rs,]$host_age <- paste0(as.character(round(this_sample$Age_mo/12, 2), " years"))
  CHOA_Replicate_attributes[rs,]$host_health_state <- this_sample$PMCA
  CHOA_Replicate_attributes[rs,]$host_sex <- this_sample$Gender
  CHOA_Replicate_attributes[rs,]$passage_history <- "Direct from clinical specimen"
  
}
View(CHOA_Replicate_attributes)

# Combine the original and replicate attributes tables
CHOA_attributes <- rbind(CHOA_Original_attributes,CHOA_Replicate_attributes)

# Save a version for the sra sample attributes table
CHOA_sra_attributes <- CHOA_attributes

# Remove the replicate column
# Remove the 'replicate' column
CHOA_sra_attributes <- CHOA_sra_attributes %>% 
  select(-replicate)


# Save then CHOA_sra_attributes table as a csv file. The entries can the be pasted into the excel template for uploading to SRA.
write.csv(x = CHOA_sra_attributes,file = file.path(SRA_info_path,"CHOA_BioSample_Attributes.csv"))
View(CHOA_attributes)

# Generate SRA sample metadata table -----------------------------------------------------

# Generate a tibble to store the sequence metadata
sra_metadata <- tibble(
  sample_name = character(),
  library_ID = character(),
  title = character(),
  library_strategy = character(),
  library_source = character(),
  library_selection = character(),
  library_layout = character(),
  platform = character(),
  instrument_model = character(),
  design_description = character(),
  filetype = character(),
  filename = character(),
  filename2 = character()
)

# Loop over the CHOA_samples entries to populate the sample metadata tibble
ms <- 1
for (ms in 1:nrow(CHOA_attributes)){ # where cs is index for the CHOA_sample

  sample_tmp <- CHOA_attributes[ms,]$sample_name; sample_tmp
  
  # Get the subject information for the current sample by extracting the first component of the name
  sample_name_tmp <- unlist(strsplit(x = sample_tmp,"_"))[1]
  sample_name_tmp
  
  # Extract the relevant entry from combined_data dataframe
  this_sample <- combined_data[combined_data$SubjectID == sample_name_tmp,]
  this_sample
  
  # Determine if this is the original sample or the replicate by determining if
  # the current sample already exists in the sra_metadata table then assign the
  # sample name of the original vs replicate
  if (CHOA_attributes[ms,]$replicate == 1){
    # This sample is the original
    sra_metadata[ms,]$sample_name <- this_sample$Subject_Original
  } else {
    sra_metadata[ms,]$sample_name <- this_sample$Subject_Replicate
  }
  
  # We use the GetFastq_name function to obtain both the sequencing date and the
  # full name of the sample for matching the fastq file
  filename_components <- GetFastq_name(sample_tmp)
  
  runDate <- filename_components[1]
  full_name <- filename_components[2]
  filename_components
  # Assign convert the string representing the sequencing date to a proper date 
  sample_date <- convert_to_date(runDate)
  
  # We set the library ID based on the study and the date the samples were run 
  sra_metadata[ms,]$library_ID <- paste0("IAV_WGS_",as.character(sample_date),"_",full_name)
  sra_metadata[ms,]$library_ID
  # We set the title for the experiment based on the method for dataset creation
  sra_metadata[ms,]$title <- "WGS of influenza A virus: clinical nasopharyngeal swabs of naturally infected pediatric patients"
  
  # We set the library parameters (strategy, source, selection, layout and platform)
  sra_metadata[ms,]$library_strategy <- "WGS"
  sra_metadata[ms,]$library_source <- "VIRAL RNA"
  sra_metadata[ms,]$library_selection <- "PCR"
  sra_metadata[ms,]$library_layout <- "paired"
  sra_metadata[ms,]$platform <- "ILLUMINA"
  
  # There were two instrument models used for this experiment (sequencer was
  # replaced mid-way during the study). Plates sequenced before 01/01/2023 used
  # the nextseq 550 whereas plates sequenced afterwards used the NextSeq 2000.
  # We use the date of sequencing to determine the Instrument model
  if (sample_date < as.Date("2023-01-01")){
    sra_metadata[ms,]$instrument_model <- "NextSeq 550"
  } else {
    sra_metadata[ms,]$instrument_model <- "NextSeq 2000"
  }
  
  # Here we add a brief description of the sequencing method
  sra_metadata[ms,]$design_description <- "RNA was extracted from remnant banked nasopharyngeal swabs from patients with symptoms of influenza virus. Influenza viral RNA was amplified and converted to cDNA using universal influenza A primers. Libraries were prepared using Illumina DNA Prep kits without enrichment."
  
  # We will upload the fastq files for the sequences and here set the type and file names.
  sra_metadata[ms,]$filetype <- "fastq"
  
  # Generate the fastq names from the full_name for the forward (R1) and reverse (R2) read pairs
  sra_metadata[ms,]$filename <- paste0(full_name,"_R1_001.fastq.gz")
  sra_metadata[ms,]$filename2 <- paste0(full_name,"_R2_001.fastq.gz")
  
  # # The reference used is based on the influenza subtype H1N1pdm of H3N2. We use this information to define the determine the reference accession numbers
  # if (this_sample$Strain == "H3N2"){
  #   sra_metadata[ms,] <- Ref_H3N2
  # } else if (this_sample$Strain == "H1N1"){
  #   sra_metadata[ms,] <- Ref_H1N1pdm
  # } else{
  #   stop(print("The strain assignement is not H3N2 or H1N1pdm. Please asssess further"))
  # }
  
    
}

View(sra_metadata)

# There was an error in the initial fastq file transfers for three samples:
# CHOA-009_D1V1_S67, CHOA-034_D2V1_S58 and CHOA-196_D2V1_S64. CHOA-009_D1V1_S67
# had been incorrectly named in the fastq files CHOA-090_D2V1_S67_R1/R2 (this
# was corrected manually for the pipeline output and verified by comparison with
# the sample sheet). Following verification, the fastq file name was corrected
# prior to upload with SRA. The CHOA-196_D2V1_S64 was incorrectly named
# CHOA-196-D2V1_S64 and also was manually corrected for upload to SRA. The fastq
# files for CHOA-034_D2V1_S58 were missing from the file containing the fastq
# files for the 050123 run. Interestingly, the pipeline output data was also
# missing for this sample. To address this, I uploadd the missing files from an
# alternate saved version of the FluPipeline processing for 050123.

# Save then sra_metadata table as a csv file. The entries can the be pasted into the excel template for uploading to SRA.
write.csv(x = sra_metadata,file = file.path(SRA_info_path,"CHOA_sra_metadata.csv"))






