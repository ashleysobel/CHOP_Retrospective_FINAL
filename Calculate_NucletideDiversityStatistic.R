# This script specifically calculates the nucleotide diversity for retrospective
# CHOP data, contrasting with other scripts by focusing on this genetic
# diversity metric. It processes and analyzes genetic sequences, assessing the
# diversity within specific segments and samples, and visualizes these
# calculations.
# Author: Ashley Sobel Leonard
# Date: 5/1/2024


# Prepare environment ----------------------------------------------------- 
# Clear the workspace to ensure a clean environment for running the script
rm(list = ls())

# Get the current working directory to use as the base for all relative paths
R_path <- getwd()

# Define paths to where data structures are stored/should be written
Figure_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","Figures")
Compiled_OutPut_path <- file.path(R_path,"PostProcessing_OutPut","CompiledOutPut")

# Load packages --------
library(tidyverse)
library(gridExtra)
library(data.table)

# Generate .bib citation file for statistical packages used in analysis
# List of packages for which you want citations
package_list <- c("tidyverse","gridExtra","data.table")

# Initialize an empty character vector to hold citations
all_citations <- character(0)

# Loop through each package and get its citation in BibTeX format
for (pkg in package_list) {
  citation_info <- utils::toBibtex(citation(pkg))
  all_citations <- c(all_citations, citation_info)
}

# Combine all citations into a single string
all_citations_text <- paste(all_citations, collapse = "\n\n")

# Write the combined citations to a .bib file
writeLines(all_citations_text, con = file.path(Compiled_OutPut_path,"NucleotideDiversity_citations.bib"))

# Load reference data -------
# Load merged data (the sample names, clinical info, and summary information)
merged_data <- read.csv(paste0(Compiled_OutPut_path,"/FINAL_CHOA_merged_data.csv"))

# Read in the csv file that contains selected samples
CHOA_Selected_Samples <-
  read.csv(paste0(Compiled_OutPut_path,"/CHOA_Selected_Samples_Revised.csv"))

# Remove mixed infectons from existing dataframes
CHOA_Selected_Samples <- CHOA_Selected_Samples %>%
  filter(!ID %in% c("CHOA-101","CHOA-117"))

# Define Functions -----

# Define a function to process pileup files and calculate nucleotide diversity statistics
get_pileup <- function(sample){
  # Load and preprocess pileup data for sample, calculating nucleotide positions and diversity
  sample_name <- sample$Sample 
  seq_date <- sample$SeqDate
  
  # Construct the path to the pileup data files.
  pileup_path <- file.path(R_path,"Pipeline_Output_Files",paste0(sample$SeqDate,"_PipelineOutput"),"Pileup")
  
  # Search for the appropriate pileup file in the given directory.
  pileup_file <- list.files(path=pileup_path,pattern = sample$Sample,full.names = TRUE)
  
  # Load the pileup file as a data frame.
  working_sample_pileup <- read.csv(pileup_file)
  working_sample_pileup$X <- NULL
  working_sample_pileup <- data.table::data.table(working_sample_pileup)
  
  # Create a new column 'seqnum' that contains the last character of each entry in the 'seqnames' colum
  working_sample_pileup[, seqnum := substr(seqnames, nchar(seqnames), nchar(seqnames))]
  
  # Partition the pileup data into 8 separate data tables based on 'seqnum'.
  PB2_working_pileup <- working_sample_pileup[seqnum == "1"]
  PB1_working_pileup <- working_sample_pileup[seqnum == "2"]
  PA_working_pileup <- working_sample_pileup[seqnum == "3"]
  HA_working_pileup <- working_sample_pileup[seqnum == "4"]
  NP_working_pileup <- working_sample_pileup[seqnum == "5"]
  NA_working_pileup <- working_sample_pileup[seqnum == "6"]
  MP_working_pileup <- working_sample_pileup[seqnum == "7"]
  NS_working_pileup <- working_sample_pileup[seqnum == "8"]
  
  # Combine the separated data tables into a list for easy processing.
  pileup_list <- list(PB2_working_pileup, PB1_working_pileup, PA_working_pileup, 
                      HA_working_pileup, NP_working_pileup, NA_working_pileup, 
                      MP_working_pileup, NS_working_pileup)
  
  # Apply a function over each data table in the list that reshapes the data,
  # replaces non-standard nucleotides with "N", and calculates the depth at each
  # position.
  pileup_list_reshaped <- lapply(pileup_list, function(pileup_df) {
    pileup_df %>%
      dplyr::mutate(nucleotide = dplyr::case_when(
        nucleotide %in% c("A", "T", "G", "C") ~ nucleotide,
        TRUE ~ "N"
      )) %>%
      dplyr::group_by(pos, nucleotide) %>%
      dplyr::summarise(count = sum(count), .groups = "drop") %>%
      pivot_wider(names_from = nucleotide, values_from = count, values_fill = 0) %>%
      dplyr::mutate(Depth = A + T + G + C)
  })
  # Return the reshaped pileup data.
  return(pileup_list_reshaped)
}


# Set run options ------ 
min_coverage <- 100 # We set minimum coverage at 100 reads 

# Calculate nucleotide diversity for each sample ------
sub <- 1
# Initialize dataframe to store results
Diversity_info <- data.frame(Sample=rep(NA,length(CHOA_Selected_Samples$Sample)),Length=rep(NA,length(CHOA_Selected_Samples$Sample)),Length_adj=rep(NA,length(CHOA_Selected_Samples$Sample)),Pi=rep(NA,length(CHOA_Selected_Samples$Sample)))

# Iterate through selected samples to calculate nucleotide diversity per segment
for (sub in 1:length(CHOA_Selected_Samples$Sample)){
  sample_tmp <- CHOA_Selected_Samples[sub,]
  
  # Load the pileup file for the current sample 
  pileup_tmp <- get_pileup(sample_tmp)
  
  # Generate dataframe to hold the nucleotide diversity information 
  D_info <- data.frame(Dseg=rep(NA, length(pileup_tmp)), L_seg=rep(NA, length(pileup_tmp)),L_adj_seg=rep(NA, length(pileup_tmp)))
  
  # Calculate the nucleotide diversity for each segment 
  seg <- 1
  for (seg in 1:length(pileup_tmp)){
    seg_tmp <- pileup_tmp[[seg]]
    
    # Calculate Dl, considering the minimum coverage requirement
    seg_tibble <- seg_tmp %>%
      mutate(Dl = ifelse(Depth > min_coverage, 
                         (Depth * (Depth - 1) - (T * (T - 1) + C * (C - 1) + A * (A - 1) + G * (G - 1))) / (Depth * (Depth - 1)), 
                         NaN))
    
    # Sum the Dl values that are not NaN. Also determine the segment length and
    # the adjusted length (L_seg_adj) which considers only those positions with
    # sufficient depth
    D_seg <- sum(seg_tibble$Dl, na.rm = TRUE)
    L_seg <- nrow(seg_tibble)
    L_seg_adj <- sum(!is.nan(seg_tibble$Dl))

    # Store values in D_info
    D_info[seg, ] <- c(D_seg, L_seg,L_seg_adj)
  }
  
  
  # Store values in D_info
  Length_tmp <- sum(D_info$L_seg)
  Length_adj_tmp <- sum(D_info$L_adj_seg)
  Total_D <- sum(D_info$Dseg)
  Pi_tmp <- Total_D/Length_adj_tmp
  print(c(sample_tmp$Sample,Pi_tmp))
  Diversity_info[sub,] <- c(sample_tmp$Sample,Length_tmp,Length_adj_tmp,Pi_tmp) 
}

# Reformat Diversity_info into a more useful framework and save the results
Diversity_info <- Diversity_info %>%
  mutate(SubjectID = str_split(Sample, "_", simplify = TRUE)[,1])
Diversity_info$Pi <- as.double(Diversity_info$Pi)
write.csv(x = Diversity_info, file = file.path(Compiled_OutPut_path,"NucleotideDiversity.csv"))

# Calulcate the avrerage diversity and save the results
Diversity_avg <- Diversity_info %>%
  group_by(SubjectID) %>%
  summarise(Avg_Pi = mean(Pi, na.rm = TRUE))
write.csv(x = Diversity_avg, file = file.path(Compiled_OutPut_path,"AverageNucleotideDiversity.csv"))

# Order SubjectIDs by average diversity
Diversity_avg <- Diversity_avg %>%
  arrange(Avg_Pi)

# Convert SubjectID to a factor ordered by Avg_Pi
Diversity_avg$SubjectID <- factor(Diversity_avg$SubjectID, levels = Diversity_avg$SubjectID)
mean_diversity <- mean(Diversity_avg$Avg_Pi)

# Plot with percentiles, y-axis limits, and updated font sizes 
Diversity_plot <- ggplot(Diversity_avg, aes(x = SubjectID, y = Avg_Pi)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  geom_hline(yintercept = mean_diversity, color = "darkgrey", linetype = "dashed", size = 0.5) + # Updated line
  theme_minimal(base_size = 8) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "darkgrey", fill = NA, size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 10)  # Adjust y-axis label font size
  ) +
  scale_y_continuous(limits = c(0, 0.004)) +  # Set y-axis limits
  labs(x = "SubjectID", y = expression("Average Nucleotide Diversity (" * pi * ")"))

# Save the output
Diversity_plot_name <- paste0("20424-NucleotideDiversity_plot.eps")
ggplot2::ggsave(file.path(Figure_OutPut_path,Diversity_plot_name), plot = Diversity_plot, width = 10.5, height = 5, device = "eps")

# Generate diversity association plots -----

# Load the diverity data
diversity_info <- read.csv(file.path(Compiled_OutPut_path,"NucleotideDiversity.csv"))

# Merge 'diversity_info' with 'merged_data' on 'SubjectID' to incorporate clinical metadata
combined_data <- merge(merged_data, diversity_info, by.x = "SubjectID", by.y = "SubjectID")

# Create a new association_data dataframe using nucleotide diversity (Pi)
association_data <- data.frame(
  SubjectID = combined_data$SubjectID,
  CT = combined_data$CT,
  Age_yo = combined_data$Age_mo / 12,
  Vaccine = combined_data$Flu_vaccine,
  Symptom_dates = combined_data$Symptom_days,
  PMCA = combined_data$PMCA,
  Setting = combined_data$Setting,
  Oseltamivir = combined_data$Prior_tamiflu,
  Admit = combined_data$Admit,
  Pi = combined_data$Pi # Use nucleotide diversity here
)


# Replace "Unknown" in "Vaccine" column with "No"
association_data$Vaccine[association_data$Vaccine == "Unknown"] <- "No"

# Fit a Negative Binomial GLM for nucleotide diversity and summarize the output
glm_model_Pi <- MASS::glm.nb(Pi ~ Vaccine + PMCA + Setting + Admit, data = association_data)
summary(glm_model_Pi)

# Extract p-values from the model and apply Bonferroni correction
pvals_Pi <- summary(glm_model_Pi)$coefficients[, "Pr(>|z|)"]
adjusted_pvals_Pi <- p.adjust(pvals_Pi, method = "bonferroni")

# Calculate r^2 for Pi vs Age_yo
model_Pi_Age <- lm(Pi ~ Age_yo, data = association_data)
summary_Pi_Age <- summary(model_Pi_Age)
r2_Pi_Age <- summary_Pi_Age$r.squared
print(paste0("R-squared for Nucleotide Diversity vs Age_yo: ", r2_Pi_Age))

# Adjusted scatter plot for Pi vs Age with filled grey circles
age_Pi_plot <- ggplot(association_data, aes(x = Age_yo, y = Pi)) +
  geom_point(color = "darkgrey", fill = "darkgrey", shape = 21, size = 1.25) +
  geom_smooth(method = "lm", linetype = "solid", color = "black", se = FALSE) +
  labs(x = "Age (years)", y = "Nucleotide Diversity (Pi)") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=12, family="Helvetica"), axis.title = element_text(size=12))

# Show the plot
print(age_Pi_plot)

# Save the plot
ggsave(filename = "Age_Pi_plot.eps", plot = age_Pi_plot, width = 5, height = 3.5, device = "eps")

# Adjusted boxplot for Pi vs Vaccine Status with custom color
Pi_vaccine_plot <- ggplot(association_data, aes(x = Vaccine, y = Pi)) +
  geom_boxplot(fill = "darkgrey") +
  labs(x = "Vaccine Status", y = "Nucleotide Diversity (Pi)") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=12, family="Helvetica"), axis.title = element_text(size=12))

# Show the plot
print(Pi_vaccine_plot)

# Save the plot
ggsave(filename = "Pi_Vaccine_Status_plot.eps", plot = Pi_vaccine_plot, width = 3.5, height = 3.5, device = "eps")

# Symptom Days vs Nucleotide Diversity Plot
SymptomDays_Pi_plot <- ggplot(association_data, aes(x = Symptom_dates, y = Pi)) +
  geom_point(color = "darkgrey", fill = "darkgrey", shape = 21, size = 1.25) +
  geom_smooth(method = "lm", linetype = "solid", color = "black", se = FALSE) +
  labs(x = "Days of Symptoms", y = "Nucleotide Diversity (Pi)") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=12, family="Helvetica"), axis.title = element_text(size=12))

# Show the plot
print(SymptomDays_Pi_plot)

# Save the plot
ggsave(filename = "SymptomDays_Pi_plot.eps", plot = SymptomDays_Pi_plot, width = 5, height = 3.5, device = "eps")

# PMCA vs Nucleotide Diversity Plot
# Ensure PMCA is treated as a factor for plotting purposes
association_data$PMCA <- as.factor(association_data$PMCA)

Pi_PMCA_plot <- ggplot(association_data, aes(x = PMCA, y = Pi)) +
  geom_boxplot(fill = "darkgrey") +
  labs(x = "PMCA Status", y = "Nucleotide Diversity (Pi)") +
  theme_bw() + 
  theme(legend.position="none", text = element_text(size=12, family="Helvetica"), axis.title = element_text(size=12))

# Show the plot
print(Pi_PMCA_plot)

# Save the plot
ggsave(filename = "Pi_PMCA_plot.eps", plot = Pi_PMCA_plot, width = 3.5, height = 3.5, device = "eps")

# Generate the final plot containing all asscociations, which will become Figure S3
# Arrange the final plots adjacent to each other
final_plot_pi <- 
  grid.arrange(
    age_Pi_plot, Pi_vaccine_plot, SymptomDays_Pi_plot, Pi_PMCA_plot,
    ncol = 2, 
    nrow = 2, 
    widths = c(1.5, 1) # Specifies the relative widths of the columns
  )
print(final_plot_pi)

# Save the plot
ggplot2::ggsave(file.path(Figure_OutPut_path,"final_pi_plot.eps"), plot = final_plot_pi, width = 8, height = 8, device = "eps")




