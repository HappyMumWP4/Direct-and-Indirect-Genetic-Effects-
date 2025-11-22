################################################################################
################## Within-study PRSes correlation analysis #####################
################################################################################

# Install and load necessary packages
packages <- c("tidyverse", "Hmisc", "openxlsx")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(tidyverse)
# library(readxl)
# library(openxlsx)
library(Hmisc)
# library(ggplot2)

# # Set working directory (update to your actual path)
# setwd("C:/path/to/your/files")  # Change this to your directory

# List of CSV files for each study (update with actual file names)
# You can use the output winsored_standardized data csv file written in line 453 of the original analysis R script
# For example the original analysis R script "Genetic nurture_analytical_models_CBCL1-5"
files <- list("PREDO_winsored_standardized_trio_dataset_cbcl1_5_fa.csv",
              "PREDO_winsored_standardized_trio_dataset_cbcl1_5_mo.csv") # change to your cohort's raw data files

# # List of xlsx files for each study (update with actual file names)
# files <- list("PREDO_Trio dataset_CBCL7-12y_fa.xlsx", 
#               "PREDO_Trio dataset_CBCL7-12y_mo.xlsx") # change to your cohort's raw data files

# Define output directory
output_dir <- "PRScmf_within-study_corr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each study file
for (file in files) {
  
  # Check if file exists before reading
  if (!file.exists(file)) {
    cat("❌ File not found:", file, "\n")
    next  # Skip this file and continue the loop
  }
  
  # Read the data
  # data <- read_csv(file)  
  data <- read_xlsx(file)
  
  # Extract all predictor column names
  predictor_columns <- grep("^prs_[cfm]_qc_", colnames(data), value = TRUE)
  
  # Extract unique trait names (e.g., ADHD, ALC, AN, etc.)
  trait_names <- unique(gsub("^prs_[cfm]_qc_|_trio$", "", predictor_columns))
  
  # Initialize an empty list to store all correlation matrices
  all_correlations <- list()
  
  # Loop through each trait
  for (trait in trait_names) {
    
    # Select predictors belonging to this trait
    subset_cols <- grep(paste0("_qc_", trait, "_trio$"), predictor_columns, value = TRUE)
    
    # If fewer than two predictors exist, skip (correlation needs at least 2 variables)
    if (length(subset_cols) < 2) next
    
    # Extract relevant predictor data
    predictor_data <- data %>% select(all_of(subset_cols))
    
    # Compute correlation matrix with rcorr
    cor_results <- rcorr(as.matrix(predictor_data), type = "pearson")
    
    # Extract r, p, and n matrices
    r_matrix <- cor_results$r
    p_matrix <- cor_results$P
    n_matrix <- cor_results$n
    
    # Convert each to long format
    r_df <- as.data.frame(as.table(r_matrix))
    p_df <- as.data.frame(as.table(p_matrix))
    n_df <- as.data.frame(as.table(n_matrix))
    
    # Combine them
    prs_cor_full_df <- r_df %>%
      rename(Predictor1 = Var1, Predictor2 = Var2, Correlation = Freq) %>%
      mutate(p_value = p_df$Freq, N = n_df$Freq)
    
    # Append to list
    all_correlations[[trait]] <- prs_cor_full_df
    
    # # Plotting correlation heatmap
    # r_plot_df <- result_df %>%
    #   filter(Var1 != Var2) %>%
    #   mutate(pair = paste(Var1, Var2, sep = " vs "))
    # 
    # p <- ggplot(r_plot_df, aes(x = pair, y = r)) +
    #   geom_bar(stat = "identity", fill = "#377EB8") +
    #   ylim(-1, 1) +
    #   theme_minimal() +
    #   ggtitle(paste("PRS Correlation -", trait)) +
    #   ylab("Pearson r") +
    #   xlab("Pair") +
    #   geom_text(aes(label = paste0("p=", signif(p, 2))), vjust = -0.5, size = 3)
    # 
    # ggsave(paste0("prs_correlation_", trait, ".png"), plot = p, width = 6, height = 4)
  }
  
  # Combine all correlations into a single dataframe
  prs_cor_final_df <- bind_rows(all_correlations)
  
  # Print the correlation matrix for this study
  cat("\nCorrelation matrix for", file, ":\n")
  print(prs_cor_final_df)
  
  # Generate output filename using the last 2 sectors of the input file name
  file_parts <- unlist(strsplit(tools::file_path_sans_ext(basename(file)), "_"))
  output_filename <- paste0("PRScmf_corr_", paste(tail(file_parts, 2), collapse = "_"), ".csv")
  csv_file <- file.path(output_dir, output_filename)
  
  # Save the correlation matrix as a CSV file with the new name
  write.csv(prs_cor_final_df, file = csv_file, row.names = FALSE)
  
  cat("✅ Saved correlation matrix to", csv_file, "\n")
}


################################################################################
######### Within-study psychiatric PRS predictors correlation analysis #########
################################################################################

# Install and load necessary packages
packages <- c("tidyverse", "Hmisc", "openxlsx")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load required libraries
library(tidyverse)
library(readxl)
library(openxlsx)
library(Hmisc)

# # Set working directory (update to your actual path)
# setwd("C:/path/to/your/files")  # Change this to your directory

# # List of CSV files for each study (update with actual file names)
# files <- list("PREDO_Trio dataset_CBCL7-12y_fa.csv",
#               "PREDO_Trio dataset_CBCL7-12y_mo.csv") # change to your cohort's raw data files

# List of xlsx files for each study (update with actual file names)
files <- list("PREDO_Trio dataset_CBCL7-12y_fa.xlsx",
              "PREDO_Trio dataset_CBCL7-12y_mo.xlsx") # change to your cohort's raw data files

# Define output directory
output_dir <- "Psy-PRSes_within-study_corr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each study file
for (file in files) {
  
  # Check if file exists before reading
  if (!file.exists(file)) {
    cat("File not found:", file, "\n")
    next  # Skip this file and continue the loop
  }
  
  # Read the data
  # data <- read_csv(file)  
  data <- read_xlsx(file)
  
  # Select only the psychiatric prs variables (replace with actual column names)
  predictor_data <- data %>% select(prs_c_qc_ADHD_trio,	prs_c_qc_ALC_trio,	prs_c_qc_AN_trio,	prs_c_qc_ANX_trio,	prs_c_qc_ASD_trio,	prs_c_qc_BD_trio,	prs_c_qc_CDG_trio,	prs_c_qc_CIG_trio,	prs_c_qc_EDU_trio,	prs_c_qc_INSOM_trio,	prs_c_qc_MDD_trio,	prs_c_qc_NEUROT_trio,	prs_c_qc_OCD_trio,	prs_c_qc_PPD_trio,	prs_c_qc_PTSD_trio,	prs_c_qc_SCZ_trio,
                                    prs_f_qc_ADHD_trio,	prs_f_qc_ALC_trio,	prs_f_qc_AN_trio,	prs_f_qc_ANX_trio,	prs_f_qc_ASD_trio,	prs_f_qc_BD_trio,	prs_f_qc_CDG_trio,	prs_f_qc_CIG_trio,	prs_f_qc_EDU_trio,	prs_f_qc_INSOM_trio,	prs_f_qc_MDD_trio,	prs_f_qc_NEUROT_trio,	prs_f_qc_OCD_trio,	prs_f_qc_PPD_trio,	prs_f_qc_PTSD_trio,	prs_f_qc_SCZ_trio,		
                                    prs_m_qc_ADHD_trio,	prs_m_qc_ALC_trio,	prs_m_qc_AN_trio,	prs_m_qc_ANX_trio,	prs_m_qc_ASD_trio,	prs_m_qc_BD_trio,	prs_m_qc_CDG_trio,	prs_m_qc_CIG_trio,	prs_m_qc_EDU_trio,	prs_m_qc_INSOM_trio,	prs_m_qc_MDD_trio,	prs_m_qc_NEUROT_trio,	prs_m_qc_OCD_trio,	prs_m_qc_PPD_trio,	prs_m_qc_PTSD_trio,	prs_m_qc_SCZ_trio,	
                                    prs_nt_qc_ADHD_trio,	prs_nt_qc_ALC_trio,	prs_nt_qc_AN_trio,	prs_nt_qc_ANX_trio,	prs_nt_qc_ASD_trio,	prs_nt_qc_BD_trio,	prs_nt_qc_CDG_trio,	prs_nt_qc_CIG_trio,	prs_nt_qc_EDU_trio,	prs_nt_qc_INSOM_trio,	prs_nt_qc_MDD_trio,	prs_nt_qc_NEUROT_trio,	prs_nt_qc_OCD_trio,	prs_nt_qc_PPD_trio,	prs_nt_qc_PTSD_trio,	prs_nt_qc_SCZ_trio
  )
  
  # Compute correlation matrix with rcorr
  cor_results <- rcorr(as.matrix(predictor_data), type = "pearson")
  
  # Extract r, p, and n matrices
  r_matrix <- cor_results$r
  p_matrix <- cor_results$P
  n_matrix <- cor_results$n
  
  # Convert each to long format
  r_df <- as.data.frame(as.table(r_matrix))
  p_df <- as.data.frame(as.table(p_matrix))
  n_df <- as.data.frame(as.table(n_matrix))
  
  # Combine them
  psychprs_cor_full_df <- r_df %>%
    rename(Predictor1 = Var1, Predictor2 = Var2, Correlation = Freq) %>%
    mutate(p_value = p_df$Freq, N = n_df$Freq)
  
  # Print the correlation matrix for this study
  cat("\nCorrelation matrix for", file, ":\n")
  print(psychprs_cor_full_df)
  
  # Generate output filename using the last 2 sectors of the input file name
  file_parts <- unlist(strsplit(tools::file_path_sans_ext(basename(file)), "_"))
  output_filename <- paste0("Psych-PRSes_corr_", paste(tail(file_parts, 2), collapse = "_"), ".csv")
  csv_file <- file.path(output_dir, output_filename)
  
  # Save the correlation matrix as a CSV file with the new name
  write.csv(psychprs_cor_full_df, file = csv_file, row.names = FALSE)
  
  cat("✅ Saved correlation matrix to", csv_file, "\n")
}


################################################################################
################# Within-study Outcome correlation analysis ####################
################################################################################

# Install and load necessary packages
packages <- c("tidyverse", "Hmisc", "openxlsx")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load required libraries
library(tidyverse)
library(readxl)
library(openxlsx)
library(Hmisc)

# # Set working directory (update to your actual path)
# setwd("C:/path/to/your/files")  # Change this to your directory

# # List of CSV files for each study (update with actual file names)
# files <- list("PREDO_Trio dataset_CBCL7-12y_fa.csv",
#               "PREDO_Trio dataset_CBCL7-12y_mo.csv") # change to your cohort's raw data files

# List of xlsx files for each study (update with actual file names)
files <- list("PREDO_Trio dataset_CBCL7-12y_fa.xlsx",
              "PREDO_Trio dataset_CBCL7-12y_mo.xlsx") # change to your cohort's raw data files

# Define output directory
output_dir <- "Outcomes_within-study_corr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop through each study file
for (file in files) {
  
  # Check if file exists before reading
  if (!file.exists(file)) {
    cat("File not found:", file, "\n")
    next  # Skip this file and continue the loop
  }
  
  # Read the data
  # data <- read_csv(file)  
  data <- read_xlsx(file)
  
  # Select only the three outcome variables (replace with actual column names)
  outcome_data <- data %>% select(h_c_extern_6_18y, h_c_intern_6_18y, h_c_total_6_18y)
  
  # Compute correlation matrix with rcorr
  cor_results <- rcorr(as.matrix(outcome_data), type = "pearson")
  
  # Extract r, p, and n matrices
  r_matrix <- cor_results$r
  p_matrix <- cor_results$P
  n_matrix <- cor_results$n
  
  # Convert each to long format
  r_df <- as.data.frame(as.table(r_matrix))
  p_df <- as.data.frame(as.table(p_matrix))
  n_df <- as.data.frame(as.table(n_matrix))
  
  # Combine them
  outcome_cor_full_df <- r_df %>%
    rename(Outcome1 = Var1, Outcome2 = Var2, Correlation = Freq) %>%
    mutate(p_value = p_df$Freq, N = n_df$Freq)
  
  # Print the correlation matrix for this study
  cat("\nCorrelation matrix for", file, ":\n")
  print(outcome_cor_full_df)
  
  # Generate output filename using the last 2 sectors of the input file name
  file_parts <- unlist(strsplit(tools::file_path_sans_ext(basename(file)), "_"))
  output_filename <- paste0("outcomes_wn_study_corr_", paste(tail(file_parts, 2), collapse = "_"), ".csv")
  csv_file <- file.path(output_dir, output_filename)
  
  # Save the correlation matrix as a CSV file with the new name
  write.csv(outcome_cor_full_df, file = csv_file, row.names = FALSE)
  
  cat("✅ Saved correlation matrix to", csv_file, "\n")
}


################################################################################
################ Between-study outcome correlation analysis ####################
################################################################################

# Install necessary packages if not already installed
packages <- c("tidyverse", "Hmisc", "openxlsx")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load required libraries
library(tidyverse)
library(readxl)
library(openxlsx)
library(Hmisc)

# Define output directory
output_dir <- "Outcomes_btw_mo-fa_corr"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# # List of CSV files for each study (update with actual file names)
# files <- list("PREDO_Trio dataset_CBCL7-12y_fa.csv",
#               "PREDO_Trio dataset_CBCL7-12y_mo.csv") # change to your cohort's raw data files

# List of xlsx files for each study (update with actual file names)
files <- list("PREDO_Trio dataset_CBCL7-12y_fa.xlsx",
              "PREDO_Trio dataset_CBCL7-12y_mo.xlsx") # change to your cohort's raw data files

# Read the datasets
# study1_data <- read_csv(files[[1]]) %>% mutate(Study = "Study_fa")
# study2_data <- read_csv(files[[2]]) %>% mutate(Study = "Study_mo")
study1_data <- read_xlsx(files[[1]]) %>% mutate(Study = "Study_fa")
study2_data <- read_xlsx(files[[2]]) %>% mutate(Study = "Study_mo")

# Merge datasets by a common ID (assuming ID exists)
# If no common ID, correlations cannot be computed this way
merged_data <- inner_join(study1_data, study2_data, by = "final_id", suffix = c("_fa", "_mo"))

# Compute correlations
cor_ext_Study1_vs_ext_Study2 <- cor(merged_data$h_c_extern_6_18y_fa, merged_data$h_c_extern_6_18y_mo, use = "complete.obs")
cor_int_Study1_vs_int_Study2 <- cor(merged_data$h_c_intern_6_18y_fa, merged_data$h_c_intern_6_18y_mo, use = "complete.obs")
cor_ext_Study1_vs_int_Study2 <- cor(merged_data$h_c_extern_6_18y_fa, merged_data$h_c_intern_6_18y_mo, use = "complete.obs")
cor_int_Study1_vs_ext_Study2 <- cor(merged_data$h_c_intern_6_18y_fa, merged_data$h_c_extern_6_18y_mo, use = "complete.obs")
cor_tot_Study1_vs_tot_Study2 <- cor(merged_data$h_c_total_6_18y_fa, merged_data$h_c_total_6_18y_mo, use = "complete.obs")

# Print results
cat("Between-study correlation for Outcome1:", cor_ext_Study1_vs_ext_Study2, "\n")
cat("Between-study correlation for Outcome2:", cor_int_Study1_vs_int_Study2, "\n")
cat("Cross-outcome between-study correlation (Outcome1 in Study_fa vs Outcome2 in Study_mo):", cor_ext_Study1_vs_int_Study2, "\n")
cat("Cross-outcome between-study correlation (Outcome2 in Study_fa vs Outcome1 in Study_mo):", cor_int_Study1_vs_ext_Study2, "\n")
cat("Between-study correlation for Outcome3:", cor_tot_Study1_vs_tot_Study2, "\n")

# Save results as a CSV file
cor_results <- data.frame(
  Comparison = c("Ext_Study_fa_vs_Ext_Study_mo", "Int_Study_fa_vs_Int_Study_mo", "Ext_Study_fa_vs_Int_Study_mo", 
                 "Int_Study_fa_vs_Ext_Study_mo", "Tot_Study_fa_vs_Tot_Study_mo"),
  Correlation = c(cor_ext_Study1_vs_ext_Study2, cor_int_Study1_vs_int_Study2, cor_ext_Study1_vs_int_Study2, 
                  cor_int_Study1_vs_ext_Study2, cor_tot_Study1_vs_tot_Study2)
)

# Generate output filename using the last 2 sectors of the input file name
output_filename <- "outcomes_btw_mo-fa_corr.csv"
csv_file <- file.path(output_dir, output_filename)

# Save the correlation matrix as a CSV file with the new name
write.csv(cor_results, file = csv_file, row.names = FALSE)

cat("✅ Saved correlation matrix to", csv_file, "\n")

