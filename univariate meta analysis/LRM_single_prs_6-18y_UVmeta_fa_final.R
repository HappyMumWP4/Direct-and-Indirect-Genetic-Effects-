################################################################################
########################## Combine Multiple CSV Files ##########################
################################################################################

# Install necessary packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("tidyverse")

# # Set the directory where your CSV files are stored
# setwd("path_to_your_directory") # Set to correct directory path 

# Load necessary libraries
library(readr)
library(dplyr)

# List all CSV files in the directory
file_list <- list.files(pattern = ".*\\.csv$")

# Initialize an empty list to store the combined data
combined_data_list <- list()

# Loop through each file
for (file in file_list) {
  # Read the current CSV file
  data <- read_csv(file, skip_empty_rows = TRUE)  # Automatically skips empty rows
  
  # Ensure the necessary columns are present
  if (all(c("CBCL_outcome", "PRS_predictor", "Beta", "SE", "N") %in% colnames(data))) {
    # Extract the study name from the last segment of the file name
    study_name <- tools::file_path_sans_ext(file)  # Remove ".csv"
    study_name <- sub(".*_", "", study_name)  # Extract the last part after the last underscore
    
    # Add the study name as a new column
    data$StudyName <- study_name
    
    # Add the modified data frame to the list
    combined_data_list[[file]] <- data
  } else {
    cat("Skipping file (missing expected columns):", file, "\n")
  }
}

# Combine all the data frames in the list into one large data frame
combined_data <- bind_rows(combined_data_list)

# View a snapshot of the combined data
head(combined_data)

cat("Cohort-level analysis data are combined.")


################################################################################
########## Meta-Analysis For Fixed Effect_Inverse-Variance-Weighted ############
################################################################################

# Install necessary packages if not already installed
if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

# # Set the directory where your CSV files are stored
# setwd("path_to_your_directory") # Set to correct directory path 

# Load necessary libraries
library(metafor)
library(dplyr)
library(ggplot2)

# Define the directory where you want to save the forest plot PNG files
output_dir <- "lrm_single_prs_6-18y_fa_UVmeta_FE"  # Change directory name according to the statistical model

# Ensure the directory exists (optional, creates the directory if it doesn't exist)
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define a function to perform fixed-effect meta-analysis using inverse-variance weights and create forest plots
perform_meta_analysis_and_plot_FE <- function(outcome, predictors_pattern, output_dir) {
  # Filter the data based on the outcome
  outcome_data <- subset(combined_data, CBCL_outcome == outcome)
  
  # Get the unique predictors that match the specified pattern
  unique_predictors <- unique(grep(predictors_pattern, outcome_data$PRS_predictor, value = TRUE))
  
  # Extract relevant parts from outcome names
  outcome_parts <- unlist(strsplit(outcome, "_"))
  short_outcome <- outcome_parts[3]  # Third sector
  
  # Initialize lists to store results
  fe_inverse_variance_results <- list()
  p_values <- numeric(length(unique_predictors))  # Store p-values for each predictor
  
  # Initialize an empty data frame to store RE-IV results
  FE_IV_results_df <- data.frame(
    Outcome = character(),# Add column for the CBCL outcome
    Predictor = character(),
    Estimate = numeric(),
    SE = numeric(),
    Z_val = numeric(),
    p_val = numeric(),
    FDR_pval = numeric(),  # Add column for FDR-adjusted p-values
    CI_lb = numeric(),
    CI_ub = numeric(),
    I2 = numeric(),
    Q = numeric(),
    Q_p_val = numeric(),
    Num_Studies = numeric(),  
    Total_N = numeric(),      
    stringsAsFactors = FALSE
  )
  
  # Loop through each predictor and perform meta-analysis, including forest plot generation
  for (i in seq_along(unique_predictors)) {
    predictor <- unique_predictors[i]
    
    # Filter the data for the current predictor
    predictor_data <- subset(outcome_data, PRS_predictor == predictor)
    
    # Perform fixed-effect meta-analysis using inverse-variance weights
    res_fe_iv <- rma.uni(yi = predictor_data$Beta, sei = predictor_data$SE, method = "FE")  # Fixed-effect model
    
    # Store results and p-value
    fe_inverse_variance_results[[predictor]] <- summary(res_fe_iv)
    p_values[i] <- res_fe_iv$pval
    
    # Extract key metrics from the FE-IV results
    estimate_val <- res_fe_iv$beta
    se_val <- res_fe_iv$se
    z_val <- res_fe_iv$zval
    p_val <- res_fe_iv$pval
    ci_lb_val <- res_fe_iv$ci.lb
    ci_ub_val <- res_fe_iv$ci.ub
    I2_val <- res_fe_iv$I2
    Q_val <- res_fe_iv$QE
    Q_pval <- res_fe_iv$QEp
    
    # Compute additional columns
    num_studies <- length(unique(predictor_data$StudyName))  # Number of unique studies
    total_n <- sum(predictor_data$N)  # Total sample size
    
    # Add the data to the data frame
    FE_IV_results_df <- rbind(FE_IV_results_df, data.frame(
      Outcome = outcome,     # Add the CBCL outcome
      Predictor = predictor,
      Estimate = estimate_val,
      SE = se_val,
      Z_val = z_val,
      p_val = p_val,
      FDR_pval = NA,  # Placeholder for FDR-adjusted p-values
      CI_lb = ci_lb_val,
      CI_ub = ci_ub_val,
      I2 = I2_val,
      Q = Q_val,
      Q_p_val = Q_pval,
      Num_Studies = num_studies,  
      Total_N = total_n,          
      stringsAsFactors = FALSE
    ))
  }
  
  # Apply multiple testing correction (FDR)
  adjusted_pvals_fdr_fe <- p.adjust(p_values, method = "fdr", n = length(unique_predictors))
  
  # Add the adjusted FDR p-values to the data frame
  FE_IV_results_df$FDR_pval <- adjusted_pvals_fdr_fe
  
  # Loop again for forest plot generation based on FDR-adjusted p-values
  for (i in seq_along(unique_predictors)) {
    predictor <- unique_predictors[i]
    
    # Extract relevant parts from predictor names
    predictor_parts <- unlist(strsplit(predictor, "_"))
    short_predictor <- paste(predictor_parts[c(1,2,4)], collapse = "_")  # First 2 and 4th sector
    
    # Generate and save the forest plot if significant (FDR-corrected p < 0.05)
    if (adjusted_pvals_fdr_fe[i] < 0.05) {
      
      # Get study-level data for this predictor & outcome
      outcome_df <- outcome_data %>% 
        filter(PRS_predictor == predictor) %>%
        select(StudyName, PRS_predictor, CBCL_outcome, Beta, CI_lower, CI_upper, N)  # Select only relevant columns
      
      if (nrow(outcome_df) < 2) next
      
      # Get UVMA result
      uvma_results <- FE_IV_results_df %>% filter(Outcome == outcome, Predictor == predictor)
      
      if (nrow(uvma_results) == 0) {
        message("No UVMA results for outcome ", outcome, ". Skipping plot.")
        next
      }
      
      # Add UVMA result as a separate row
      uvma_row <- data.frame(
        StudyName = "UVMA Summary",
        PRS_predictor = uvma_results$Predictor,
        CBCL_outcome = uvma_results$Outcome,
        Beta = uvma_results$Estimate,
        CI_lower = uvma_results$CI_lb,
        CI_upper = uvma_results$CI_ub,
        N = uvma_results$Total_N
      )
      
      # Combine study-level data with UVMA summary
      plot_data <- bind_rows(outcome_df, uvma_row) 
      
      # Get all study names EXCEPT "UVMA Summary" and sort them alphabetically
      sorted_studies <- sort(setdiff(unique(plot_data$StudyName), "UVMA Summary"))
      
      # Reassign factor levels with "UVMA Summary" first, then sorted studies
      plot_data$StudyName <- factor(
        plot_data$StudyName,
        levels = c("UVMA Summary", sorted_studies)  # Ensures summary is on top, others in A-Z order
      )
      
      # Define output file
      file_name <- paste0("UVFE_lrm_", short_predictor, "_", short_outcome, "_fa", ".png") # Change file name according to the statistical model
      file_path <- file.path(output_dir, file_name)
      
      # Prepare labels for effect sizes and confidence intervals
      plot_data$label <- sprintf("%.2f (%.2f, %.2f)", plot_data$Beta, plot_data$CI_lower, plot_data$CI_upper)
      
      # Save plot
      g <- ggplot(plot_data, aes(y = StudyName, x = Beta, xmin = CI_lower, xmax = CI_upper)) +
        geom_point(aes(size = N, color = StudyName == "UVMA Summary"), alpha = 0.8) +  # Point size reflects sample size
        geom_errorbarh(aes(color = StudyName == "UVMA Summary"), height = 0.4) +  # Horizontal error bars
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +  # Zero-line
        geom_text(aes(label = label), hjust = -0.2, vjust = 0, size = 4) +  # Adjusted label position
        scale_color_manual(values = c("black", "red"), guide = "none") +  # Red for UVMA
        scale_size_continuous(
          range = c(3, 8),  # Adjust point size range
          breaks = c(100, 200, 500, 1000, 2000, 4000),  # Define custom breaks
          labels = c("100", "200", "500", "1K", "2K", "4K"),  # Format labels
          limits = c(100, 4000)  # Ensure the full range is considered
        ) + 
        
        labs(x = "Beta Coefficient", y = "Study", 
             title = paste("Forest Plot (FE) of Single-LRM for", outcome, "and", predictor),
             size = "Sample Size") +
        
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 14, hjust = 0.5, vjust = 1, face = "bold"),  # Ensure title is centered and bold
          plot.margin = margin(10, 10, 10, 10)  # Increase margins to prevent cropping
        )
      
      # Save the plot
      ggsave(file_path, plot = g + theme(plot.title.position = "plot"), width = 8, height = 6, dpi = 300)
      
      cat("Forest plot saved for Variable:", predictor, "as", file_path, "\n")
    } else {
      cat("No significant effect after FDR correction for Variable:", predictor, ". Skipping plot.\n")
    }
  }
  
  # Save all the FE-IV results (including FDR-adjusted p-values) as a CSV file
  print(head(FE_IV_results_df))
  
  csv_file_path <- file.path(output_dir, paste0("lrm_", short_outcome, "_", predictors_pattern, "UVmeta_FE_fa.csv")) # Change file name according to the statistical model
  write.csv(FE_IV_results_df, csv_file_path, row.names = FALSE)
  
  cat("FE-IV results (including FDR-adjusted p-values) saved to:", csv_file_path, "\n")
  
  return(FE_IV_results_df)
}

# Initialize a list to store results for all (outcome, predictor pattern) combinations
all_results <- list()

# Get unique CBCL_outcome values
unique_outcomes <- unique(combined_data$CBCL_outcome)

# Define predictor patterns
patterns <- c("prs_c_qc_", "prs_m_qc_", "prs_f_qc_", "prs_nt_qc_", "pTDTd_qc_")

# Loop through each outcome and predictor pattern
for (outcome in unique_outcomes) {
  for (pattern in patterns) {
    print(paste("Processing:", outcome, "with predictor pattern:", pattern))
    
    # Perform meta-analysis and forest plotting
    result_df <- perform_meta_analysis_and_plot_FE(outcome, pattern, output_dir)  
    
    # Store results in a list using a unique key
    key <- paste(outcome, pattern, sep = "_")
    all_results[[key]] <- result_df  
  }
}

# Combine all results into a single data frame
FE_final_results_df <- bind_rows(all_results, .id = "Outcome_PredictorPattern")

cat("FE-IV results (including FDR-adjusted p-values) saved")


################################################################################
########## Meta-Analysis For Random Effect_Inverse-Variance-Weighted ###########
################################################################################

# Install necessary packages if not already installed
if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

# # Set the directory where your CSV files are stored
# setwd("path_to_your_directory") # Set to correct directory path 

# Load necessary libraries
library(metafor)
library(dplyr)
library(ggplot2)

# Define the directory where you want to save the forest plot PNG files
output_dir <- "lrm_single_prs_6-18y_fa_UVmeta_RE"  # Change directory name according to the statistical model

# Ensure the directory exists (optional, creates the directory if it doesn't exist)
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define a function to perform random-effect meta-analysis using inverse-variance weights and create forest plots
perform_meta_analysis_and_plot_RE <- function(outcome, predictors_pattern, output_dir) {
  # Filter the data based on the outcome
  outcome_data <- subset(combined_data, CBCL_outcome == outcome)
  
  # Get the unique predictors that match the specified pattern
  unique_predictors <- unique(grep(predictors_pattern, outcome_data$PRS_predictor, value = TRUE))
  
  # Extract relevant parts from outcome names
  outcome_parts <- unlist(strsplit(outcome, "_"))
  short_outcome <- outcome_parts[3]  # Third sector
  
  # Initialize lists to store results
  re_inverse_variance_results <- list()
  p_values <- numeric(length(unique_predictors))  # Store p-values for each predictor
  
  # Initialize an empty data frame to store RE-IV results
  RE_IV_results_df <- data.frame(
    Outcome = character(),# Add column for the CBCL outcome
    Predictor = character(),
    Estimate = numeric(),
    SE = numeric(),
    Z_val = numeric(),
    p_val = numeric(),
    FDR_pval = numeric(),  # Add column for FDR-adjusted p-values
    CI_lb = numeric(),
    CI_ub = numeric(),
    I2 = numeric(),
    Q = numeric(),
    Q_p_val = numeric(),
    Num_Studies = numeric(),  
    Total_N = numeric(),      
    stringsAsFactors = FALSE
  )
  
  # Loop through each predictor and perform meta-analysis, including forest plot generation
  for (i in seq_along(unique_predictors)) {
    predictor <- unique_predictors[i]
    
    # Filter the data for the current predictor
    predictor_data <- subset(outcome_data, PRS_predictor == predictor)
    
    # Perform random-effect meta-analysis using inverse-variance weights
    res_re_iv <- rma.uni(yi = predictor_data$Beta, sei = predictor_data$SE, method = "REML")  # Random-effect model
    
    # Store results and p-value
    re_inverse_variance_results[[predictor]] <- summary(res_re_iv)
    p_values[i] <- res_re_iv$pval
    
    # Extract key metrics from the RE-IV results
    estimate_val <- res_re_iv$beta
    se_val <- res_re_iv$se
    z_val <- res_re_iv$zval
    p_val <- res_re_iv$pval
    ci_lb_val <- res_re_iv$ci.lb
    ci_ub_val <- res_re_iv$ci.ub
    I2_val <- res_re_iv$I2
    Q_val <- res_re_iv$QE
    Q_pval <- res_re_iv$QEp
    
    # Compute additional columns
    num_studies <- length(unique(predictor_data$StudyName))  # Number of unique studies
    total_n <- sum(predictor_data$N)  # Total sample size
    
    # Add the data to the data frame
    RE_IV_results_df <- rbind(RE_IV_results_df, data.frame(
      Outcome = outcome,     # Add the CBCL outcome
      Predictor = predictor,
      Estimate = estimate_val,
      SE = se_val,
      Z_val = z_val,
      p_val = p_val,
      FDR_pval = NA,  # Placeholder for FDR-adjusted p-values
      CI_lb = ci_lb_val,
      CI_ub = ci_ub_val,
      I2 = I2_val,
      Q = Q_val,
      Q_p_val = Q_pval,
      Num_Studies = num_studies,  
      Total_N = total_n,          
      stringsAsFactors = FALSE
    ))
  }
  
  # Apply multiple testing correction (FDR)
  adjusted_pvals_fdr_re <- p.adjust(p_values, method = "fdr", n = length(unique_predictors))
  
  # Add the adjusted FDR p-values to the data frame
  RE_IV_results_df$FDR_pval <- adjusted_pvals_fdr_re
  
  # Loop again for forest plot generation based on FDR-adjusted p-values
  for (i in seq_along(unique_predictors)) {
    predictor <- unique_predictors[i]
    
    # Extract relevant parts from predictor names
    predictor_parts <- unlist(strsplit(predictor, "_"))
    short_predictor <- paste(predictor_parts[c(1,2,4)], collapse = "_")  # First 2 and 4th sector
    
    # Generate and save the forest plot if significant (FDR-corrected p < 0.05)
    if (adjusted_pvals_fdr_re[i] < 0.05) {
      
      # Get study-level data for this predictor & outcome
      outcome_df <- outcome_data %>% 
        filter(PRS_predictor == predictor) %>%
        select(StudyName, PRS_predictor, CBCL_outcome, Beta, CI_lower, CI_upper, N)  # Select only relevant columns
      
      if (nrow(outcome_df) < 2) next
      
      # Get UVMA result
      uvma_results <- RE_IV_results_df %>% filter(Outcome == outcome, Predictor == predictor)
      
      if (nrow(uvma_results) == 0) {
        message("No UVMA results for outcome ", outcome, ". Skipping plot.")
        next
      }
      
      # Add UVMA result as a separate row
      uvma_row <- data.frame(
        StudyName = "UVMA Summary",
        PRS_predictor = uvma_results$Predictor,
        CBCL_outcome = uvma_results$Outcome,
        Beta = uvma_results$Estimate,
        CI_lower = uvma_results$CI_lb,
        CI_upper = uvma_results$CI_ub,
        N = uvma_results$Total_N
      )
      
      # Combine study-level data with UVMA summary
      plot_data <- bind_rows(outcome_df, uvma_row) 
      
      # Get all study names EXCEPT "UVMA Summary" and sort them alphabetically
      sorted_studies <- sort(setdiff(unique(plot_data$StudyName), "UVMA Summary"))
      
      # Reassign factor levels with "UVMA Summary" first, then sorted studies
      plot_data$StudyName <- factor(
        plot_data$StudyName,
        levels = c("UVMA Summary", sorted_studies)  # Ensures summary is on top, others in A-Z order
      )
      
      # Define output file
      file_name <- paste0("UVRE_lrm_", short_predictor, "_", short_outcome, "_fa", ".png") # Change file name according to the statistical model
      file_path <- file.path(output_dir, file_name)
      
      # Prepare labels for effect sizes and confidence intervals
      plot_data$label <- sprintf("%.2f (%.2f, %.2f)", plot_data$Beta, plot_data$CI_lower, plot_data$CI_upper)
      
      # Save plot
      g <- ggplot(plot_data, aes(y = StudyName, x = Beta, xmin = CI_lower, xmax = CI_upper)) +
        geom_point(aes(size = N, color = StudyName == "UVMA Summary"), alpha = 0.8) +  # Point size reflects sample size
        geom_errorbarh(aes(color = StudyName == "UVMA Summary"), height = 0.4) +  # Horizontal error bars
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +  # Zero-line
        geom_text(aes(label = label), hjust = -0.2, vjust = 0, size = 4) +  # Adjusted label position
        scale_color_manual(values = c("black", "red"), guide = "none") +  # Red for UVMA
        scale_size_continuous(
          range = c(3, 8),  # Adjust point size range
          breaks = c(100, 200, 500, 1000, 2000, 4000),  # Define custom breaks
          labels = c("100", "200", "500", "1K", "2K", "4K"),  # Format labels
          limits = c(100, 4000)  # Ensure the full range is considered
        ) + 
        
        labs(x = "Beta Coefficient", y = "Study", 
             title = paste("Forest Plot (RE) of Single-LRM for", outcome, "and", predictor),
             size = "Sample Size") +
        
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 14, hjust = 0.5, vjust = 1, face = "bold"),  # Ensure title is centered and bold
          plot.margin = margin(10, 10, 10, 10)  # Increase margins to prevent cropping
        )
      
      # Save the plot
      ggsave(file_path, plot = g + theme(plot.title.position = "plot"), width = 8, height = 6, dpi = 300)
      
      cat("Forest plot saved for Variable:", predictor, "as", file_path, "\n")
    } else {
      cat("No significant effect after FDR correction for Variable:", predictor, ". Skipping plot.\n")
    }
  }
  
  # Save all the RE-IV results (including FDR-adjusted p-values) as a CSV file
  print(head(RE_IV_results_df))
  
  csv_file_path <- file.path(output_dir, paste0("lrm_", short_outcome, "_", predictors_pattern, "UVmeta_RE_fa.csv")) # Change file name according to the statistical model
  write.csv(RE_IV_results_df, csv_file_path, row.names = FALSE)
  
  cat("RE-IV results (including FDR-adjusted p-values) saved to:", csv_file_path, "\n")
  
  return(RE_IV_results_df)
}

# Initialize a list to store results for all (outcome, predictor pattern) combinations
all_results <- list()

# Get unique CBCL_outcome values
unique_outcomes <- unique(combined_data$CBCL_outcome)

# Define predictor patterns
patterns <- c("prs_c_qc_", "prs_m_qc_", "prs_f_qc_", "prs_nt_qc_", "pTDTd_qc_")

# Loop through each outcome and predictor pattern
for (outcome in unique_outcomes) {
  for (pattern in patterns) {
    print(paste("Processing:", outcome, "with predictor pattern:", pattern))
    
    # Perform meta-analysis and forest plotting
    result_df <- perform_meta_analysis_and_plot_RE(outcome, pattern, output_dir)  
    
    # Store results in a list using a unique key
    key <- paste(outcome, pattern, sep = "_")
    all_results[[key]] <- result_df  
  }
}

# Combine all results into a single data frame
RE_final_results_df <- bind_rows(all_results, .id = "Outcome_PredictorPattern")

cat("RE-IV results (including FDR-adjusted p-values) saved")


################################################################################
################## Save meta-analysis objects for future use ##################
################################################################################

# Save the combined input data, meta-analysis models and output statistics for future use
meta_results <- list(combined_data, FE_final_results_df, RE_final_results_df)
save(meta_results, file = "LRM_single_prs_6-18y_fa_UVmeta.RData")

cat("Combined data and output results saved to 'LRM_single_prs_6-18y_fa_UVmeta.RData'.\n") # Change file name according to the statistical model


################################################################################
################### Save combined input data for future use ####################
################################################################################

# Save the combined data to a new CSV file
write_csv(combined_data, "lrm_single_prs_6-18y_fa_combined_studies.csv") # Change file name according to the statistical model

cat("Combined data saved to 'lrm_single_prs_6-18y_fa_combined_studies.csv'.\n") # Change file name according to the statistical model

