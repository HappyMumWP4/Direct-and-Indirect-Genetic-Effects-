################################################################################
########################## Combine Multiple CSV Files ##########################
################################################################################

# Install necessary packages if not installed
packages <- c("Matrix", "metafor", "tidyverse")
installed_packages <- rownames(installed.packages())

for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# # Set the directory where your CSV files are stored
# setwd("path_to_your_directory") # Set to correct directory path 

# Load necessary libraries
library(dplyr)

# List all CSV files in the directory
csv_files <- list.files(pattern = "*.csv", full.names = TRUE)

# Read and combine all CSV files
meta_data_list <- lapply(csv_files, function(file) {
  # Read CSV file
  data <- read.csv(file)
  
  # Ensure the necessary columns are present
  if (all(c("CBCL_outcome", "PRS_predictor", "Beta", "SE", "N") %in% colnames(data))) {
    # Extract StudyName from the file name
    file_name <- basename(file)  # Get file name without directory
    file_name <- sub("\\.csv$", "", file_name)  # Remove ".csv"
    parts <- unlist(strsplit(file_name, "_"))  # Split by underscores
    study_name <- paste(parts[(length(parts)-1):length(parts)], collapse = "_")  # Last 2 parts
    
    # Add StudyName column
    data$StudyName <- study_name
    
  } else {
    cat("Skipping file (missing expected columns):", file, "\n")
  }
  
  return(data)
})

# Combine all datasets into one
meta_data <- do.call(rbind, meta_data_list)

# Check the first few rows
head(meta_data)

cat("Cohort-level analysis data are combined.")


################################################################################
################### MVmeta_FE m+f: Ext., Int., Tot. Prob. ####################
################################################################################

# Load necessary libraries
library(metafor)
library(dplyr)
library(Matrix)
library(ggplot2)

# Define output directory
output_dir <- "MVMA_FE_1step_m+f_lrm_trio_all_6-18y"

# Ensure the directory is created, with recursive = TRUE in case of nested folders
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define between-study correlations per outcome
rho_values <- list(
  "h_c_total_6_18y" = c("GenR3" = 0.402, "GenR4" = 0.38, "NTR" = 0.698, "PREDO" = 0.432),
  "h_c_extern_6_18y" = c("GenR3" = 0.558, "GenR4" = 0.53, "NTR" = 0.675, "PREDO" = 0.452),
  "h_c_intern_6_18y" = c("GenR3" = 0.512, "GenR4" = 0.477, "NTR" = 0.616, "PREDO" = 0.376)
)

# Define correlated study pairs
correlated_pairs <- list(
  "GenR3" = c("fa_GenR3", "mo_GenR3"),
  "GenR4" = c("fa_GenR4", "mo_GenR4"),
  "NTR" = c("fa_NTR", "mo_NTR"),
  "PREDO" = c("fa_PREDO", "mo_PREDO")
)

# Define study lists
full_studies <- c("mo_ALSPAC", "fa_GenR3", "mo_GenR3", "fa_GenR4", "mo_GenR4", "mix_MCS", 
                  "mo_MoBa", "fa_NTR", "mo_NTR", "mo_POSEIDON", "fa_PREDO", "mo_PREDO")
anx_studies <- setdiff(full_studies, c("fa_NTR", "mo_NTR")) # Exclude NTR studies for ANX

# Get PRS predictors
prs_predictors <- unique(meta_data$PRS_predictor)

# Prepare output lists
output_list <- list()
meta_results <- list()
fdr_list <- list()

# Loop over each outcome
for (outcome in names(rho_values)) {
  rho_between <- rho_values[[outcome]]
  outcome_output_list <- list()
  
  for (prs in prs_predictors) {
    # Subset data for given PRS and outcome
    meta_subset <- subset(meta_data, PRS_predictor == prs & CBCL_outcome == outcome)
    
    # Determine which study list to use
    if (grepl("_qc_ANX_trio", prs)) {
      study_list <- anx_studies
      required_studies <- 10 # Define selected number studies for certain prs!
    } else {
      study_list <- full_studies
      required_studies <- 12 # Define full number of studies! 
    }
    
    # Ensure order matches expected study layout
    meta_subset <- meta_subset[match(study_list, meta_subset$StudyName), ]
    
    # Check if required number of studies are present
    if (nrow(meta_subset) < required_studies || any(is.na(meta_subset$SE))) next
    
    # Initialize variance-covariance matrix
    V_matrix <- diag(meta_subset$SE^2)
    
    # Insert covariances for correlated study pairs (only for available studies)
    for (group in names(correlated_pairs)) {
      pair <- correlated_pairs[[group]]
      # Skip if either study in pair is missing
      if (!all(pair %in% meta_subset$StudyName)) next
      
      rho <- rho_between[group]
      if (!is.na(rho)) {
        idx1 <- which(meta_subset$StudyName == pair[1])
        idx2 <- which(meta_subset$StudyName == pair[2])
        cov_val <- rho * meta_subset$SE[idx1] * meta_subset$SE[idx2]
        V_matrix[idx1, idx2] <- cov_val
        V_matrix[idx2, idx1] <- cov_val
      }
    }
    
    # Print V_matrix for the first PRS-outcome combo
    if (prs == prs_predictors[1]) {
      cat("\n--- Variance-Covariance Matrix for:", prs, "|", outcome, "---\n")
      print(V_matrix)
    }
    
    # Fit FE MVMA
    mv_model <- rma.mv(
      yi = meta_subset$Beta,
      V = V_matrix,
      random = NULL,
      data = meta_subset
    )
    
    # Calculate I²
    I2 <- 0
    
    # Store result with dynamic study info
    res_df <- data.frame(
      CBCL_outcome = outcome,
      PRS_predictor = prs,
      Beta = mv_model$b,
      SE = mv_model$se,
      CI_Lower = mv_model$ci.lb,
      CI_Upper = mv_model$ci.ub,
      P_Value = mv_model$pval,
      I2 = I2,
      Num_Studies = nrow(meta_subset), # Use actual number of studies
      N = sum(meta_subset$N),
      Analysis_Type = ifelse(grepl("_qc_ANX_trio", prs), "-NTR", "Full") # Note on difference in study numbers for certain prs
    )
    
    outcome_output_list[[paste(prs, outcome, sep = "_")]] <- res_df
    
    # Store model & summary
    meta_results[[paste(prs, outcome, sep = "_")]] <- list(
      model = mv_model,
      summary = summary(mv_model),
      variance_components = list(
        tau2 = sum(mv_model$tau2),  # Sum of all tau² components
        sigma2 = mean(diag(mv_model$V)),  # Mean sampling variance
        I2 = I2
      ),
      correlation_param = rho_between
    )
    
    # Organize for FDR
    prs_set <- substr(prs, 1, 5)
    fdr_key <- paste(prs_set, outcome, sep = "_")
    fdr_list[[fdr_key]] <- rbind(fdr_list[[fdr_key]], res_df)
  }
  
  # Save results for each outcome
  output_list[[outcome]] <- outcome_output_list
}

# Apply FDR correction within each PRS group per outcome and create forest plot for significant result
fdr_corrected_results <- list()

for (key in names(fdr_list)) {
  fdr_df <- fdr_list[[key]]
  fdr_df$FDR_p <- p.adjust(fdr_df$P_Value, method = "fdr")
  fdr_corrected_results[[key]] <- fdr_df
}

# Combine FDR-adjusted results back to a single data frame
final_results <- do.call(rbind, fdr_corrected_results)
  
# Generate and save the forest plot if significant (FDR-corrected p < 0.05)
for (key in names(fdr_corrected_results)) {
  fdr_df <- fdr_corrected_results[[key]]
  
  # Filter rows with FDR p < 0.05
  sig_df <- fdr_df %>% filter(FDR_p < 0.05)
  
  if (nrow(sig_df) == 0) next
  
  for (i in seq_len(nrow(sig_df))) {
    row <- sig_df[i, ]
    prs <- row$PRS_predictor
    outcome <- row$CBCL_outcome
    
    # Short names for plot filename
    predictor_parts <- unlist(strsplit(prs, "_"))
    short_predictor <- paste(predictor_parts[c(1,2,4)], collapse = "_")
    
    outcome_parts <- unlist(strsplit(outcome, "_"))
    short_outcome <- outcome_parts[3]
    
    # Get study-level data
    outcome_df <- meta_data %>% 
      filter(PRS_predictor == prs, CBCL_outcome == outcome) %>%
      select(StudyName, PRS_predictor, CBCL_outcome, Beta, CI_lower, CI_upper, N)
    
    if (nrow(outcome_df) < 2) next
    
    # Extract MVMA summary result
    res_key <- paste(prs, outcome, sep = "_")
    mvma_res <- output_list[[outcome]][[res_key]][1, ]
    
    mvma_row <- data.frame(
      StudyName = "MVMA Summary",
      PRS_predictor = mvma_res$PRS_predictor,
      CBCL_outcome = mvma_res$CBCL_outcome,
      Beta = mvma_res$Beta,
      CI_lower = mvma_res$CI_Lower,
      CI_upper = mvma_res$CI_Upper,
      N = mvma_res$N
    )
    
    # Combine study data and MVMA
    plot_data <- bind_rows(mvma_row, outcome_df)
    
    # Extract cohort names (e.g., "GenR3" from "fa_GenR3")
    plot_data <- plot_data %>%
      mutate(
        cohort = gsub("^[^_]*_", "", StudyName),  # Removes "fa_" / "mo_" prefix
        reporter = ifelse(grepl("^fa_", StudyName), "fa", "mo")
      )
    
    # Define custom sorting order for cohorts (adjust as needed)
    cohort_order <- c(
      "ALSPAC", "GenR3", "GenR4", "MCS", "MoBa", 
      "NTR", "POSEIDON", "PREDO"
    )                                        # Change study names if necessary!
    
    # Sort studies: cohort order -> reporter (fa before mo)
    sorted_studies <- plot_data %>%
      filter(StudyName != "MVMA Summary") %>%
      mutate(
        cohort = factor(cohort, levels = cohort_order),
        reporter = factor(reporter, levels = c("fa", "mo"))
      ) %>%
      arrange(cohort, reporter) %>%
      pull(StudyName) %>%
      unique()
    
    # Apply ordering (MVMA Summary first, then grouped studies)
    plot_data$StudyName <- factor(
      plot_data$StudyName,
      levels = c("MVMA Summary", sorted_studies)
    )
    
    # Labels
    plot_data$label <- sprintf("%.2f (%.2f, %.2f)", plot_data$Beta, plot_data$CI_lower, plot_data$CI_upper)
    
    # Output file
    file_name <- paste0("MVMA_FE_lrm_trio_", short_predictor, "_", short_outcome, ".png")
    file_path <- file.path(output_dir, file_name)
    
    # Create plot
    g <- ggplot(plot_data, aes(y = StudyName, x = Beta, xmin = CI_lower, xmax = CI_upper)) +
      geom_point(aes(size = N, color = StudyName == "MVMA Summary"), alpha = 0.8) +
      geom_errorbarh(aes(color = StudyName == "MVMA Summary"), height = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_text(aes(label = label), hjust = -0.2, vjust = 0, size = 4) +
      scale_color_manual(values = c("black", "red"), guide = "none") +
      scale_size_continuous(
        range = c(3, 8),
        breaks = c(100, 200, 500, 1000, 2000, 5000, 11000),
        labels = c("100", "200", "500", "1K", "2K", "5K", "11K"),
        limits = c(100, 12000)
      ) +
      labs(x = "Beta Coefficient", y = "Study",
           title = paste("Forest Plot (FE) of Trio-LRM for", outcome, "and", prs),
           size = "Sample Size") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
            plot.margin = margin(10, 10, 10, 10))
    
    # Save plot
    ggsave(file_path, plot = g + theme(plot.title.position = "plot"), width = 8, height = 6, dpi = 300)
    
    cat("Forest plot saved for:", prs, "and", outcome, "\n")
  }
}

# Save results to CSV
output_file <- file.path(output_dir, "MVMA_FE_1step_lrm_trio_all_6-18y.csv")
write.csv(final_results, file = output_file, row.names = FALSE)

cat("Final results saved to", output_file, "\n")

# Save the meta-analysis models & statistics for future use
save(meta_results, file = file.path(output_dir, "MVMA_FE_1step_lrm_trio_all_6-18y.RData")) # Change file name according to the model and cohort
 
cat("Meta-analysis objects saved as 'MVMA_FE_1step_lrm_trio_all_6-18y.RData'.\n")


################################################################################
################### MVmeta_RE m+f: Ext., Int., Tot. Prob. ####################
################################################################################

# # Install orchaRd from Source (Advanced)
# if (!require("remotes")) install.packages("remotes")
# remotes::install_github("daniel1noble/orchaRd")

# Load necessary libraries
library(metafor)
library(dplyr)
library(Matrix)
library(orchaRd)
library(ggplot2)

# Define output directory
output_dir <- "MVMA_RE_1step_m+f_lrm_trio_all_6-18y"

# Ensure the directory is created, with recursive = TRUE in case of nested folders
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define between-study correlations per outcome
rho_values <- list(
  "h_c_total_6_18y" = c("GenR3" = 0.402, "GenR4" = 0.38, "NTR" = 0.698, "PREDO" = 0.432),
  "h_c_extern_6_18y" = c("GenR3" = 0.558, "GenR4" = 0.53, "NTR" = 0.675, "PREDO" = 0.452),
  "h_c_intern_6_18y" = c("GenR3" = 0.512, "GenR4" = 0.477, "NTR" = 0.616, "PREDO" = 0.376)
)

# Define correlated study pairs
correlated_pairs <- list(
  "GenR3" = c("fa_GenR3", "mo_GenR3"),
  "GenR4" = c("fa_GenR4", "mo_GenR4"),
  "NTR" = c("fa_NTR", "mo_NTR"),
  "PREDO" = c("fa_PREDO", "mo_PREDO")
)

# Define study lists
full_studies <- c("mo_ALSPAC", "fa_GenR3", "mo_GenR3", "fa_GenR4", "mo_GenR4", "mix_MCS", 
                  "mo_MoBa", "fa_NTR", "mo_NTR", "mo_POSEIDON", "fa_PREDO", "mo_PREDO")
anx_studies <- setdiff(full_studies, c("fa_NTR", "mo_NTR")) # Exclude NTR studies for ANX

# Get PRS predictors
prs_predictors <- unique(meta_data$PRS_predictor)

# Prepare output lists
output_list <- list()
meta_results <- list()
fdr_list <- list()

# Loop over each outcome
for (outcome in names(rho_values)) {
  rho_between <- rho_values[[outcome]]
  outcome_output_list <- list()
  
  for (prs in prs_predictors) {
    # Subset data for given PRS and outcome
    meta_subset <- subset(meta_data, PRS_predictor == prs & CBCL_outcome == outcome)
    
    # Determine which study list to use
    if (grepl("_qc_ANX_trio", prs)) {
      study_list <- anx_studies
      required_studies <- 10 # Define selected number studies for certain prs!
    } else {
      study_list <- full_studies
      required_studies <- 12 # Define full number of studies! 
    }
    
    # Ensure order matches expected study layout
    meta_subset <- meta_subset[match(study_list, meta_subset$StudyName), ]
    
    # Check if required number of studies are present
    if (nrow(meta_subset) < required_studies || any(is.na(meta_subset$SE))) next
    
    # Initialize variance-covariance matrix
    V_matrix <- diag(meta_subset$SE^2)
    
    # Insert covariances for correlated study pairs (only for available studies)
    for (group in names(correlated_pairs)) {
      pair <- correlated_pairs[[group]]
      # Skip if either study in pair is missing
      if (!all(pair %in% meta_subset$StudyName)) next
      
      rho <- rho_between[group]
      if (!is.na(rho)) {
        idx1 <- which(meta_subset$StudyName == pair[1])
        idx2 <- which(meta_subset$StudyName == pair[2])
        cov_val <- rho * meta_subset$SE[idx1] * meta_subset$SE[idx2]
        V_matrix[idx1, idx2] <- cov_val
        V_matrix[idx2, idx1] <- cov_val
      }
    }
    
    # Print V_matrix for the first PRS-outcome combo
    if (prs == prs_predictors[1]) {
      cat("\n--- Variance-Covariance Matrix for:", prs, "|", outcome, "---\n")
      print(V_matrix)
    }
    
    # Fit RE MVMA
    mv_model <- tryCatch({
      rma.mv(
        yi = meta_subset$Beta,
        V = V_matrix,
        random = ~ 1 | StudyName,
        method = "REML",
        control = list(iter.max = 1000, rel.tol = 1e-8),
        data = meta_subset
      )
    }, error = function(e) {
      cat("Error fitting model for", prs, "|", outcome, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(mv_model)) next
    
    # Calculate I² (I2 = (sum(tau2) / (sum(tau2) + mean(diag(V))) * 100)
    # total_var <- sum(mv_model$tau2) + mean(diag(mv_model$V))
    # I2 <- (sum(mv_model$tau2) / total_var) * 100
    # 
    # Multivariate I² calculated from the package ‘orchaRd’
    i2 <- i2_ml(mv_model)
    I2_orchaRd <- i2["I2_Total"]
    
    # Store result with dynamic study info
    res_df <- data.frame(
      CBCL_outcome = outcome,
      PRS_predictor = prs,
      Beta = mv_model$b,
      SE = mv_model$se,
      CI_Lower = mv_model$ci.lb,
      CI_Upper = mv_model$ci.ub,
      P_Value = mv_model$pval,
      I2_orchaRd = I2_orchaRd,
      Num_Studies = nrow(meta_subset), # Use actual number of studies
      N = sum(meta_subset$N),
      Analysis_Type = ifelse(grepl("_qc_ANX_trio", prs), "-NTR", "Full") # Note on difference in study numbers for certain prs
    )
    
    outcome_output_list[[paste(prs, outcome, sep = "_")]] <- res_df
    
    # Store model & summary
    meta_results[[paste(prs, outcome, sep = "_")]] <- list(
      model = mv_model,
      summary = summary(mv_model),
      variance_components = list(
        tau2 = sum(mv_model$tau2),  # Sum of all tau² components
        sigma2 = mean(diag(mv_model$V)),  # Mean sampling variance
        I2_orchaRd = I2_orchaRd
      ),
      correlation_param = rho_between
    )
    
    # Organize for FDR
    prs_set <- substr(prs, 1, 5)
    fdr_key <- paste(prs_set, outcome, sep = "_")
    fdr_list[[fdr_key]] <- rbind(fdr_list[[fdr_key]], res_df)
  }
  
  # Save results for each outcome
  output_list[[outcome]] <- outcome_output_list
}

# Apply FDR correction within each PRS group per outcome and create forest plot for significant result
fdr_corrected_results <- list()

for (key in names(fdr_list)) {
  fdr_df <- fdr_list[[key]]
  fdr_df$FDR_p <- p.adjust(fdr_df$P_Value, method = "fdr")
  fdr_corrected_results[[key]] <- fdr_df
}

# Combine FDR-adjusted results back to a single data frame
final_results <- do.call(rbind, fdr_corrected_results)

# Generate and save the forest plot if significant (FDR-corrected p < 0.05)
for (key in names(fdr_corrected_results)) {
  fdr_df <- fdr_corrected_results[[key]]
  
  # Filter rows with FDR p < 0.05
  sig_df <- fdr_df %>% filter(FDR_p < 0.05)
  
  if (nrow(sig_df) == 0) next
  
  for (i in seq_len(nrow(sig_df))) {
    row <- sig_df[i, ]
    prs <- row$PRS_predictor
    outcome <- row$CBCL_outcome
    
    # Short names for plot filename
    predictor_parts <- unlist(strsplit(prs, "_"))
    short_predictor <- paste(predictor_parts[c(1,2,4)], collapse = "_")
    
    outcome_parts <- unlist(strsplit(outcome, "_"))
    short_outcome <- outcome_parts[3]
    
    # Get study-level data
    outcome_df <- meta_data %>% 
      filter(PRS_predictor == prs, CBCL_outcome == outcome) %>%
      select(StudyName, PRS_predictor, CBCL_outcome, Beta, CI_lower, CI_upper, N)
    
    if (nrow(outcome_df) < 2) next
    
    # Extract MVMA summary result
    res_key <- paste(prs, outcome, sep = "_")
    mvma_res <- output_list[[outcome]][[res_key]][1, ]
    
    mvma_row <- data.frame(
      StudyName = "MVMA Summary",
      PRS_predictor = mvma_res$PRS_predictor,
      CBCL_outcome = mvma_res$CBCL_outcome,
      Beta = mvma_res$Beta,
      CI_lower = mvma_res$CI_Lower,
      CI_upper = mvma_res$CI_Upper,
      N = mvma_res$N
    )
    
    # Combine study data and MVMA
    plot_data <- bind_rows(mvma_row, outcome_df)
    
    # Extract cohort names (e.g., "GenR3" from "fa_GenR3")
    plot_data <- plot_data %>%
      mutate(
        cohort = gsub("^[^_]*_", "", StudyName),  # Removes "fa_" / "mo_" prefix
        reporter = ifelse(grepl("^fa_", StudyName), "fa", "mo")
      )
    
    # Define custom sorting order for cohorts (adjust as needed)
    cohort_order <- c(
      "ALSPAC", "GenR3", "GenR4", "MCS", "MoBa", 
      "NTR", "POSEIDON", "PREDO"
    )                                        # Change study names if necessary!
    
    # Sort studies: cohort order -> reporter (fa before mo)
    sorted_studies <- plot_data %>%
      filter(StudyName != "MVMA Summary") %>%
      mutate(
        cohort = factor(cohort, levels = cohort_order),
        reporter = factor(reporter, levels = c("fa", "mo"))
      ) %>%
      arrange(cohort, reporter) %>%
      pull(StudyName) %>%
      unique()
    
    # Apply ordering (MVMA Summary first, then grouped studies)
    plot_data$StudyName <- factor(
      plot_data$StudyName,
      levels = c("MVMA Summary", sorted_studies)
    )
    
    # Labels
    plot_data$label <- sprintf("%.2f (%.2f, %.2f)", plot_data$Beta, plot_data$CI_lower, plot_data$CI_upper)
    
    # Output file
    file_name <- paste0("MVMA_RE_lrm_trio_", short_predictor, "_", short_outcome, ".png")
    file_path <- file.path(output_dir, file_name)
    
    # Create plot
    g <- ggplot(plot_data, aes(y = StudyName, x = Beta, xmin = CI_lower, xmax = CI_upper)) +
      geom_point(aes(size = N, color = StudyName == "MVMA Summary"), alpha = 0.8) +
      geom_errorbarh(aes(color = StudyName == "MVMA Summary"), height = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_text(aes(label = label), hjust = -0.2, vjust = 0, size = 4) +
      scale_color_manual(values = c("black", "red"), guide = "none") +
      scale_size_continuous(
        range = c(3, 8),
        breaks = c(100, 200, 500, 1000, 2000, 5000, 11000),
        labels = c("100", "200", "500", "1K", "2K", "5K", "11K"),
        limits = c(100, 12000)
      ) +
      labs(x = "Beta Coefficient", y = "Study",
           title = paste("Forest Plot (RE) of Trio-LRM for", outcome, "and", prs),
           size = "Sample Size") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
            plot.margin = margin(10, 10, 10, 10))
    
    # Save plot
    ggsave(file_path, plot = g + theme(plot.title.position = "plot"), width = 8, height = 6, dpi = 300)
    
    cat("Forest plot saved for:", prs, "and", outcome, "\n")
  }
}

# Save results to CSV
output_file <- file.path(output_dir, "MVMA_RE_1step_lrm_trio_all_6-18y.csv")
write.csv(final_results, file = output_file, row.names = FALSE)

cat("Final results saved to", output_file, "\n")

# Save the meta-analysis models & statistics for future use
save(meta_results, file = file.path(output_dir, "MVMA_RE_1step_lrm_trio_all_6-18y.RData")) # Change file name according to the model and cohort

cat("Meta-analysis objects saved as 'MVMA_RE_1step_lrm_trio_all_6-18y.RData'.\n")


################################################################################
################### Save combined input data for future use ####################
################################################################################

# Save the combined data to a new CSV file
write.csv(meta_data, "lrm_trio_6-18y_m+f_combined.csv") # Change file name according to the model and cohort

cat("Combined data saved to 'lrm_trio_6-18y_m+f_combined.csv'.\n")

