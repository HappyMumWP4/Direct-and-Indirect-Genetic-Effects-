
###################################################################################
############################ LOAD NECESSARY LIBRARIES #############################
###################################################################################
# Install necessary packages if not already installed
# install.packages(c("car", "dplyr", "tidyr", "ggplot2"))

# Load required libraries
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)



###################################################################################
############## LOAD DATASET AND FILTER MISSING CASES IN THE DATASET ###############
###################################################################################

# Load necessary libraries
#library(dplyr)

#load data from fulfilled template xlsx file
filename <- "Template_PREDO_Trio dataset_PRScmf_WP4.xlsx" # change this to your own file
predo_dataset <- readxl::read_excel(filename, sheet = "dataset_3_wp4") # change this to your own datasheet name if necessary
head(predo_dataset) # change this to your own datasheet

# Filter the dataset to include only cases with no missing SDQ or CBCL scores
filtered_data <- predo_dataset %>%
  filter(!is.na("prs_c_") & !is.na("prs_m_") & !is.na("prs_f_"))
head(filtered_data)



###################################################################################
########### CALCULATE pTDT DEVIATION (pTDTd) USING PRSc, PRSm, AND PRSf ###########
###################################################################################

# calculate pTDTd
calculate_pTDT <- function(child_col, mother_col, father_col) {
  
  PRS_MP <- (filtered_data[[mother_col]] + filtered_data[[father_col]]) / 2
  pTDT_deviation <- (filtered_data[[child_col]] - PRS_MP) / sd(PRS_MP, na.rm = TRUE)
  
  return(pTDT_deviation)
}

# Create a new dataframe to store the pTDTds
filtered_data <- filtered_data %>%
  mutate(
    pTDTd_qc_AN_trio = calculate_pTDT("prs_c_qc_AN_trio", "prs_m_qc_AN_trio",
                                      "prs_f_qc_AN_trio"),
    pTDTd_qc_ADHD_trio = calculate_pTDT("prs_c_qc_ADHD_trio", "prs_m_qc_ADHD_trio",
                                        "prs_f_qc_ADHD_trio"),
    pTDTd_qc_ANX_trio = calculate_pTDT("prs_c_qc_ANX_trio", "prs_m_qc_ANX_trio",
                                       "prs_f_qc_ANX_trio"),
    pTDTd_qc_ASD_trio = calculate_pTDT("prs_c_qc_ASD_trio", "prs_m_qc_ASD_trio",
                                       "prs_f_qc_ASD_trio"),
    pTDTd_qc_BD_trio = calculate_pTDT("prs_c_qc_BD_trio", "prs_m_qc_BD_trio",
                                      "prs_f_qc_BD_trio"),
    pTDTd_qc_CDG_trio = calculate_pTDT("prs_c_qc_CDG_trio", "prs_m_qc_CDG_trio",
                                       "prs_f_qc_CDG_trio"),
    pTDTd_qc_EDU_trio = calculate_pTDT("prs_c_qc_EDU_trio", "prs_m_qc_EDU_trio",
                                       "prs_f_qc_EDU_trio"),
    pTDTd_qc_INSOM_trio = calculate_pTDT("prs_c_qc_INSOM_trio",
                                         "prs_m_qc_INSOM_trio",
                                         "prs_f_qc_INSOM_trio"),
    pTDTd_qc_MDD_trio = calculate_pTDT("prs_c_qc_MDD_trio", "prs_m_qc_MDD_trio",
                                       "prs_f_qc_MDD_trio"),
    pTDTd_qc_NEUROT_trio = calculate_pTDT("prs_c_qc_NEUROT_trio",
                                          "prs_m_qc_NEUROT_trio",
                                          "prs_f_qc_NEUROT_trio"),
    pTDTd_qc_OCD_trio = calculate_pTDT("prs_c_qc_OCD_trio", "prs_m_qc_OCD_trio",
                                       "prs_f_qc_OCD_trio"),
    pTDTd_qc_PPD_trio = calculate_pTDT("prs_c_qc_PPD_trio", "prs_m_qc_PPD_trio",
                                       "prs_f_qc_PPD_trio"),
    pTDTd_qc_PTSD_trio = calculate_pTDT("prs_c_qc_PTSD_trio", "prs_m_qc_PTSD_trio",
                                        "prs_f_qc_PTSD_trio"),
    pTDTd_qc_SCZ_trio = calculate_pTDT("prs_c_qc_SCZ_trio", "prs_m_qc_SCZ_trio",
                                       "prs_f_qc_SCZ_trio"),
    pTDTd_qc_ALC_trio = calculate_pTDT("prs_c_qc_ALC_trio", "prs_m_qc_ALC_trio",
                                       "prs_f_qc_ALC_trio"),
    pTDTd_qc_CIG_trio = calculate_pTDT("prs_c_qc_CIG_trio", "prs_m_qc_CIG_trio",
                                       "prs_f_qc_CIG_trio")
  )

# Rename the variable 'filtered_data'
data <- filtered_data
head(data)



###################################################################################
########################## OUTLIER WINSORIZATION ################################
###################################################################################

# Function to winsorize outliers
winsor_cols <- grep("prs_c_|prs_m_|prs_f_|pTDTd_", names(data), value = TRUE)
i = 0
winsorize_outlier <- function(data, columns = winsor_cols) {
  data_clean <- data
  for (col in columns) {
    Quantile1 <- quantile(data_clean[[col]], probs = .25, na.rm = TRUE)
    Quantile3 <- quantile(data_clean[[col]], probs = .75, na.rm = TRUE)
    IQR <- Quantile3 - Quantile1
    lower_bound <- Quantile1 - (IQR * 3)
    upper_bound <- Quantile3 + (IQR * 3)
    
    # Winsorize the column values in-place
    data_clean[[col]][data_clean[[col]] < lower_bound] <- lower_bound
    data_clean[[col]][data_clean[[col]] > upper_bound] <- upper_bound
    
    i = i+1
  }
  return(data_clean)
}

# The winsorized data is now stored in the variable 'data'
data <- winsorize_outlier(data)

# View the updated dataset
head(data)



###################################################################################
############################# pTDTd ONE-SAMPLE T-TSET #############################
###################################################################################

# Please notice that one-sample t-test should be done with nonstandardized pTDTd variables.
# List of nonstandardized pTDTd variables
prs_vars_pTDTd <- grep("pTDTd_qc_", names(data), value = TRUE)

# Function to perform the t-test and Shapiro-Wilk test
perform_tests <- function(data, variable) {
  t_test <- t.test(data[[variable]], mu = 0)
  shapiro_test <- shapiro.test(data[[variable]])
  
  results <- data.frame(
    Variable = variable,
    N = length(data[[variable]]),
    Shapiro_Wilk_p = shapiro_test$p.value,
    Mean = mean(data[[variable]], na.rm = TRUE),
    SD = sd(data[[variable]], na.rm = TRUE),
    t = t_test$statistic,
    df = t_test$parameter,
    Two_tailed_p = t_test$p.value,
    CI_lower = t_test$conf.int[1],
    CI_upper = t_test$conf.int[2]
  )
  
  return(results)
}

# Apply function to all pTDT variables and combine results
results_list <- lapply(prs_vars_pTDTd, function(var) perform_tests(data, var))
results_pTDTd <- do.call(rbind, results_list)

# View all results
head(results_pTDTd)

# Filter for significant results
significant_results_pTDTd <- results_pTDTd %>% filter(Two_tailed_p < 0.05)

# Print significant results
print(significant_results_pTDTd)

# Save the results to CSV files
write.csv(results_pTDTd, "one-sample_t_test_pTDTd_PREDO.csv", row.names = FALSE) # replace "PREDO" with your own cohort name
write.csv(significant_results_pTDTd, "one-sample_t_test_significant_pTDTd_PREDO.csv", row.names = FALSE) # replace "PREDO" with your own cohort name

# Function to generate and save box plots for significant results
generate_boxplots <- function(variable, data) {
  
  if (!dir.exists("pTDTd_Boxplots")) {
    dir.create("pTDTd_Boxplots")
  }
  
  # Create the box plot 
  box_plot <- ggplot(data, aes_string(x = "factor(1)", y = variable)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    ggtitle(paste("Box plot for", variable)) +
    theme_minimal() +
    labs(x = "", y = variable)  
  
  print(box_plot)
  
  # Define the filename
  plot_filename <- paste("pTDTd_Boxplots", paste(variable, "boxplot_PREDO.png", sep="_"), sep="/") # replace "PREDO" with your own cohort name
  
  # Save the box plot to a png file
  ggsave(plot_filename, plot = box_plot, width = 7, height = 5, dpi = 300)
}

# Apply the function 
lapply(significant_results_pTDTd$Variable, generate_boxplots, data = data)

# Completion message 
print("One-sample Ttest csv file and box plots for pTDTd have been saved.")