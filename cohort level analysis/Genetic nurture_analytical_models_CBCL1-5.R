## title: "Genetic nurture_analytical_models_Child_1-5y"
## authors: "Mina Shahisavandi: m.shahisavandi@erasmusmc.nl; Li Tian: li.tian@helsinki.fi"
## date: "2024-08-18"

## Note
# Scripts in this R file are meant for analyzing a CBCL1_5 cohort. They can be adapted to analyze other similar cohort after changing definitions of corresponding covariates and child behavioral outcome variables.

# Scripts for exporting csv files in "JOINT LRM", "JOINT LOGISTIC" and "ANCOVA" sections are "commented-out" currently and replaced by a joint xlsx file, respectively.


###################################################################################
############################ LOAD NECESSARY LIBRARIES #############################
###################################################################################
# List of required packages
required_packages <- c("car", "dplyr", "tidyr", "psych", "FactoMineR", "factoextra",
                       "nortest", "glmnet", "ppcor", "purrr", "lmtest", "tidyverse",
                       "broom", "caret", "pROC", "multcomp", "ggplot2", "writexl")

# Install missing packages if needed
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# Load required libraries
library(car)
library(dplyr)
library(tidyr)
library(psych)
library(FactoMineR)
library(factoextra)
library(nortest)
library(glmnet)
library(ppcor)
library(purrr)
library(lmtest)
library(tidyverse)
library(broom)
library(caret)
library(pROC)
library(multcomp)
library(writexl)
library(ggplot2)



###################################################################################
############## LOAD DATASET AND FILTER MISSING CASES IN THE DATASET ###############
###################################################################################

# Load necessary libraries
#library(dplyr)

#load data from fulfilled template xlsx file
filename <- "Template_PREDO_Trio dataset_CBCL1-6y_WP4.xlsx" # change this to your own file
predo_dataset <- readxl::read_excel(filename, sheet = "dataset_1_wp4") # change this to your own datasheet
head(predo_dataset)

# Filter the dataset to include only cases with no missing SDQ or CBCL scores
filtered_data <- predo_dataset %>%
  filter(!is.na(h_c_intern_1_5y) & !is.na(h_c_extern_1_5y) & !is.na(h_c_total_1_5y))

# View the updated dataset
head(filtered_data)





###################################################################################
########################## OUTLIER WINSORIZATION_#1 ###############################
###################################################################################

# Load necessary libraries
#library(dplyr)

# Function to winsorize outliers in CBCL1_5 (adjust these to your dataset's variable names if necessary)
winsor_cols <- grep("h_c_sex|h_c_total_cat_1_5y", names(data), value = TRUE, invert = TRUE)
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

# Replace with the winsorized variables in the dataset
data <- winsorize_outlier(data)

# View the updated dataset
head(data)




###################################################################################
################### MEAN SCORE OF 13 PSYCHIATRIC PRS VARIABLES ####################
###################################################################################

# Load necessary libraries
#library(dplyr)
#library(psych)
#library(FactoMineR)
#library(factoextra)

# Define the 13 psychiatric PRS variables for CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_13 <- c("MDD", "BD", "SCZ", "ADHD", "AN", "ASD", "NEUROT", "ANX", "OCD", "CDG", "PPD", "PTSD", "INSOM")
prs_vars_13 <- paste0("prs_c_qc_", prs_13, "_trio")

# Function to calculate mean values and add to dataset
calculate_and_add_mean <- function(data, prs_vars_13, suffix) {
  mean_value <- rowMeans(data[prs_vars_13], na.rm = TRUE)
  data[[paste0("", suffix)]] <- mean_value
  data
}

# Calculate mean values and add them for PRSc, PRSm, PRSf, PRSnt, and pTDTd for CBCL1_5

# PRSc
data <- calculate_and_add_mean(data, prs_vars_13, "prs_c_qc_psych_mean_trio") # adjust this variable name if necessary

# PRSm
prs_vars_13 <- sub("prs_c_qc_", "prs_m_qc_", prs_vars_13)
data <- calculate_and_add_mean(data, prs_vars_13, "prs_m_qc_psych_mean_trio") # adjust this variable name if necessary

# PRSf
prs_vars_13 <- sub("prs_m_qc_", "prs_f_qc_", prs_vars_13)
data <- calculate_and_add_mean(data, prs_vars_13, "prs_f_qc_psych_mean_trio") # adjust this variable name if necessary


# View the updated dataset with new variables added
head(data)



###################################################################################
######################## MEAN SCORE OUTLIER WINSORIZATION #########################
###################################################################################

# Load necessary libraries
#library(dplyr)

# Function to winsorize outliers in mean scores (adjust to your dataset's variable names if necessary)
winsor_cols_2 <- grep("psych_mean", names(data), value = TRUE)
i = 0
winsorize_outlier <- function(data, columns = winsor_cols_2) {
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

# Replace with the winsorized mean scores in the dataset
data <- winsorize_outlier(data)

# View the updated dataset
head(data)



###################################################################################
###################### STANDARDIZE VARIABLES IN THE DATASET #######################
###################################################################################

# Define a function to standardize all the continuous variables
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Identify the columns to be standardized in CBCL1_5 (adjust these to your dataset's variable names if necessary)
stand_cols <- grep("h_c_sex|h_c_total_cat_1_5y", names(data), value = TRUE, invert = TRUE)

# Apply the standardization function to the identified columns
data[stand_cols] <- lapply(data[stand_cols], standardize) 

# View the standardized dataset
head(data)



###################################################################################
########## PRINCIPAL COMPONENT ANALYSIS OF 13 PSYCHIATRIC PRS VARIABLES ###########
###################################################################################

# Load necessary libraries
#library(dplyr)
#library(psych)
#library(FactoMineR)
#library(factoextra)

# Define the 13 psychiatric PRS variables for CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_13 <- c("MDD", "BD", "SCZ", "ADHD", "AN", "ASD", "NEUROT", "ANX", "OCD", "CDG", "PPD", "PTSD", "INSOM")
prs_vars_13 <- paste0("prs_c_qc_", prs_13, "_trio")

# Function to perform PCA and extract the first principal component
perform_pca_and_add_pc1 <- function(data, prs_vars_13, suffix) {
  prs_data_13 <- data[prs_vars_13]
  
  # Perform PCA
  pca_result <- PCA(prs_data_13, graph = FALSE)
  
  # Extract the first principal component
  pc1 <- pca_result$ind$coord[, 1]
  
  # Add the first principal component to the dataset
  data[[paste0("", suffix)]] <- pc1
  
  # Return the updated dataset
  data
}

# Perform PCA and add the first PC for PRSc, PRSm, PRSf, PRSnt, and pTDTd for CBCL1_5

# PRSc
data <- perform_pca_and_add_pc1(data, prs_vars_13, "prs_c_qc_psych_pc1_trio") # adjust this variable name if necessary

# PRSm
prs_vars_13 <- sub("prs_c_qc_", "prs_m_qc_", prs_vars_13)
data <- perform_pca_and_add_pc1(data, prs_vars_13, "prs_m_qc_psych_pc1_trio") # adjust this variable name if necessary

# PRSf
prs_vars_13 <- sub("prs_m_qc_", "prs_f_qc_", prs_vars_13)
data <- perform_pca_and_add_pc1(data, prs_vars_13, "prs_f_qc_psych_pc1_trio") # adjust this variable name if necessary


# View the updated dataset
head(data)



###################################################################################
################### PC1 OUTLIER WINSORIZATION & STANDARDIZATION ###################
###################################################################################

# Load necessary libraries
#library(dplyr)

# Function to winsorize outliers in PC1 scores (adjust to your dataset's variable names if necessary)
winsor_cols_3 <- grep("psych_pc1", names(data), value = TRUE)
i = 0
winsorize_outlier <- function(data, columns = winsor_cols_3) {
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

# Replace with the winsorized PC1 variables in the dataset
data <- winsorize_outlier(data)

# Define a function to standardize the PC1 variables
standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Identify the columns to be standardized (adjust to your dataset's variable names if necessary)
stand_cols_2 <- grep("psych_pc1", names(data), value = TRUE)

# Apply the standardization function to the identified columns
data[stand_cols_2] <- lapply(data[stand_cols_2], standardize)

# View the updated dataset
head(data)

# Write data to a csv file
# write.csv(data, "winsored_standardized_trio_dataset_cbcl1_5_wp4.csv") # adjust file name to specify your own cohort



###################################################################################
################# PRINCIPAL COMPONENT ANALYSIS OF ALL 16 PRSs #####################
###################################################################################

# Load necessary libraries
#library(dplyr)
#library(psych)
#library(FactoMineR)
#library(factoextra)
#library(ggplot2)

# Define the 16 PRS variables for CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_16 <- c("MDD", "BD", "SCZ", "ADHD", "AN", "ASD", "NEUROT", "ANX", "OCD", "CDG", "PPD", "PTSD", "INSOM", "EDU", "ALC", "CIG")
prs_vars_16 <- paste0("prs_c_qc_", prs_16, "_trio")

# Function to perform PCA
perform_pca <- function(data, prs_vars_16) {
  prs_data_16 <- data[prs_vars_16]
  
  # Perform PCA
  pca_results <- PCA(prs_data_16, graph = FALSE)
  
  # Return PCA results
  return(pca_results)
}

# Perform PCA for PRSc, PRSm, PRSf, PRSnt, and pTDTd for CBCL1_5 and get PCA results

all_pca_results <- list()

# PRSc
pca_results <- perform_pca(data, prs_vars_16)
all_pca_results$PRSc <- pca_results$var$contrib

# PRSm
prs_vars_16 <- sub("prs_c_qc_", "prs_m_qc_", prs_vars_16) # adjust these to your dataset's variable names if necessary
pca_results <- perform_pca(data, prs_vars_16)
all_pca_results$PRSm <- pca_results$var$contrib

# PRSf
prs_vars_16 <- sub("prs_m_qc_", "prs_f_qc_", prs_vars_16) # adjust these to your dataset's variable names if necessary
pca_results <- perform_pca(data, prs_vars_16)
all_pca_results$PRSf <- pca_results$var$contrib


# View results
all_pca_results 

# View and write PCA results to a csv file
write.csv(all_pca_results, "pca_16_prs_contrib_cbcl1_5.csv") # replace cbcl1_5 to specify your own cohort

# Create a directory for the plots
if (!dir.exists("PRS_16_PCA_Barplots_CBCL1_5")) {
    dir.create("PRS_16_PCA_Barplots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Print and save PC1 plots
print("Barplots of PCA_16-PRSs:")

for (result_name in names(all_pca_results)) {
  result <- all_pca_results[[result_name]]
  dim1_col <- data.frame(result[,1])
  dim1_col <- tibble::rownames_to_column(dim1_col, "name")
  colnames(dim1_col) <- c('name', 'value')
  plot <- ggplot((dim1_col), aes(x = reorder(name, value), y = value)) +
    geom_bar(stat = "identity") +
    labs(title = paste0("bar_plot_of_16_prs_pc1_", result_name),
         x = "PRS contribution") +
    coord_flip()
  
  # Print plots
  print(plot)
  
  # Save PC1 plots to the folder
  plot_filename <- paste0("PRS_16_PCA_Barplots_CBCL1_5/", plot$labels$title, "_cbcl1_5", ".png") # replace cbcl1_5 to specify your own cohort
  ggsave(plot_filename, plot = plot, width = 8, height = 6, dpi = 300)
}

# Completion message
print("Bar plots for PC1 load of 16 PRSs in PCA have been saved.")



###################################################################################
####################### CHECK NORMALITY AND CREATE QQ PLOTS #######################
###################################################################################

# Load necessary libraries
#library(car)
#library(nortest)
#library(tidyr)
#library(dplyr)
#library(glmnet)

# Initialize a dataframe where to store Shapiro–Wilk's test results
norm_test <- data.frame(variable = "", pval = "")

# Initialize a variable used to iterate through the "norm_test" dataframe in CBCL1_5
a=1
data_cols <- grep("h_c_sex|h_c_total_cat_1_5y", names(data), value = TRUE, invert = TRUE) # adjust these to your dataset's variable names if necessary

# Initialize a directory for Q-Q plots
plots_dir <- "Normality_QQ_Plots_CBCL1_5"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir) # replace CBCL1_5 in the directory name to specify your own cohort
}

for(i in data_cols) {
  # Print(paste0("Processing variable ", i))
  # Define continuous variables as numeric
  data[[i]] = as.numeric(data[[i]])
  
  # create and save QQ plots
  jpeg(paste(plots_dir, "/","qq_plot_normality_cbcl1_5_", toString(i),".jpeg", sep = ""), width = 400, height = 400) # replace cbcl1_5 to specify your own cohort
  qqnorm(data[[i]])
  grid()
  qqline(data[[i]],lwd=2, col="red")
  dev.off()

  # Test normality with Shapiro-Wilk's test
    sw <- shapiro.test(data[[i]])
    norm_result <- sw
  
  # Store normality test results to in "norm_test" dataframe
  norm_test[a, "variable"] = toString(i)
  norm_test[a, "pval"] = norm_result$p.value
  
  a=a+1
}

# save normality results to a csv file
write.csv(norm_test, "normality_results_cbcl1_5.csv") # replace cbcl1_5 to specify your own cohort

# View normality results
head(norm_test)

# Completion message
print("Normality QQ plots for continuous variables have been saved.")



###################################################################################
############################# CHECK MULTICOLLINEARITY #############################
###################################################################################

# Load necessary libraries
#library(car)
#library(tidyr)
#library(dplyr)
#library(glmnet)
#library(writexl)

# Define covariates for CBCL1_5 (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Define predictors for CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs = c("MDD", "BD", "SCZ", "ADHD", "AN", "ASD", "NEUROT", "ANX", "OCD", "CDG", "PPD", "PTSD", "INSOM", "EDU", "ALC", "CIG")
vif_c <- paste0("prs_c_qc_", prs, "_trio")
vif_m <- paste0("prs_m_qc_", prs, "_trio")
vif_f <- paste0("prs_f_qc_", prs, "_trio")

vifs = list(vif_c, vif_m, vif_f)

vif_trio = c("prs_c_qc_", "prs_m_qc_", "prs_f_qc_")

prs = c("MDD", "BD", "SCZ", "ADHD", "AN", "ASD", "NEUROT", "ANX", "OCD", "CDG", "PPD", "PTSD", "INSOM", "EDU", "ALC", "CIG", "psych_pc1", "psych_mean")

for (i in prs){
  vifs[[length(vifs) +1]] = (paste0(vif_trio, i, "_trio"))
  
}

names(vifs) = c("Joint_child", "Joint_mother", "Joint_father", "Joint_nt", "Joint_pTDTd", prs)

# Initialize a directory for Q-Q plots
plots_dir <- "Multicollinearity_VIFs_CBCL1_5"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir) # replace CBCL1_5 in the directory name to specify your own cohort
}

# Make multicollinearity analysis of PRSs
i = 1
for (predictors_qc in vifs) {
  predictors_qc = c(covariates, predictors_qc)
  name = names(vifs)[i]
  i = i+1
  
  sub_HM_T4 <- data[, predictors_qc]
  # sub_HM_T4 <- sub_HM_T4[complete.cases(sub_HM_T4), ] 
  
  vif_dataframe_qc <- data.frame(var = character(0), predictor = character(0), VIF = numeric(0))
  
  # Loop through each variable
  for (var in predictors_qc) {
    # print(name)
    # print(var)
    # 
    # Exclude the current variable from the list of predictors and adapt the variables names
    predictors <- predictors_qc[predictors_qc != var]
    # print(predictors)
    predictors <- paste0("sub_HM_T4$`", predictors, '`')
    var_name <- paste0("sub_HM_T4$`", var, "`")
    
    # Construct the formula for linear regression
    formula <- as.formula(paste(var_name, '~', paste(predictors, collapse = "+")))
    
    # Fit the linear regression model
    model <- lm(formula, data = sub_HM_T4)
    
    # Calculate VIF
    vif_results <- vif(model)
    
    # Create dataframe for current variable
    var_dataframe_qc <- data.frame(var = var, predictor = names(vif_results), VIF = vif_results)
    
    # Bind the dataframe to the main dataframe
    vif_dataframe_qc <- rbind(vif_dataframe_qc, var_dataframe_qc)
  }
  
  # Spread the dataframe to pivot it and write to a csv file
  vif_matrix_qc <- spread(vif_dataframe_qc, var, VIF)
  vif_matrix_qc$predictor <- NULL
  rownames(vif_matrix_qc) <- colnames(vif_matrix_qc)
  write.table(vif_matrix_qc, paste0("Multicollinearity_VIFs_CBCL1_5/vif_qc_prs_trio_cbcl1_5_", name, ".csv"), sep = '\t') # replace cbcl1_5 to specify your own cohort
  
  # If you get the following error "there are aliased coefficients in the model", it is possible that some variables are really highly or perfectly correlated
  # Therefore, try the following code to get the correlation matrix and check whether it is true
  # corr_matrix <- cor(sub_HM_T4, use="complete.obs")
  # write.table(corr_matrix, paste0("corr_matrix_results_", name, ".csv"), quote= FALSE, row.names=FALSE, sep=";")
}

# Completion message
print("Multicollinearity csv files of PRSs have been saved.")



###################################################################################
############################### PARTIAL CORRELATION ###############################
###################################################################################

# Load necessary libraries
#library(dplyr)
#library(psych)
#library(ppcor)
#library(ggplot2)

# Define covariates for CBCL1_5 (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Define prs predictors for CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_vars_c <- grep("prs_c_qc_", names(data), value = TRUE)
prs_vars_m <- grep("prs_m_qc_", names(data), value = TRUE)
prs_vars_f <- grep("prs_f_qc_", names(data), value = TRUE)
prs_vars <- c(prs_vars_c, prs_vars_m, prs_vars_f)

# Define CBCL Tscore variables for CBCL1_5 (adjust these to your dataset's variable names if necessary)
cbcl_vars <- c("h_c_intern_1_5y", "h_c_extern_1_5y", "h_c_total_1_5y")

# Function to perform partial correlation and calculate descriptive statistics
perform_partial_correlation <- function(data, prs_vars, cbcl_vars, covariates) {
  formula <- as.formula(paste(cbcl_vars, "~", prs_vars, "+", paste(covariates, collapse = "+")))
  
  # Calculate partial correlation
  partial_corr <- pcor.test(data[[prs_vars]], data[[cbcl_vars]], data[covariates], method = "spearman")
  
  # Calculate descriptive statistics
  n <- sum(!is.na(data[[prs_vars]]) & !is.na(data[[cbcl_vars]]))
  shapiro_p_prs <- shapiro.test(data[[prs_vars]])$p.value
  shapiro_p_cbcl <- shapiro.test(data[[cbcl_vars]])$p.value
  mean_prs <- mean(data[[prs_vars]], na.rm = TRUE)
  sd_prs <- sd(data[[prs_vars]], na.rm = TRUE)
  mean_cbcl <- mean(data[[cbcl_vars]], na.rm = TRUE)
  sd_cbcl <- sd(data[[cbcl_vars]], na.rm = TRUE)
  
  list(
    n = n,
    shapiro_p_prs = shapiro_p_prs,
    shapiro_p_cbcl = shapiro_p_cbcl,
    mean_prs = mean_prs,
    sd_prs = sd_prs,
    mean_cbcl = mean_cbcl,
    sd_cbcl = sd_cbcl,
    r = partial_corr$estimate,
    p_value = partial_corr$p.value
  )
}

# Initialize an empty list to store results
results <- list()

# Loop through each combination of PRS and CBCL variables and perform partial correlation
for (prs_var in prs_vars) {
  for (cbcl_var in cbcl_vars) {
    result <- perform_partial_correlation(data, prs_var, cbcl_var, covariates)
    results <- append(results, list(list(prs_var = prs_var, cbcl_var = cbcl_var, result = result)))
  }
}

# Extract significant correlations and create scatter plots
significant_results <- results %>% 
  purrr::keep(~ .x$result$p_value < 0.05)

# Create scatterplots for significant correlations
plot_scatter <- function(data, prs_var, cbcl_var) {
  ggplot(data, aes_string(x = prs_var, y = cbcl_var)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste("scatter_plot_of_", prs_var, "vs", cbcl_var),
         x = prs_var, y = cbcl_var) +
    theme_minimal()
}

plots <- lapply(significant_results, function(x) {
  plot_scatter(data, x$prs_var, x$cbcl_var)
})

# Print plots
print("Scatterplots for Significant Partial Correlations in CBCL1_5:")
plots

# Create a directory for the plots
if (!dir.exists("Partial_Correlation_PRS_Scatterplots_CBCL1_5")) {
    dir.create("Partial_Correlation_PRS_Scatterplots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Save scatter plots for significant results
for (i in plots) {
  plot_filename <- paste0("Partial_Correlation_PRS_Scatterplots_CBCL1_5/", i$labels$title, ".png") # replace cbcl1_5 to specify your own cohort
  ggsave(plot_filename, plot = i, width = 6, height = 4, dpi = 300)
}

# Save significant correlation results to a CSV file
significant_results_partcorr <- do.call(rbind, lapply(significant_results, function(x) {
  data.frame(
    PRS_Var = x$prs_var,
    CBCL_Var = x$cbcl_var,
    N = x$result$n,
    Shapiro_p_PRS = x$result$shapiro_p_prs,
    Shapiro_p_CBCL = x$result$shapiro_p_cbcl,
    Mean_PRS = x$result$mean_prs,
    SD_PRS = x$result$sd_prs,
    Mean_CBCL = x$result$mean_cbcl,
    SD_CBCL = x$result$sd_cbcl,
    Spearman_r = x$result$r,
    Spearman_p_value = x$result$p_value
  )
}))
write.csv(significant_results_partcorr, "partial_correlation_significant_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# View significant correlation results
head(significant_results_partcorr)

# Completion message
print("Partial correlation csv file and scatter plots for significant variables have been saved.")



###################################################################################
############ ADJUSTED LINEAR REGRESSION MODEL (LRM) WITH PRS, PC1, MEAN ###########
###################################################################################

# Load necessary libraries
#library(car)
#library(lmtest)
#library(purrr)
#library(ggplot2)

# Define PRS variables in CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_vars_c <- grep("prs_c_qc_", names(data), value = TRUE)
prs_vars_m <- grep("prs_m_qc_", names(data), value = TRUE)
prs_vars_f <- grep("prs_f_qc_", names(data), value = TRUE)
rs_vars <- c(prs_vars_c, prs_vars_m, prs_vars_f)

# Define CBCL1_5 outcome variables (adjust these to your dataset's variable names if necessary)
cbcl_vars <- c("h_c_intern_1_5y", "h_c_extern_1_5y", "h_c_total_1_5y")

# Define covariates in CBCL1_5 (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Function to perform linear regression and calculate diagnostics
perform_regression <- function(data, prs_var, cbcl_var, covariates) {
  formula <- as.formula(paste(cbcl_var, "~", prs_var, "+", paste(covariates, collapse = "+")))
  
  # Fit the linear model
  model <- lm(formula, data = data)
  
  # Extract coefficients
  coef_matrix <- summary(model)$coefficients
  coef_df <- as.data.frame(t(coef_matrix[2, ]))  # Transpose to ensure a single row data frame
  colnames(coef_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coef_df$lower_ci <- coef_df$Estimate - 1.96 * coef_df$`Std. Error`
  coef_df$upper_ci <- coef_df$Estimate + 1.96 * coef_df$`Std. Error`
  coef_df$Predictor <- rownames(coef_matrix)[2]

  # Create Q-Q plot of standardized residuals
  qq_plot <- ggplot(data = data.frame(residuals = scale(residuals(model))), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = paste("Normal Q-Q Plot of Standardized Residuals for ", prs_var, " on ", cbcl_var),
         x = "Theoretical Quantiles", y = "Standardized Residuals") +
    theme_minimal()

  # Create coefficient plot
  coef_plot <- ggplot(coef_df, aes(x = Predictor, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    theme_minimal() +
    coord_flip() +  # Flip the coordinates for better readability
    ggtitle(paste("Coefficient Plot of LRM for Single", prs_var, "on", cbcl_var)) +
    ylab("Estimate") +
    xlab("Predictor")
  
  # Create a list of statistic results
  list(
    model_summary = summary(model),
    n = nobs(model),
    f_value = summary(model)$fstatistic[1],
    p_value_total = pf(summary(model)$fstatistic[1], summary(model)$fstatistic[2], summary(model)$fstatistic[3], lower.tail = FALSE),
    mean_value = mean(data[[cbcl_var]], na.rm = TRUE),
    se = coef_df$`Std. Error`,
    df = model$df.residual,
    beta = coef_df$Estimate,
    p_value = coef_df$`Pr(>|t|)`,
    ci_lower = coef_df$lower_ci,
    ci_upper = coef_df$upper_ci,
    r_squared = summary(model)$r.squared,
    vif_values = car::vif(model),
    constant_variance = bptest(model)$p.value,
    normality_residuals = shapiro.test(residuals(model))$p.value,
    durbin_watson = lmtest::dwtest(model)$statistic[1],
    qq_plot = qq_plot,
    coef_plot = coef_plot
  )
}

# Perform linear regression for prs in CBCL1_5
results_lrm <- list()
for (prs_var in prs_vars) {
  for (cbcl_var in cbcl_vars) {
    result <- perform_regression(data, prs_var, cbcl_var, covariates)
    results_lrm <- append(results_lrm, list(list(prs_var = prs_var, cbcl_var = cbcl_var, result = result)))
    }
  }

# Save all linear regression results into a data frame
results_df_lrm <- do.call(rbind, lapply(results_lrm, function(x) {
  data.frame(
    CBCL_outcome = x$cbcl_var,
    PRS_predictor = x$prs_var,
    N = x$result$n,
    F_value = x$result$f_value,
    P_value_total = x$result$p_value_total,
    Mean_value = x$result$mean_value,
    SE = x$result$se,
    df = x$result$df,
    Beta = x$result$beta,
    P_value = x$result$p_value,
    CI_lower = x$result$ci_lower,
    CI_upper = x$result$ci_upper,
    R_squared = x$result$r_squared,
    # vif_values = paste(x$result$vif_values, collapse = ";"), # convert VIF values to a single string
    Durbin_Watson = x$result$durbin_watson,
    Constant_variance = x$result$constant_variance,
    Normality_residuals = x$result$normality_residuals,
    
    row.names = NULL  # This ensures that row names are not automatically set
  )
}))

# View all linear regression results
head(results_df_lrm)

# Save all linear regression results to a CSV file
write.csv(results_df_lrm, "lrm_prs_summary_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Extract and print significant linear regression results
significant_results_df_lrm <- subset(results_df_lrm, P_value < 0.05)
print(significant_results_df_lrm)

# Save significant linear regression results to a CSV file
write.csv(significant_results_df_lrm, "lrm_prs_significant_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Create directories for Q-Q and Coefficient plots
if (!dir.exists("LRM_Single_PRS_QQ_Plots_CBCL1_5")) {
    dir.create("LRM_Single_PRS_QQ_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LRM_Single_PRS_Coeff_Plots_CBCL1_5")) {
  dir.create("LRM_Single_PRS_Coeff_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Print and save Q-Q plots for significant results
print("Q-Q Plots of LRM Single PRS Significant Results for CBCL1_5:")
print("Coefficient Plots of LRM Single PRS Significant Results for CBCL1_5:")

# Plot and save png files
plot_and_save_qq <- function(predictor_name, outcome_name, result) {
  png(paste0("LRM_Single_PRS_QQ_Plots_CBCL1_5/qq_single_", predictor_name, "_", outcome_name, ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(result$qq_plot)
  dev.off()
  plot(result$qq_plot)
}

plot_and_save_coef <- function(predictor_name, outcome_name, result) {
  png(paste0("LRM_Single_PRS_Coeff_Plots_CBCL1_5/lrm_coeff_single_", predictor_name, "_", outcome_name, ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(result$coef_plot)
  dev.off()
  plot(result$coef_plot)
}

# Loop through results and plot only significant ones
for (i in seq_along(results_lrm)) {
  if (results_lrm[[i]]$result$model_summary$coefficients[2, 4] < 0.05) {
    # Define predictor and outcome variable names in png files
    predictor_name <- gsub(" ", "_", results_lrm[[i]]$prs_var)
    outcome_name <- gsub(" ", "_", results_lrm[[i]]$cbcl_var)
    
    # Perform plotting function
    plot_and_save_qq(predictor_name, outcome_name, results_lrm[[i]]$result)
    plot_and_save_coef(predictor_name, outcome_name, results_lrm[[i]]$result)
  }
}

# Completion message
print("LRM CSV files, QQ and Coefficient plots for significant variables have been saved.")




###################################################################################
########################### TRIO-PRS, -PC1, -MEAN LRM #############################
###################################################################################

# Load required packages
#library(car)
#library(purrr)
#library(lmtest)
#library(ggplot2)

# Define CBCL1_5 outcome variables (adjust these to your dataset's variable names if necessary)
cbcl_vars <- c("h_c_intern_1_5y", "h_c_extern_1_5y", "h_c_total_1_5y")
  
# Define covariates in CBCL1_5 (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Define variables for PRS child, mother, father in CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_vars_c <- grep("prs_c_qc_", names(data), value = TRUE)
prs_vars_m <- grep("prs_m_qc_", names(data), value = TRUE)
prs_vars_f <- grep("prs_f_qc_", names(data), value = TRUE)

# Function to perform Trio-LRM analysis and calculate diagnostics
perform_trio_regression <- function(data, prs_child_var, prs_mother_var, prs_father_var, cbcl_var, covariates) {
  formula <- as.formula(paste(cbcl_var, "~", prs_child_var, "+", prs_mother_var, "+", prs_father_var, "+", paste(covariates, collapse = "+")))
  
  # Fit the linear model
  model <- lm(formula, data = data)
  
  # Calculate diagnostics
  n <- nobs(model)
  f_value <- summary(model)$fstatistic[1]
  df1 <- summary(model)$fstatistic[2]
  df2 <- summary(model)$fstatistic[3]
  p_value_total <- pf(f_value, df1, df2, lower.tail = FALSE)
  mean_value <- mean(data[[cbcl_var]], na.rm = TRUE)
  se <- summary(model)$coefficients[2, 2]
  df <- model$df.residual
  betas <- summary(model)$coefficients[2:4, 1]
  p_values <- summary(model)$coefficients[2:4, 4]
  ci <- confint(model, level = 0.95)[2:4, ]
  r_squared <- summary(model)$r.squared
  vif_values <- vif(model)
  residuals <- residuals(model)
  durbin_watson <- dwtest(model)$statistic
  constant_variance <- bptest(model)$p.value
  normality_residuals <- shapiro.test(residuals)$p.value
  
  # Create Q-Q plot of standardized residuals
  qq_plot <- ggplot(data = data.frame(residuals = scale(residuals)), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = paste("Normal Q-Q Plot of Standardized Residuals for", prs_child_var, prs_mother_var, prs_father_var, "on", cbcl_var),
         x = "Theoretical Quantiles", y = "Standardized Residuals") +
    theme_minimal()

  # Create coefficients plot
  coef_df <- as.data.frame(summary(model)$coefficients[2:4,])
  colnames(coef_df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  coef_df$lower_ci <- coef_df$Estimate - 1.96 * coef_df$`Std. Error`
  coef_df$upper_ci <- coef_df$Estimate + 1.96 * coef_df$`Std. Error`
  coef_df$Predictor <- rownames(coef_df)
  coef_plot <- ggplot(coef_df, aes(x = Predictor, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    theme_minimal() +
    coord_flip() +  # Flip the coordinates for better readability
    ggtitle(paste("Coefficient Plot of Trio-LRM for", prs_child_var, prs_mother_var, prs_father_var, "on", cbcl_var)) +
    ylab("Estimate") +
    xlab("Predictors")

  # Create a list of statistic results
  list(
    model_summary = summary(model),
    n = n,
    f_value = f_value,
    p_value_total = p_value_total,
    mean_value = mean_value,
    se = se,
    df = df,
    betas = betas,
    p_values = p_values,
    ci = ci,
    r_squared = r_squared,
    vif_values = vif_values,
    constant_variance = constant_variance,
    normality_residuals = normality_residuals,
    durbin_watson = durbin_watson,
    qq_plot = qq_plot,
    coef_plot = coef_plot
  )
}

# Perform Trio-LRM analysis for CBCL1_5
results_lrm_trio <- list()
for (cbcl_var in cbcl_vars) {
  for (i in 1:length(prs_vars_c)) {
    result <- perform_trio_regression(data, prs_vars_c[i], prs_vars_m[i], prs_vars_f[i], cbcl_var, covariates)
    results_lrm_trio <- append(results_lrm_trio, list(list(cbcl_var = cbcl_var, prs_child_var = prs_vars_c[i], prs_mother_var = prs_vars_m[i], prs_father_var = prs_vars_f[i], result = result)))
  }
}

# Save all linear regression results to a data frame
results_df_lrm_trio <- do.call(rbind, lapply(results_lrm_trio, function(x) {
  data.frame(
    CBCL_outcome = x$cbcl_var,
    PRS_predictor = names(x$result$p_values),
    N = x$result$n,
    F_value = x$result$f_value,
    P_value_total = x$result$p_value_total,
    Mean_value = x$result$mean_value,
    SE = x$result$se,
    df = x$result$df,
    Beta = x$result$betas,
    P_value = x$result$p_values,
    CI_lower = x$result$ci[, 1],
    CI_upper = x$result$ci[, 2],
    R_squared = x$result$r_squared,
    Durbin_Watson = x$result$durbin_watson,
    Constant_variance = x$result$constant_variance,
    Normality_residuals = x$result$normality_residuals,
    
    row.names = NULL  # This ensures that row names are not automatically set
  )
}))

# View all linear regression results
head(results_df_lrm_trio)

# Save all linear regression results to a CSV file
write.csv(results_df_lrm_trio, "lrm_trio_PRScmf_summary_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Extract and print significant linear regression results
significant_results_df_lrm_trio <- subset(results_df_lrm_trio, P_value < 0.05)
print(significant_results_df_lrm_trio)

# Save significant linear regression results to a CSV file
write.csv(significant_results_df_lrm_trio, "lrm_trio_PRScmf_significant_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Create directories for Q-Q and Coefficient plots
if (!dir.exists("LRM_Trio_PRScmf_QQ_Plots_CBCL1_5")) {
  dir.create("LRM_Trio_PRScmf_QQ_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LRM_Trio_PRScmf_Coeff_Plots_CBCL1_5")) {
  dir.create("LRM_Trio_PRScmf_Coeff_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Print and save Q-Q and Coefficient plots for significant results
print("Q-Q Plots of LRM Trio PRScmf Significant Results for CBCL1_5:")
print("Coefficient Plots of LRM Trio PRScmf Significant Results for CBCL1_5:")

# Plot and save png files
plot_and_save_qq <- function(predictor_name, outcome_name, result) {
  png(paste0("LRM_Trio_PRScmf_QQ_Plots_CBCL1_5/qq_PRScmf_", predictor_name, "_", outcome_name, ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(result$qq_plot)
  dev.off()
  plot(result$qq_plot)
}

plot_and_save_coef <- function(predictor_name, outcome_name, result) {
  png(paste0("LRM_Trio_PRScmf_Coeff_Plots_CBCL1_5/lrm_coeff_PRScmf_", predictor_name, "_", outcome_name, ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(result$coef_plot)
  dev.off()
  plot(result$coef_plot)
}

# Loop through results and plot only significant ones
for (i in seq_along(results_lrm_trio)) {
  # Extract the current result item
  current_result <- results_lrm_trio[[i]]
  
  # Extract the p-values for child, mother, and father PRS
  p_values <- current_result$result$model_summary$coefficients[2:4, 4]

  # List of PRS variables
  prs_vars <- c(current_result$prs_child_var, current_result$prs_mother_var, current_result$prs_father_var)
  
  # Loop through each PRS variable (child, mother, father) and check significance
  for (j in seq_along(prs_vars)) {
    prs_var <- prs_vars[j]
    
    # Check if this PRS variable is significant
    if (p_values[j] < 0.05) {
      # Define predictor and outcome variable names in png files
      predictor_name <- prs_var
      outcome_name <- gsub(" ", "_", current_result$cbcl_var)
      
      # Perform plotting function using only rows 2 to 4 of the coefficients matrix
      plot_and_save_qq(predictor_name, outcome_name, current_result$result)
      plot_and_save_coef(predictor_name, outcome_name, current_result$result)
    }
  }
}

# Completion message
print("Trio_LRM csv files, QQ and Coefficient plots for significant variables have been saved.")



###################################################################################
#################### LOGISTIC REGRESSION WITH PRS, PC1, MEAN ######################
###################################################################################

# # Load required libraries
# #library(tidyverse)
# #library(broom)
# #library(caret)
# #library(pROC)
# #library(ggplot2)
# 
# Define PRS predictors in CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_vars_c <- grep("prs_c_qc_", names(data), value = TRUE)
prs_vars_m <- grep("prs_m_qc_", names(data), value = TRUE)
prs_vars_f <- grep("prs_f_qc_", names(data), value = TRUE)
rs_vars <- c(prs_vars_c, prs_vars_m, prs_vars_f)

# Define covariates in CBCL1_5 (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Function to perform logistic regression
logistic_regression_singular <- function(data, predictor) {
  formula <- as.formula(paste("h_c_total_cat_1_5y ~", paste(c(predictor, covariates), collapse = "+")))

  # Fit the logistic regression model
  model <- glm(formula, data = data, family = binomial)

  # Model evaluation: Wald test Coefficients, ROC_AUC, Confusion Matrix, and Likelihood Ratio Test for singular model
  # Summarize Wald test results
  summary <- tidy(model)
  summary <- summary %>%
    mutate(
      OR = exp(estimate),
      CI_lower = exp(estimate - 1.96 * std.error),
      CI_upper = exp(estimate + 1.96 * std.error)
    )
  
  # Create ROC curve and calculate AUC
  predicted_probs <- predict(model, type = "response")
  actual_classes <- factor(data$h_c_total_cat_1_5y, levels = c(0, 1))
  roc_obj <- roc(actual_classes, predicted_probs)
  auc <- auc(roc_obj)
  summary$AUC <- auc
  
  # Create confusion matrix and perform Fisher test
  predicted_classes <- factor(ifelse(predicted_probs > 0.5, 1, 0), levels = c(0, 1))
  confusion_matrix <- confusionMatrix(predicted_classes, actual_classes)
  confusion_table <- confusion_matrix$table
  fisher_test <- fisher.test(confusion_table)
  summary$ConfM_Fisher_p_value <- fisher_test$p.value

  # Perform Likelihood Ratio Test and calculate overall model fitness
  lrt_p_value <- pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual, lower.tail = FALSE)
  summary$LRT_p_value <- lrt_p_value

  # Calculate and save AIC and BIC
  summary$AIC <- AIC(model)
  summary$BIC <- BIC(model)

  # Extract the number of observations (N) from the fitted glm model
  N_model <- model$df.residual + model$rank  # df.residual + rank gives the total N

  # Check if any result is significant
  significant <- any(summary$p.value[2] < 0.05)
  
  # Return the model and summary as part of a list
  return(list(model = model, summary = summary[2,], confusion_matrix = confusion_table, roc_obj = roc_obj, N_model = N_model, significant = significant))
}

# Ensure directories exist
if (!dir.exists("LogR_Single_PRS_CM_Plots_CBCL1_5")) {
  dir.create("LogR_Single_PRS_CM_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LogR_Single_PRS_ROC_Plots_CBCL1_5")) {
  dir.create("LogR_Single_PRS_ROC_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LogR_Single_PRS_Coeff_Plots_CBCL1_5")) {
  dir.create("LogR_Single_PRS_Coeff_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Print plots for significant results 
print("Plots of LogR Single PRS Significant Results for CBCL1_5:")

# Function to process and plot significant logistic regression results
process_and_plot_results <- function(predictor_name, results) {
  # Plot and save confusion matrix
  cm_data <- as.data.frame(results$confusion_matrix)
  png(paste0("LogR_Single_PRS_CM_Plots_CBCL1_5/cm_single_", predictor_name, "_h_c_total_cat_cbcl1_5y", ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(as.factor(cm_data$Prediction), as.factor(cm_data$Reference), col = cm_data$Freq, main = paste("Confusion Matrix for", predictor_name))
  dev.off()
  plot(as.factor(cm_data$Prediction), as.factor(cm_data$Reference), col = cm_data$Freq, main = paste("Confusion Matrix for", predictor_name))
  
  # Plot and save ROC Curve
  png(paste0("LogR_Single_PRS_ROC_Plots_CBCL1_5/roc_single_", predictor_name, "_h_c_total_cat_cbcl1_5y", ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(results$roc_obj, main = paste("ROC Curve for", predictor_name))
  dev.off()
  plot(results$roc_obj, main = paste("ROC Curve for", predictor_name))
  
  # Plot and save model coefficients
  coef_matrix <- summary(results$model)$coefficients
  coef_df <- as.data.frame(t(coef_matrix[2, ]))  # Transpose to ensure a single row data frame
  colnames(coef_df) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  coef_df$lower_ci <- coef_df$Estimate - 1.96 * coef_df$`Std. Error`
  coef_df$upper_ci <- coef_df$Estimate + 1.96 * coef_df$`Std. Error`
  coef_df$Predictor <- rownames(coef_matrix)[2]
  
  coef_plot <- ggplot(coef_df, aes(x = Predictor, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    theme_minimal() +
    coord_flip() +  # Flip the coordinates for better readability
    ggtitle(paste("Coefficient Plot of LogR for Single", predictor_name, "on h_c_total_cat_cbcl1_5y")) +
    ylab("Estimate") +
    xlab("Predictor")
  
  # Save the coefficient plot
  png(paste0("LogR_Single_PRS_Coeff_Plots_CBCL1_5/logr_coeff_single_", predictor_name, "_h_c_total_cat_cbcl1_5y", ".png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  print(coef_plot)
  dev.off()
  print(coef_plot)
}

# Main loop to run logistic regression analysis and plot significant results
results_logr_single <- lapply(prs_vars, function(predictor_name) {
  results <- logistic_regression_singular(data, predictor_name)
  
  if (results$significant) {
    process_and_plot_results(predictor_name, results)
  }
  
  return(data.frame(
    CBCL_outcome = "h_c_total_cat_1_5y",
    PRS_predictor = predictor_name,
    Beta = results$summary$estimate,
    SE = results$summary$std.error,
    Wald_Z = results$summary$statistic,
    Wald_p_value = results$summary$p.value,
    Odds_ratio = results$summary$OR,
    CI_lower = results$summary$CI_lower,
    CI_upper = results$summary$CI_upper,
    AUC = results$summary$AUC,
    ConfM_Fisher_p_value = results$summary$ConfM_Fisher_p_value,
    LRT_chisqr_p_value = results$summary$LRT_p_value,
    AIC = results$summary$AIC,
    BIC = results$summary$BIC,
    N = results$N_model
  ))
})

# Combine all logistic regression results into a single data frame
results_logr_single_df <- do.call(rbind, results_logr_single)

# View all logistic regression results
head(results_logr_single_df)

# Save all logistic regression results to a CSV file
write.csv(results_logr_single_df, "logr_prs_summary_cbcl1_5.csv", row.names=FALSE)

# Extract and print significant logistic regression results
significant_results_logr_single_df <- subset(results_logr_single_df, Wald_p_value < 0.05)
print(significant_results_logr_single_df)

# Save significant logistic regression results to a csv file
write.csv(significant_results_logr_single_df, "logr_prs_significant_cbcl1_5.csv", row.names = FALSE)

# Completion message
print("LogR csv file and plots for significant variables have been saved.")





###################################################################################
#################### TRIO-PRS, -PC1, -MEAN LOGISTIC REGRESSION ####################
###################################################################################

# # Load required libraries
# #library(tidyverse)
# #library(broom)
# #library(caret)
# #library(pROC)
# #library(ggplot2)
# 
# Define PRS variables for child, mother, father in CBCL1_5 (adjust these to your dataset's variable names if necessary)
prs_vars_c <- grep("prs_c_qc_", names(data), value = TRUE)
prs_vars_m <- grep("prs_m_qc_", names(data), value = TRUE)
prs_vars_f <- grep("prs_f_qc_", names(data), value = TRUE)

# Define CBCL1_5 outcome variables (adjust these to your dataset's variable names if necessary)
covariates <- c("h_c_age_1_5y", "h_c_sex", "h_m_age", "h_m_edu_1_5y", paste0("anc_pc_", 1:10))

# Function to perform logistic regression
logistic_regression_trio <- function(data, prs_child_var, prs_mother_var, prs_father_var, covariates) {
  # Constructing the model formula
  formula <- as.formula(paste("h_c_total_cat_1_5y", "~", paste(c(prs_child_var, prs_mother_var, prs_father_var, covariates), collapse = "+")))

  # Fit the logistic regression model
  model <- glm(formula, data = data, family = binomial)

  # Model evaluation: Wald test Coefficients, ROC_AUC, Confusion Matrix, and Likelihood Ratio Test for singular model
  # Summarize Wald test results
  summary <- tidy(model)
  summary <- summary %>%
    mutate(
      OR = exp(estimate),
      CI_lower = exp(estimate - 1.96 * std.error),
      CI_upper = exp(estimate + 1.96 * std.error)
    )

  # Create ROC curve and calculate AUC
  predicted_probs <- predict(model, type = "response")
  actual_classes <- factor(data$h_c_total_cat_1_5y, levels = c(0, 1))
  roc_obj <- roc(actual_classes, predicted_probs)
  auc <- auc(roc_obj)
  summary$AUC <- auc

  # Create confusion matrix and perform Fisher test
  predicted_classes <- factor(ifelse(predicted_probs > 0.5, 1, 0), levels = c(0, 1))
  confusion_matrix <- confusionMatrix(predicted_classes, actual_classes)
  confusion_table <- confusion_matrix$table
  fisher_test <- fisher.test(confusion_table)
  summary$ConfM_Fisher_p_value <- fisher_test$p.value

  # Perform Likelihood Ratio Test and calculate overall model fitness
  lrt_p_value <- pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual, lower.tail = FALSE)
  summary$LRT_p_value <- lrt_p_value

  # Calculate and save AIC and BIC
  summary$AIC <- AIC(model)
  summary$BIC <- BIC(model)

  # Extract the number of observations (N) from the fitted glm model
  N_model <- model$df.residual + model$rank  # df.residual + rank gives the total N
  
  # Create a list of statistic results
  return(list(model = model, summary = summary[2:4,], confusion_matrix = confusion_table, roc_obj = roc_obj, N_model = N_model))
}

# Initialize results list
results_logr_trio <- list()

# Perform trio logistic regression analysis for CBCL1_5
for (i in 1:length(prs_vars_c)) {
  results <- logistic_regression_trio(data, prs_vars_c[i], prs_vars_m[i], prs_vars_f[i], covariates)
  results_logr_trio <- append(results_logr_trio, list(list(prs_child_var = prs_vars_c[i], prs_mother_var = prs_vars_m[i], prs_father_var = prs_vars_f[i], results = results)))
}

# Compile all logistic regression results into a data frame
results_logr_trio_df <- do.call(rbind, lapply(results_logr_trio, function(res) {
  data.frame(
    CBCL_outcome = "h_c_total_cat_1_5y",
    PRS_predictor = res$results$summary$term,
    Beta = res$results$summary$estimate,
    SE = res$results$summary$std.error,
    Wald_Z = res$results$summary$statistic,
    Wald_p_value = res$results$summary$p.value,
    Odds_ratio = res$results$summary$OR,
    CI_lower = res$results$summary$CI_lower,
    CI_upper = res$results$summary$CI_upper,
    AUC = res$results$summary$AUC,
    ConfM_Fisher_p_value = res$results$summary$ConfM_Fisher_p_value,
    LRT_chisqr_p_value = res$results$summary$LRT_p_value,
    AIC = res$results$summary$AIC,
    BIC = res$results$summary$BIC,
    N = results$N_model
  )
}))

# View all logistic regression results
head(results_logr_trio_df)

# Save all logistic regression results to a CSV file
write.csv(results_logr_trio_df, "logr_trio_PRScmf_summary_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Extract and print significant logistic regression results
significant_results_logr_trio <- subset(results_logr_trio_df, results_logr_trio_df[,6] < 0.05)
print(significant_results_logr_trio)

# Save significant logistic regression results to a csv file
write.csv(significant_results_logr_trio, "logr_trio_PRScmf_significant_cbcl1_5.csv", row.names = FALSE) # replace cbcl1_5 to specify your own cohort

# Ensure directories exist
if (!dir.exists("LogR_Trio_PRScmf_CM_Plots_CBCL1_5")) {
  dir.create("LogR_Trio_PRScmf_CM_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LogR_Trio_PRScmf_ROC_Plots_CBCL1_5")) {
  dir.create("LogR_Trio_PRScmf_ROC_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}
if (!dir.exists("LogR_Trio_PRScmf_Coeff_Plots_CBCL1_5")) {
  dir.create("LogR_Trio_PRScmf_Coeff_Plots_CBCL1_5") # replace CBCL1_5 in the directory name to specify your own cohort
}

# Print plots for significant results
print("Plots of LogR Trio PRS Significant Results for CBCL1_5:")

# Function to plot and save confusion matrix
plot_and_save_confusion_matrix <- function(cm, predictor_name) {
  cm_data <- as.data.frame(cm)
  png(paste0("LogR_Trio_PRScmf_CM_Plots_CBCL1_5/cm_PRScmf_", predictor_name, "_h_c_total_cat_cbcl1_5y.png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(as.factor(cm_data$Prediction), as.factor(cm_data$Reference), col = cm_data$Freq, main = paste("Confusion Matrix for", predictor_name))
  dev.off()
  plot(as.factor(cm_data$Prediction), as.factor(cm_data$Reference), col = cm_data$Freq, main = paste("Confusion Matrix for", predictor_name))
}

# Function to plot and save ROC curve
plot_and_save_roc_curve <- function(roc_obj, predictor_name) {
  png(paste0("LogR_Trio_PRScmf_ROC_Plots_CBCL1_5/roc_PRScmf_", predictor_name, "_h_c_total_cat_cbcl1_5y.png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  plot(roc_obj, main = paste("ROC Curve for", predictor_name))
  dev.off()
  plot(roc_obj, main = paste("ROC Curve for", predictor_name))
}

# Function to plot and save model coefficients
plot_and_save_model_coefficients <- function(model, predictor_name) {
  # Extract the coefficients and standard errors (only rows 2 to 4)
  coef_matrix <- summary(model)$coefficients
  coef_matrix <- coef_matrix[2:4, ]
  
  # Convert the selected rows to a data frame for plotting
  coef_df <- as.data.frame(coef_matrix)
  coef_df$Predictor <- rownames(coef_df)
  coef_df$lower_ci <- coef_df$Estimate - 1.96 * coef_df$`Std. Error`
  coef_df$upper_ci <- coef_df$Estimate + 1.96 * coef_df$`Std. Error`
  
  # Create the plot
  coef_plot <- ggplot(coef_df, aes(x = Predictor, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    theme_minimal() +
    coord_flip() +
    ggtitle(paste("Coefficient Plot of LogR for PRScmf:", predictor_name, "on h_c_total_cat_cbcl1_5y")) +
    ylab("Estimate") +
    xlab("Predictor")
  
  # Save the coefficient plot
  png(paste0("LogR_Trio_PRScmf_Coeff_Plots_CBCL1_5/logr_coeff_PRScmf_", predictor_name, "_h_c_total_cat_cbcl1_5y.png"), width = 700, height = 500) # replace cbcl1_5 to specify your own cohort
  print(coef_plot)
  dev.off()
  print(coef_plot)
}

# Loop through results and plot only significant ones
for (i in seq_along(results_logr_trio)) {
  # Extract the current result item
  current_result <- results_logr_trio[[i]]
  
  # Extract the p-values for child, mother, and father PRS
  p_values <- current_result$results$summary$p.value
  
  # List of PRS variables
  prs_vars <- c(current_result$prs_child_var, current_result$prs_mother_var, current_result$prs_father_var)
  
  # Loop through each PRS variable (child, mother, father) and check significance
  for (j in seq_along(prs_vars)) {
    prs_var <- prs_vars[j]
    
    # Check if this PRS variable is significant
    if (p_values[j] < 0.05) {
      # Perform plotting function using only rows 2 to 4 of the coefficients matrix
      plot_and_save_confusion_matrix(current_result$results$confusion_matrix, prs_var)
      plot_and_save_roc_curve(current_result$results$roc_obj, prs_var)
      plot_and_save_model_coefficients(current_result$results$model, prs_var)
    }
  }
}

# Completion message
print("Trio_LogR csv file and plots for significant variables have been saved.")



