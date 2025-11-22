# Load necessary library
library(dplyr)

# Read the CSV file
data <- read.csv("lrm_trio_PRScmf_summary_sdq6_18_mo_ALSPAC.csv")

# Remove rows where PRS_predictor starts with "prs_f_"
filtered_data <- data %>% filter(!grepl("^prs_f_", PRS_predictor))

# Save the cleaned data to a new CSV file
write.csv(filtered_data, "lrm_trio_PRScmf_-PRSf_summary_sdq6_18_mo_ALSPAC.csv", row.names = FALSE)
