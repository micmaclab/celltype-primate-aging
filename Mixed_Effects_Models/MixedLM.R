# Load required libraries
library(lme4)
library(lmerTest)
library(here)
library(dplyr)

# Set working directory to the project root (one level above current folder)
project_root <- dirname(here::here())
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# ------------------------
# Load and Prepare Data
# ------------------------

# Read in regional total similarity strength data
data_path <- file.path(project_root, "MIND_Network", "similarity_strength_subject_data.csv")
data <- read.csv(data_path, header = TRUE)

# Initialize list to store results
result_list <- list(
  region = character(),
  p_value = numeric(),
  t_value = numeric(),
  intercept = numeric(),
  slope = numeric()
)

# ------------------------
# Fit Region-Wise Models
# ------------------------

for (region_id in unique(data$region)) {
  
  # Filter data for current region
  region_data <- filter(data, region == region_id)
  
  # Fit linear mixed effects model
  model <- lmer(value ~ age + (1 | hemi) + (1 | sex), 
                data = region_data, REML = TRUE)
  
  # Extract model coefficients
  res_df <- as.data.frame(lmerTest:::get_coefmat(model))
  
  # Append results to list
  result_list$region <- c(result_list$region, region_id)
  result_list$p_value <- c(result_list$p_value, res_df$`Pr(>|t|)`[2])
  result_list$t_value <- c(result_list$t_value, res_df$`t value`[2])
  result_list$intercept <- c(result_list$intercept, res_df$Estimate[1])
  result_list$slope <- c(result_list$slope, res_df$Estimate[2])
}

# ------------------------
# Save Results
# ------------------------

# Convert result list to data frame
df <- as.data.frame(result_list)

# Export results to CSV
output_path <- file.path(project_root, "Mixed_Effects_Models", "regionwise_age_effects_MixedLM.csv")
write.csv(df, file = output_path, row.names = FALSE)


print(df)
