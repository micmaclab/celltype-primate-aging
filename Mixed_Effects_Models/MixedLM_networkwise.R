# Load required packages
library(lme4)
library(lmerTest)
library(merTools)
library(ggplot2)
library(dplyr)
library(here)

# Set working directory to the project root (one level above this script)
project_root <- dirname(here::here())
setwd(project_root)
cat("Working directory set to:", getwd(), "\n")

# ------------------------
# Load and Prepare Data
# ------------------------

# Load regional totaly similarity strength data
data_path <- file.path(project_root, "MIND_Network", "similarity_strength_subject_data.csv")
data <- read.csv(data_path, header = TRUE)

# Initialize a list to store model results
result_list <- list(
  yeo_label = character(),
  p_value = numeric(),
  t_value = numeric(),
  intercept = numeric(),
  slope = numeric()
)

# ------------------------
# Fit Models by Network
# ------------------------

for (network in unique(data$yeo_label)) {
  
  # Subset data for the current network
  network_data <- filter(data, yeo_label == network)
  
  # Fit mixed effects model
  model <- lmer(value ~ age + (1 | hemi) + (1 | subject) + (1 | sex),
                data = network_data, REML = TRUE)
  
  # Plot smoothed relationship between age and Total Similarity Strength
  avg_by_subject <- network_data %>%
    group_by(subject, age) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  gg <- ggplot(avg_by_subject, aes(x = age, y = value)) +
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray70") +
    labs(title = paste0(network, ": Total Similarity Strength vs Age"),
         x = "Age", y = "Total Similarity Strength") +
    theme_minimal(base_size = 14)
  
  # (Optional) Save smoothed data — uncomment to use
  # smoothed_data <- ggplot_build(gg)$data[[2]][, c("x", "y", "ymin", "ymax", "se")]
  # write.csv(smoothed_data, file = file.path(".../ggplot_output", paste0("predictions_", gsub(" ", "_", network), ".csv")),
  #           row.names = FALSE)
  
  # Extract model results
  res_df <- as.data.frame(lmerTest:::get_coefmat(model))
  
  # Append to results list
  result_list$yeo_label <- c(result_list$yeo_label, network)
  result_list$p_value <- c(result_list$p_value, res_df$`Pr(>|t|)`[2])
  result_list$t_value <- c(result_list$t_value, res_df$`t value`[2])
  result_list$intercept <- c(result_list$intercept, res_df$Estimate[1])
  result_list$slope <- c(result_list$slope, res_df$Estimate[2])
}

# ------------------------
# Fit Whole-Brain Model
# ------------------------

whole_model <- lmer(value ~ age + (1 | hemi) + (1 | subject) + (1 | sex),
                    data = data, REML = TRUE)

avg_by_subject <- data %>%
  group_by(subject, age) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

gg <- ggplot(avg_by_subject, aes(x = age, y = value)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray70") +
  labs(title = "Whole Brain: Total Similarity Strength vs Age",
       x = "Age", y = "Total Similarity Strength") +
  theme_minimal(base_size = 14)

# (Optional) Save smoothed data — uncomment to use
# smoothed_data <- ggplot_build(gg)$data[[2]][, c("x", "y", "ymin", "ymax", "se")]
# write.csv(smoothed_data, file = file.path(".../ggplot_output", "predictions_Whole_Brain.csv"),
#           row.names = FALSE)

res_df <- as.data.frame(lmerTest:::get_coefmat(whole_model))

# Append whole-brain results to result list
result_list$yeo_label <- c(result_list$yeo_label, "Whole Brain")
result_list$p_value <- c(result_list$p_value, res_df$`Pr(>|t|)`[2])
result_list$t_value <- c(result_list$t_value, res_df$`t value`[2])
result_list$intercept <- c(result_list$intercept, res_df$Estimate[1])
result_list$slope <- c(result_list$slope, res_df$Estimate[2])

# ------------------------
# Save Results
# ------------------------

final_results <- as.data.frame(result_list)

# Save to CSV (change path as needed)
write.csv(final_results,
          file = file.path(project_root, "Mixed_Effects_Models", "networkwise_age_effects_MixedLM.csv"),
          row.names = FALSE)

# Optional: print results
print(final_results)
