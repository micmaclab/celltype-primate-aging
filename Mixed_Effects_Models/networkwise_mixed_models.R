
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(dplyr)
library(here)

project_root <- here()
data_path    <- file.path(project_root,"MIND_Network", "similarity_strength_subject_data.csv")
output_dir   <- file.path(project_root, "ggplot_output")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

data <- read.csv(data_path, header = TRUE)

# Demean age so intercepts are interpretable at the sample mean
mean_age          <- mean(data$age, na.rm = TRUE)
data$age_demeaned <- data$age - mean_age

# -----------------------------------------------------------------------------
# Network-wise Models
# -----------------------------------------------------------------------------

networks <- unique(data$yeo_label)

results <- data.frame(
  yeo_label          = character(),
  p_value            = numeric(),
  t_value            = numeric(),
  intercept          = numeric(),
  slope              = numeric(),
  intercept_demeaned = numeric(),
  slope_demeaned     = numeric(),
  stringsAsFactors   = FALSE
)

for (network in networks) {
  
  message("Fitting model for network: ", network)
  
  network_data <- filter(data, yeo_label == network)
  
  # Per-subject averages for plotting observed data
  avg_by_subject <- network_data %>%
    group_by(subject, age) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  # --- Model: raw age ---
  model  <- lmer(value ~ age + sex + (1 | hemi) + (1 | subject),
                 data = network_data, REML = TRUE)
  res_df <- as.data.frame(lmerTest:::get_coefmat(model))
  
  # Marginal predictions from model
  pred_df <- ggpredict(model, terms = "age")
  write.csv(as.data.frame(pred_df),
            file = file.path(output_dir, paste0("model_predictions_", gsub(" ", "_", network), ".csv")),
            row.names = FALSE)
  
  # Plot observed averages + model predictions
  gg_pred <- ggplot() +
    geom_point(data = avg_by_subject, aes(x = age, y = value),
               color = "black", alpha = 0.5) +
    geom_line(data = pred_df, aes(x = x, y = predicted),
              color = "blue", linewidth = 1) +
    geom_ribbon(data = pred_df, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "blue", alpha = 0.2) +
    labs(title = paste0(network, ": Predicted Strength vs Age"),
         x = "Age", y = "Total Similarity Strength") +
    theme_minimal(base_size = 14)
  print(gg_pred)
  
  # ggplot smooth and save underlying smooth data
  gg_smooth <- ggplot(avg_by_subject, aes(x = age, y = value)) +
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray70") +
    labs(title = paste0(network, ": Total Similarity Strength vs Age"),
         x = "Age", y = "Total Similarity Strength") +
    theme_minimal(base_size = 14)
  
  smooth_data <- ggplot_build(gg_smooth)$data[[2]][, c("x", "y", "ymin", "ymax", "se")]
  write.csv(smooth_data,
            file = file.path(output_dir, paste0("predictions_", gsub(" ", "_", network), "2.csv")),
            row.names = FALSE)
  
  # --- Model: demeaned age (intercept = value at mean age) ---
  model_demeaned <- lmer(value ~ age_demeaned + sex + (1 | hemi) + (1 | subject),
                         data = network_data, REML = TRUE)
  res_demeaned   <- as.data.frame(summary(model_demeaned)$coefficients)
  
  # Collect results
  results <- rbind(results, data.frame(
    yeo_label          = network,
    p_value            = res_df$`Pr(>|t|)`[2],
    t_value            = res_df$`t value`[2],
    intercept          = res_df$Estimate[1],
    slope              = res_df$Estimate[2],
    intercept_demeaned = res_demeaned$Estimate[1],
    slope_demeaned     = res_demeaned$Estimate[2],
    stringsAsFactors   = FALSE
  ))
}

# -----------------------------------------------------------------------------
# Whole-Brain Model
# -----------------------------------------------------------------------------

message("Fitting whole-brain model...")

avg_whole <- data %>%
  group_by(subject, age) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

whole_model <- lmer(value ~ age + sex + (1 | hemi) + (1 | subject),
                    data = data, REML = TRUE)
res_whole   <- as.data.frame(lmerTest:::get_coefmat(whole_model))

pred_whole <- ggpredict(whole_model, terms = "age")
write.csv(as.data.frame(pred_whole),
          file = file.path(output_dir, "model_predictions_Whole_Brain2.csv"),
          row.names = FALSE)

gg_whole_pred <- ggplot() +
  geom_point(data = avg_whole, aes(x = age, y = value),
             color = "black", alpha = 0.5) +
  geom_line(data = pred_whole, aes(x = x, y = predicted),
            color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_whole, aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "blue", alpha = 0.2) +
  labs(title = "Whole Brain: Predicted Strength vs Age",
       x = "Age", y = "Total Similarity Strength") +
  theme_minimal(base_size = 14)
print(gg_whole_pred)

gg_whole_smooth <- ggplot(avg_whole, aes(x = age, y = value)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray70") +
  labs(title = "Whole Brain: Total Similarity Strength vs Age",
       x = "Age", y = "Total Similarity Strength") +
  theme_minimal(base_size = 14)

smooth_whole <- ggplot_build(gg_whole_smooth)$data[[2]][, c("x", "y", "ymin", "ymax", "se")]
write.csv(smooth_whole,
          file = file.path(output_dir, "predictions_Whole_Brain.csv"),
          row.names = FALSE)

whole_model_demeaned <- lmer(value ~ age_demeaned + sex + (1 | hemi) + (1 | subject),
                             data = data, REML = TRUE)
res_whole_demeaned   <- as.data.frame(summary(whole_model_demeaned)$coefficients)

results <- rbind(results, data.frame(
  yeo_label          = "Whole Brain",
  p_value            = res_whole$`Pr(>|t|)`[2],
  t_value            = res_whole$`t value`[2],
  intercept          = res_whole$Estimate[1],
  slope              = res_whole$Estimate[2],
  intercept_demeaned = res_whole_demeaned$Estimate[1],
  slope_demeaned     = res_whole_demeaned$Estimate[2],
  stringsAsFactors   = FALSE
))

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

write.csv(results,
          file = file.path(project_root, "celltype-primate-aging",
                           "networkwise_age_effects_MixedLM.csv"),
          row.names = FALSE)

print(results)
message("Done. Results saved to: ", file.path(project_root, "celltype-primate-aging"))