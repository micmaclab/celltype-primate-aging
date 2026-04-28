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


result_list <- list(
  region       = character(),
  age_t        = numeric(),
  age_p        = numeric(),
  age_estimate = numeric(),
  age_ci_low   = numeric(),
  age_ci_high  = numeric(),
  sex_t        = numeric(),
  sex_p        = numeric(),
  sex_estimate = numeric(),
  sex_ci_low   = numeric(),
  sex_ci_high  = numeric(),
  intercept    = numeric()
)

for (region_id in unique(data$region)) {
  
  region_data <- filter(data, region == region_id)
  model <- lmer(value ~ age + sex + (1|hemi), data = region_data, REML = TRUE)
  
  coefs <- as.data.frame(lmerTest:::get_coefmat(model))
  ci    <- confint(model, parm = "beta_", method = "profile")
  
  result_list$region       <- c(result_list$region,       region_id)
  result_list$intercept    <- c(result_list$intercept,    coefs$Estimate[1])
  result_list$age_estimate <- c(result_list$age_estimate, coefs["age",  "Estimate"])
  result_list$age_t        <- c(result_list$age_t,        coefs["age",  "t value"])
  result_list$age_p        <- c(result_list$age_p,        coefs["age",  "Pr(>|t|)"])
  result_list$age_ci_low   <- c(result_list$age_ci_low,   ci["age",  1])
  result_list$age_ci_high  <- c(result_list$age_ci_high,  ci["age",  2])
  result_list$sex_estimate <- c(result_list$sex_estimate, coefs["sexM", "Estimate"])
  result_list$sex_t        <- c(result_list$sex_t,        coefs["sexM", "t value"])
  result_list$sex_p        <- c(result_list$sex_p,        coefs["sexM", "Pr(>|t|)"])
  result_list$sex_ci_low   <- c(result_list$sex_ci_low,   ci["sexM", 1])
  result_list$sex_ci_high  <- c(result_list$sex_ci_high,  ci["sexM", 2])
}

df_model <- as.data.frame(result_list)
write.csv(df_model,
          file = file.path(project_root, "regionwise_age_effects_MixedLM.csv"),
          row.names = FALSE)
print(df_model)
