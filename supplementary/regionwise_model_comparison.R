# Load required libraries
library(lme4)
library(lmerTest)
library(dplyr)

project_root <- here()

data_path <- file.path(project_root, "similarity_strength_subject_data.csv")
data <- read.csv(data_path, header = TRUE)

#---------------------------------MODEL1----------------------------------------
# total similarity strength ~ age + (1|sex) + (1|hemisphere)

result_list <- list(
  region       = character(),
  age_t        = numeric(),
  age_p        = numeric(),
  age_estimate = numeric(),
  age_ci_low   = numeric(),
  age_ci_high  = numeric(),
  intercept    = numeric()
)

for (region_id in unique(data$region)) {
  
  region_data <- filter(data, region == region_id)
  model <- lmer(value ~ age + (1|sex) + (1|hemi), data = region_data, REML = TRUE)
  
  coefs <- as.data.frame(lmerTest:::get_coefmat(model))
  ci    <- confint(model, parm = "beta_", method = "profile")
  
  result_list$region       <- c(result_list$region,       region_id)
  result_list$intercept    <- c(result_list$intercept,    coefs$Estimate[1])
  result_list$age_estimate <- c(result_list$age_estimate, coefs["age", "Estimate"])
  result_list$age_t        <- c(result_list$age_t,        coefs["age", "t value"])
  result_list$age_p        <- c(result_list$age_p,        coefs["age", "Pr(>|t|)"])
  result_list$age_ci_low   <- c(result_list$age_ci_low,   ci["age", 1])
  result_list$age_ci_high  <- c(result_list$age_ci_high,  ci["age", 2])
}

df_model1 <- as.data.frame(result_list)
write.csv(df_model1,
          file = file.path(project_root, "regional_model_testing/lmer_MIND_fixed_age_random_sex_hemi.csv"),
          row.names = FALSE)
print(df_model1)

#---------------------------------MODEL2----------------------------------------
# total similarity strength ~ age + sex + (1|hemisphere)

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

df_model2 <- as.data.frame(result_list)
write.csv(df_model2,
          file = file.path(project_root, "regional_model_testing/lmer_MIND_fixed_age_sex_random_hemi.csv"),
          row.names = FALSE)
print(df_model2)

#---------------------------------MODEL3----------------------------------------
# total similarity strength ~ age + sex + age*sex + (1|hemisphere)

result_list <- list(
  region           = character(),
  age_t            = numeric(),
  age_p            = numeric(),
  age_estimate     = numeric(),
  age_ci_low       = numeric(),
  age_ci_high      = numeric(),
  sex_t            = numeric(),
  sex_p            = numeric(),
  sex_estimate     = numeric(),
  sex_ci_low       = numeric(),
  sex_ci_high      = numeric(),
  age_sex_t        = numeric(),
  age_sex_p        = numeric(),
  age_sex_estimate = numeric(),
  age_sex_ci_low   = numeric(),
  age_sex_ci_high  = numeric(),
  intercept        = numeric()
)

for (region_id in unique(data$region)) {
  
  region_data <- filter(data, region == region_id)
  model <- lmer(value ~ age + sex + age:sex + (1|hemi), data = region_data, REML = TRUE)
  
  coefs <- as.data.frame(lmerTest:::get_coefmat(model))
  ci    <- confint(model, parm = "beta_", method = "profile")
  
  result_list$region           <- c(result_list$region,           region_id)
  result_list$intercept        <- c(result_list$intercept,        coefs$Estimate[1])
  result_list$age_estimate     <- c(result_list$age_estimate,     coefs["age",      "Estimate"])
  result_list$age_t            <- c(result_list$age_t,            coefs["age",      "t value"])
  result_list$age_p            <- c(result_list$age_p,            coefs["age",      "Pr(>|t|)"])
  result_list$age_ci_low       <- c(result_list$age_ci_low,       ci["age",      1])
  result_list$age_ci_high      <- c(result_list$age_ci_high,      ci["age",      2])
  result_list$sex_estimate     <- c(result_list$sex_estimate,     coefs["sexM",     "Estimate"])
  result_list$sex_t            <- c(result_list$sex_t,            coefs["sexM",     "t value"])
  result_list$sex_p            <- c(result_list$sex_p,            coefs["sexM",     "Pr(>|t|)"])
  result_list$sex_ci_low       <- c(result_list$sex_ci_low,       ci["sexM",     1])
  result_list$sex_ci_high      <- c(result_list$sex_ci_high,      ci["sexM",     2])
  result_list$age_sex_estimate <- c(result_list$age_sex_estimate, coefs["age:sexM", "Estimate"])
  result_list$age_sex_t        <- c(result_list$age_sex_t,        coefs["age:sexM", "t value"])
  result_list$age_sex_p        <- c(result_list$age_sex_p,        coefs["age:sexM", "Pr(>|t|)"])
  result_list$age_sex_ci_low   <- c(result_list$age_sex_ci_low,   ci["age:sexM", 1])
  result_list$age_sex_ci_high  <- c(result_list$age_sex_ci_high,  ci["age:sexM", 2])
}

df_model3 <- as.data.frame(result_list)
write.csv(df_model3,
          file = file.path(project_root, "regional_model_testing/lmer_MIND_fixed_age_sex_ageBYsex_random_hemi.csv"),
          row.names = FALSE)
print(df_model3)

#---------------------------------MODEL4----------------------------------------
# total similarity strength ~ age + sex + hemisphere

result_list <- list(
  region        = character(),
  age_t         = numeric(),
  age_p         = numeric(),
  age_estimate  = numeric(),
  age_ci_low    = numeric(),
  age_ci_high   = numeric(),
  sex_t         = numeric(),
  sex_p         = numeric(),
  sex_estimate  = numeric(),
  sex_ci_low    = numeric(),
  sex_ci_high   = numeric(),
  hemi_t        = numeric(),
  hemi_p        = numeric(),
  hemi_estimate = numeric(),
  hemi_ci_low   = numeric(),
  hemi_ci_high  = numeric(),
  intercept     = numeric()
)

for (region_id in unique(data$region)) {
  
  region_data <- filter(data, region == region_id)
  model <- lm(value ~ age + sex + hemi, data = region_data)
  
  coefs <- as.data.frame(summary(model)$coefficients)
  ci    <- confint(model)
  
  result_list$region        <- c(result_list$region,        region_id)
  result_list$intercept     <- c(result_list$intercept,     coefs$Estimate[1])
  result_list$age_estimate  <- c(result_list$age_estimate,  coefs["age",       "Estimate"])
  result_list$age_t         <- c(result_list$age_t,         coefs["age",       "t value"])
  result_list$age_p         <- c(result_list$age_p,         coefs["age",       "Pr(>|t|)"])
  result_list$age_ci_low    <- c(result_list$age_ci_low,    ci["age",       1])
  result_list$age_ci_high   <- c(result_list$age_ci_high,   ci["age",       2])
  result_list$sex_estimate  <- c(result_list$sex_estimate,  coefs["sexM",      "Estimate"])
  result_list$sex_t         <- c(result_list$sex_t,         coefs["sexM",      "t value"])
  result_list$sex_p         <- c(result_list$sex_p,         coefs["sexM",      "Pr(>|t|)"])
  result_list$sex_ci_low    <- c(result_list$sex_ci_low,    ci["sexM",      1])
  result_list$sex_ci_high   <- c(result_list$sex_ci_high,   ci["sexM",      2])
  result_list$hemi_estimate <- c(result_list$hemi_estimate, coefs["hemiright", "Estimate"])
  result_list$hemi_t        <- c(result_list$hemi_t,        coefs["hemiright", "t value"])
  result_list$hemi_p        <- c(result_list$hemi_p,        coefs["hemiright", "Pr(>|t|)"])
  result_list$hemi_ci_low   <- c(result_list$hemi_ci_low,   ci["hemiright", 1])
  result_list$hemi_ci_high  <- c(result_list$hemi_ci_high,  ci["hemiright", 2])
}

df_model4 <- as.data.frame(result_list)
write.csv(df_model4,
          file = file.path(project_root, "regional_model_testing/lm_MIND_fixed_age_sex_hemi.csv"),
          row.names = FALSE)
print(df_model4)
