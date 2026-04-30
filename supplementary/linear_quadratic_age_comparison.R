library(lme4)
library(lmerTest)
library(here)
library(dplyr)


data_path <- here("similarity_strength_subject_data.csv")

data <- read.csv(data_path, header = TRUE)

# ------------------------------------------------------------------------------
# Initialize results 
# ------------------------------------------------------------------------------
result_list <- list(
  region = character(),
  
  # Model 1: Linear
  lme_intercept      = numeric(),
  lme_age_estimate   = numeric(),
  lme_age_t          = numeric(),
  lme_age_p          = numeric(),
  aic_lme_linear     = numeric(),
  bic_lme_linear     = numeric(),
  
  # Model 2: Quadratic
  lme2_intercept      = numeric(),
  lme2_age_estimate   = numeric(),
  lme2_age_t          = numeric(),
  lme2_age_p          = numeric(),
  lme2_age2_estimate  = numeric(),
  lme2_age2_t         = numeric(),
  lme2_age2_p         = numeric(),
  aic_lme_quadratic   = numeric(),
  bic_lme_quadratic   = numeric(),
  
  # Comparison
  f_test_pval = numeric()
)

# ------------------------------------------------------------------------------
# Fit Region-Wise Models
# ------------------------------------------------------------------------------
regions <- unique(data$region)

for (region_id in regions) {
  
  # Subset data for the current region
  region_data <- data[data$region == region_id, ]
  
  # ---- Model 1: Linear ----
  # Using REML = FALSE for Likelihood Ratio Testing
  lme_linear <- lmer(value ~ age + sex + (1 | hemi), data = region_data, REML = FALSE)
  lme_linear_sum <- coef(summary(lme_linear))
  
  # ---- Model 2: Quadratic ----
  lme_quad <- lmer(value ~ age + I(age^2) + sex + (1 | hemi), data = region_data, REML = FALSE)
  lme_quad_sum <- coef(summary(lme_quad))
  
  # ---- Append results ----
  result_list$region <- c(result_list$region, region_id)
  
  # Linear stats
  result_list$lme_intercept    <- c(result_list$lme_intercept,    lme_linear_sum["(Intercept)", "Estimate"])
  result_list$lme_age_estimate <- c(result_list$lme_age_estimate, lme_linear_sum["age", "Estimate"])
  result_list$lme_age_t        <- c(result_list$lme_age_t,        lme_linear_sum["age", "t value"])
  result_list$lme_age_p        <- c(result_list$lme_age_p,        lme_linear_sum["age", "Pr(>|t|)"])
  result_list$aic_lme_linear   <- c(result_list$aic_lme_linear,   AIC(lme_linear))
  result_list$bic_lme_linear   <- c(result_list$bic_lme_linear,   BIC(lme_linear))
  
  # Quadratic stats
  result_list$lme2_intercept     <- c(result_list$lme2_intercept,     lme_quad_sum["(Intercept)", "Estimate"])
  result_list$lme2_age_estimate  <- c(result_list$lme2_age_estimate,  lme_quad_sum["age", "Estimate"])
  result_list$lme2_age_t         <- c(result_list$lme2_age_t,         lme_quad_sum["age", "t value"])
  result_list$lme2_age_p         <- c(result_list$lme2_age_p,         lme_quad_sum["age", "Pr(>|t|)"])
  result_list$lme2_age2_estimate <- c(result_list$lme2_age2_estimate, lme_quad_sum["I(age^2)", "Estimate"])
  result_list$lme2_age2_t        <- c(result_list$lme2_age2_t,        lme_quad_sum["I(age^2)", "t value"])
  result_list$lme2_age2_p        <- c(result_list$lme2_age2_p,        lme_quad_sum["I(age^2)", "Pr(>|t|)"])
  result_list$aic_lme_quadratic  <- c(result_list$aic_lme_quadratic,  AIC(lme_quad))
  result_list$bic_lme_quadratic  <- c(result_list$bic_lme_quadratic,  BIC(lme_quad))
  
  # Likelihood Ratio Test
  lrt <- anova(lme_linear, lme_quad)
  result_list$f_test_pval <- c(result_list$f_test_pval, lrt$Pr[2])
}

# ------------------------------------------------------------------------------
# Post-Processing & FDR 
# ------------------------------------------------------------------------------
df <- as.data.frame(result_list)

# Multiple comparison corrections (FDR)
df$lme_age_fdr   <- p.adjust(df$lme_age_p, method = "fdr")
df$lme2_age_fdr  <- p.adjust(df$lme2_age_p, method = "fdr")
df$lme2_age2_fdr <- p.adjust(df$lme2_age2_p, method = "fdr")
df$f_test_fdr    <- p.adjust(df$f_test_pval, method = "fdr")

# Calculate Delta BIC (Positive = Quadratic has lower/better BIC)
df$delta_bic <- df$bic_lme_linear - df$bic_lme_quadratic

output_path <- here("regionwise_final_model_decisions.csv")
write.csv(df, output_path, row.names = FALSE)
cat("\nResults saved to:", output_path, "\n")