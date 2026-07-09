# Linear GAM vs Nonlinear (smooth) GAM, per region -- BIC as the sole decision rule.
# Both models fit identically in mgcv except how age enters:
#   linear:    value ~ age          + sex + s(hemi, re)   (straight line)
#   nonlinear: value ~ s(age, k=5)  + sex + s(hemi, re)   (penalized smooth)
# BIC is the most conservative criterion against added nonlinearity.

library(dplyr)
library(mgcv)
library(ggplot2)
library(here)

here::i_am("Supplementary/gam_model_comparison.R")  # match your actual filename

data <- read.csv(here("MIND_Network", "similarity_strength_subject_data.csv"))


res <- list(
  region        = character(),
  region_name   = character(),
  edf_age       = numeric(),   # kept for reference only; NOT used in the verdict
  bic_linear    = numeric(),
  bic_nonlinear = numeric()
)

# collect fitted curves for plotting
pred_list <- list()

for (region_id in unique(data$region)) {
  
  rd <- filter(data, region == region_id)
  rd$hemi <- as.factor(rd$hemi)
  rd$sex  <- as.factor(rd$sex)
  
  gam_linear    <- gam(value ~ age         + sex + s(hemi, bs = "re"), data = rd, method = "ML")
  gam_nonlinear <- gam(value ~ s(age, k=5) + sex + s(hemi, bs = "re"), data = rd, method = "ML")
  
  res$region        <- c(res$region,        region_id)
  res$region_name   <- c(res$region_name, rd$D99_abbr[1])
  res$edf_age        <- c(res$edf_age,        summary(gam_nonlinear)$s.table["s(age)", "edf"])
  res$bic_linear     <- c(res$bic_linear,     BIC(gam_linear))
  res$bic_nonlinear  <- c(res$bic_nonlinear,  BIC(gam_nonlinear))
  
  # --- fitted curves: population age trend, sex at reference, RE excluded ---
  age_grid <- seq(min(rd$age), max(rd$age), length.out = 100)
  nd <- data.frame(age  = age_grid,
                   sex  = factor(levels(rd$sex)[1], levels = levels(rd$sex)),
                   hemi = factor(levels(rd$hemi)[1], levels = levels(rd$hemi)))
  
  pl <- predict(gam_linear,    nd, exclude = "s(hemi)", se.fit = TRUE)
  pn <- predict(gam_nonlinear, nd, exclude = "s(hemi)", se.fit = TRUE)
  
  pred_list[[as.character(region_id)]] <- rbind(
    data.frame(region = region_id, age = age_grid, fit = pl$fit, se = pl$se.fit, model = "Linear"),
    data.frame(region = region_id, age = age_grid, fit = pn$fit, se = pn$se.fit, model = "Nonlinear")
  )
}

df <- as.data.frame(res)

# ---- BIC-only verdict ----
df$delta_bic <- df$bic_linear - df$bic_nonlinear        # positive => smooth preferred
df$verdict   <- ifelse(df$delta_bic > 2, "Nonlinear", "Linear")

cat("Verdict by BIC alone:\n"); print(table(df$verdict))
cat(sprintf("\nRegions where BIC prefers nonlinear: %d / %d   |   median delta_BIC: %.2f\n",
            sum(df$verdict == "Nonlinear"), nrow(df), median(df$delta_bic)))

write.csv(df,  here("Supplementary", "regionwise_bic_linear_vs_nonlinear_gam.csv"), row.names = FALSE)

# ------------------------
# Plot fitted curves for selected regions
# ------------------------
preds_all <- bind_rows(pred_list)

# pick the regions where BIC MOST favors the nonlinear fit (largest delta_bic)
n_show <- sum(df$delta_bic > 2)
top <- df %>% arrange(desc(delta_bic)) %>% slice(1:n_show)
target_regions <- as.character(top$region)

# make region the same type on both sides so the filter matches all of them
preds_all$region <- as.character(preds_all$region)
data$region      <- as.character(data$region)

# look up D99_abbr per region (assumes one abbr per region)
region_abbr <- data %>%
  distinct(region, D99_abbr) %>%
  filter(region %in% target_regions)

# join abbr onto top, preserving the delta_bic ordering
top <- top %>%
  mutate(region = as.character(region)) %>%
  left_join(region_abbr, by = "region")

# facet labels show D99_abbr + delta_BIC so you can read the strength
region_labeller <- setNames(
  sprintf("%s  (dBIC=%.1f)", top$D99_abbr, top$delta_bic),
  top$region
)

cat("\nTop regions by delta_BIC (BIC prefers nonlinear when > 0):\n")
print(top %>% select(region, D99_abbr, edf_age, delta_bic))
cat(sprintf("\nRegions selected: %d | with predictions available: %d\n",
            length(target_regions), length(intersect(target_regions, unique(preds_all$region)))))

plot_preds <- filter(preds_all, region %in% target_regions)
plot_pts   <- filter(data,      region %in% target_regions)

p <- ggplot() +
  geom_point(data = plot_pts, aes(age, value), color = "gray60", alpha = 0.5) +
  geom_ribbon(data = plot_preds,
              aes(age, ymin = fit - 1.96*se, ymax = fit + 1.96*se, fill = model),
              alpha = 0.15) +
  geom_line(data = plot_preds, aes(age, fit, color = model), size = 1) +
  facet_wrap(~region, scales = "free_y", ncol = 5, nrow = 9,
             labeller = labeller(region = region_labeller)) +
  scale_color_manual(name = "GAM fit", values = c("Linear" = "#1f77b4", "Nonlinear" = "#2ca02c")) +
  scale_fill_manual(name  = "GAM fit", values = c("Linear" = "#1f77b4", "Nonlinear" = "#2ca02c")) +
  theme_minimal(base_size = 12) +
  labs(x = "Age (years)", y = "Similarity strength") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(here("Supplementary", "gam_linear_vs_nonlinear_fits.png"),
       p, width = 12, height = 6, dpi = 200)
print(p)