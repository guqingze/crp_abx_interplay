# mortality_early_crp.R: R script for analysis of the impact of early CRP centile changes on 5-30 day all-cause mortality.

# Load required packages
library(data.table)  # for efficient data manipulation
library(survival)    # for survival analysis
library(splines)     # for natural splines in regression
library(ggplot2)     # for plotting
library(viridis)     # for color scales
library(ggsurvfit)   # for enhanced survival plots

# Load Data
covariate_data <- fread("~/Data/covariate_data.csv")
crp_centile_data <- fread("~/Data/crp_centile_data.csv")

# Preprocess CRP Centile Data
crp_centile_data <- crp_centile_data[Time >= 1 & Time <= 4, .(EpisodeID, Time, CRP, Centile)]
crp_centile_data[, N := .N, by = EpisodeID]
crp_centile_data <- crp_centile_data[N >= 2][, N := NULL]

# Calculate Daily Changes using Linear Regression
calculate_daily_change <- function(data, var) {
  coef(lm(var ~ Time, data = data))["Time"]
}

# Calculate the centile change per day for each EpisodeID and add it to the dataframe
crp_centile_data[, CentileChangePerDay := round(calculate_daily_change(.SD, Centile), 2), by = EpisodeID]

# Classify centile changes into ranges
crp_centile_data[, CentileChangeRange := fcase(
  CentileChangePerDay < -5, "< -5",
  CentileChangePerDay >= -5 & CentileChangePerDay < 0, "[-5, 0)",
  CentileChangePerDay >= 0 & CentileChangePerDay < 10, "[0, 10)",
  CentileChangePerDay >= 10 & CentileChangePerDay < 20, "[10, 20)",
  CentileChangePerDay >= 20, "≥ 20"
)]

# Convert the centile change range to a factor with specified levels
crp_centile_data$CentileChangeRange <- factor(crp_centile_data$CentileChangeRange, 
                                          levels = c("< -5", "[-5, 0)", "[0, 10)", "[10, 20)", "≥ 20"))

# Preprocess Covariate Data
covariate_data <- covariate_data[
  , `:=`(
    Age_scaled = Age / 10,
    Mortality30d = fifelse(is.na(DeathTime) | DeathTime > 30, 0L, 1L)
  )
]
covariate_data <- covariate_data[is.na(DeathTime) | DeathTime > 4]

# Merge Data
crp_centile_data <- covariate_data[crp_centile_data, on = "EpisodeID"]
crp_centile_data <- unique(crp_centile_data)

# Censor DeathTime at 30 days
crp_centile_data[, DeathTime_censored := fifelse(is.na(DeathTime) | DeathTime > 30, 30, DeathTime)]

# Fit survival curves based on centile change range
surv_fit <- survfit2(Surv(DeathTime_censored, Mortality30d) ~ CentileChangeRange, data = crp_centile_data)

# Plot Kaplain-Meier survival curves
p_km_surv <- 
  ggsurvfit(surv_fit, linewidth = 1) +
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({n.event})", size = 4, risktable_height = 0.14, times = c(0, 5, 10, 15, 20, 25, 30)) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 15) +
  scale_color_viridis(discrete = TRUE, direction = -1) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  scale_x_continuous(limits = c(0, 30), breaks = c(4, seq(0, 30, 5)), minor_breaks = seq(7.5, 27.5, 2.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  geom_vline(xintercept = 4, linetype = "dashed", color = "red") +
  guides(color = guide_legend(title = "Centile Change Per Day During Days 1-4: "),
         fill = guide_legend(title = "Centile Change Per Day During Days 1-4: ")) +
  labs(x = "Time (days after blood culture collection)", y = "Unadjusted Survival Probability") + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 18))

print(p_km_surv)

# Scale centile change per day for spline modeling
crp_centile_data[, CentileChangePerDay_scaled := CentileChangePerDay / 10]

# Create natural spline terms for centile change per day
CentileChangePerDaySpline <- ns(
  crp_centile_data$CentileChangePerDay_scaled,
  knots = c(-1.5, -0.75, 0, 0.75, 1.5),
  Boundary.knots = c(-2.5, 2.5))

# Fit logistic regression model for 30-day mortality
lr_mortality_centile <- glm(
  Mortality30d ~ 
  CentileChangePerDaySpline +
  Age_scaled + Sex + Charlson + Elixhauser + EndRenalDialysis + Diabetes +
  NEWS_base + ImmunoSuppression + PalliativeCare + CommunityOnset + EpisodeBugGroup + SourceByIndication,
  data = crp_centile_data, family = binomial)

# Predict mortality probability based on the centile change spline
Mortality30d_CentileChange <- ggpredict(lr_mortality_centile, terms = "CentileChangePerDaySpline")

# Plot relationship between centile change and mortality probability
plot_mortality_prediction <- 
  ggplot(Mortality30d_CentileChange, aes(x = x, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1), labels = seq(-30, 30, 10)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 1, 0.01), labels = scales::percent) +
  labs(x = "CRP centile change per day during days 1-4", y = "Adjusted Probability of 4-30-day mortality") +
  theme_bw() +
  theme(text = element_text(size = 18))

print(plot_mortality_prediction)