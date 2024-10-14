library(data.table)
library(splines)  # For splines
library(lme4)  # For linear mixed-effects models
library(ggeffects)  # For creating predicted plots
library(ggplot2)
library(viridis)

### Load Data ###

covariate_data <- fread("~/Data/covariate_data.csv") # Load covariate data (e.g., demographics, comorbidities)
crp_centile_data <- fread("~/Data/crp_centile_data.csv") # Load CRP centile data (e.g., CRP values, converted CRP centile values)
abx_rx_data <- fread("~/Data/abx_rx_data.csv") # Load antibiotic prescription data (e.g., antibiotic rank, route, number of antibiotics, time of prescription)

### Compute Rx Decisions ###

source("computeRxDecision.R")
computeRxDecision(abx_rx_data, AbxRank, Route, NoOfAbx)
abx_rxdecision <- abx_rx_data[, .(EpisodeID, RxNo, RxDay = time_start, time_start, time_end, RxDecision)]
# Retain Non-NA Rx Decisions
abx_rxdecision <- abx_rxdecision[!is.na(RxDecision)]

### Match CRP Centile Data around Rx Window ###

abx_rxdecision[, `:=`(join_time_start = time_start - 0.5, join_time_end = time_end)]
rxdecision_crp <- crp_centile_data[abx_rxdecision, on = .(EpisodeID, Time >= join_time_start, Time < join_time_end)]
rxdecision_crp <- rxdecision_crp[!is.na(CRP)]
rxdecision_crp[, TimeFromRxDay := Time - time_start]

### Join with Covariates ###

rxdecision_crp <- rxdecision_crp[, .(EpisodeID, RxNo, RxDay, TimeFromRxDay, Centile, RxDecision)]
rxdecision_crp <- covariate_data[rxdecision_crp, on = .(EpisodeID)]

# Reassign RxIDs
rxdecision_crp[, RxID := .GRP, by = .(EpisodeID, RxNo)]

### Apply Exclusions ###

rxdecision_crp <- rxdecision_crp[TimeFromRxDay <= 2 & RxDay >= 0]
# only keep episodes with any(TimeFromRxDay > 0), keep all records for those episodes
rxdecision_crp <- rxdecision_crp[, .SD[any(TimeFromRxDay > 0)], by = RxID]

### Modelling and Prediction ###

# Define Splines
TimeFromRxDayBoundary <- quantile(rxdecision_crp$TimeFromRxDay, c(0.05, 0.95))
TimeFromRxDaySpline <- ns(rxdecision_crp$TimeFromRxDay, knots = c(1), Boundary.knots = TimeFromRxDayBoundary)

RxDaySplineBoundary <- quantile(rxdecision_crp$RxDay, c(0.05, 0.95))
RxDaySpline <- ns(rxdecision_crp$RxDay, knots = c(1, 3, 5), Boundary.knots = RxDaySplineBoundary)

# Fit Linear Mixed-Effects Model
RE_Centile <- lmer(
    Centile ~
        TimeFromRxDaySpline * RxDecision + RxDaySpline * TimeFromRxDaySpline + RxDaySpline * RxDecision +
        Age_scaled + Sex_int + Charlson + Elixhauser + EndRenalDialysis + Diabetes + NEWS_base +
        ImmunoSuppression + PalliativeCare + CommunityOnset + EpisodeBugGroup + SourceByIndication +
        (1 + TimeFromRxDaySpline | RxID),
        data = rxdecision_crp,
        REML = TRUE)


# Create Plot of CRP Centile Predictions over Time
Centile_Time_Rx_RE <- ggpredict(RE_Centile, terms = c("TimeFromRxDay [all]", "RxDecision", "RxDay [1, 2, 3, 4, 5]"))
Centile_Time_Rx_RE$group <- factor(Centile_Time_Rx_RE$group, levels = c("Escalation", "Unchanged", "De-escalation", "Stop"))

# Assign Facet Categories
Centile_Time_Rx_RE$facet_cat <- dplyr::case_when(
  facet == 1~ "Prescribing Day: 1",
  facet == 2~ "Prescribing Day: 2",
  facet == 3~ "Prescribing Day: 3",
  facet == 4~ "Prescribing Day: 4",
  facet == 5~ "Prescribing Day: 5"
)

p_rx_centile_time <- ggplot(Centile_Time_Rx_RE, aes(x = x, y = predicted, color = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  scale_x_continuous(breaks = seq(-0.5, 2, 0.5), labels = seq(-12, 48, 12)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Time (hours since prescribing decision)",
       y = "CRP Centile",
       color = "Prescribing Decision",
       fill = "Prescribing Decision") +
  theme_bw() +
  facet_grid(~facet_cat)

print(p_rx_centile_time)

