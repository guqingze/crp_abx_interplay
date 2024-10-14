library(data.table)
library(nnet)
library(splines)
library(ggplot2)
library(patchwork)
library(viridis)

### Load Data ###

covariate_data <- fread("~/Data/covariate_data.csv") # Load covariate data (e.g., demographics, comorbidities)
crp_centile_data <- fread("~/Data/crp_centile_data.csv") # Load CRP centile data (e.g., CRP values, converted CRP centile values)
abx_rx_data <- fread("~/Data/abx_rx_data.csv") # Load antibiotic prescription data (e.g., antibiotic rank, route, number of antibiotics, time of prescription)

### Preprocess CRP Centile Data ###

# Generate CRP Window Data
crp_windows <- unique(crp_centile_data[, .(EpisodeID)])
crp_windows <- crp_windows[rep(1:.N, each = 8)][, `:=`(
  time_start = rep(0:7, .N / 8),
  time_end = rep(1:8, .N / 8)
)]

# Join CRP Data with Created Window
crp_centile_data[, join_time := Time]
crp_centile_data <- crp_centile_data[crp_windows, on = .(EpisodeID, join_time >= time_start, join_time < time_end)]

# Organize into Daily Data with Calculations
crp_centile_daily <- crp_centile_data[, .(
  EpisodeID, Day = time_start, CRP, Centile, Time, time_start, time_end
)]

# Calculate Daily First and Last Values
crp_centile_daily[, `:=`(
  Centile_FirstInDay = first(Centile),
  Time_FirstInDay = first(Time),
  Centile_LastInDay = last(Centile),
  Time_LastInDay = last(Time)
), by = .(EpisodeID, Day)]

# Ensure Unique Records and Calculate Changes
crp_centile_daily <- unique(crp_centile_daily, by = c("ClusterID", "EpisodeID", "Day"))
crp_centile_daily[, `:=`(
  CentileChange = Centile_LastInDay - shift(Centile_FirstInDay, type = "lag"),
  Time_diff = Time_LastInDay - shift(Time_FirstInDay, type = "lag")
), by = EpisodeID]

# Filter Out NAs and Determine Recovery Status
crp_centile_daily <- crp_centile_daily[!is.na(CentileChange)]
crp_centile_daily[, RecoveryStatus := fcase(
  abs(CentileChange) <= 15, "Recovering as expected",
  CentileChange > 15, "Sub-optimal recovery",
  CentileChange < -15, "Recovering faster than expected"
)]
crp_centile_daily[, RecoveryStatus := factor(RecoveryStatus, levels = c("Recovering faster than expected", "Recovering as expected", "Sub-optimal recovery"))]

### Join CRP Centile Data and Antibiotic Rx Data ###

crp_centile_daily <- crp_centile_daily[, .(
  EpisodeID, CRPNo = seq_len(.N), CentileChange, RecoveryStatus, 
  Time_crp2 = Time, Time_crp2_1d = Time + 1
), by = .(EpisodeID)]

# Joins: Before and After CRP Event
crp_centile_daily[, join_time := Time_crp2]
crp_rx_before <- abx_rx_data[crp_centile_daily, on = .(EpisodeID,
                                                 time_start <= join_time,
                                                 time_end > join_time)][, RxOder := "Before"]
crp_centile_daily[, join_time := Time_crp2_1d]
crp_rx_after <- abx_rx_data[crp_centile_daily, on = .(EpisodeID,
                                                time_start < join_time,
                                                time_end >= join_time)][, RxOder := "After"]

crp_rx_data <- rbind(crp_rx_before, crp_rx_after)
crp_rx_data <- crp_rx_data[order(EpisodeID, CRPNo)]

### Compute Rx Decisions ###

source("computeRxDecision.R")
computeRxDecision(crp_rx_data, AbxRank, Route, NoOfAbx)
crp_rxdecision <- crp_rx_data[!is.na(RxDecision)]

### Join Covariates ###

crp_rxdecision <- crp_rxdecision[, .(EpisodeID, CentileChange, RxDecision, RxDay = Time_crp2_1d, RecoveryStatus]
crp_rxdecision <- covariate_data[crp_rxdecision, on = "EpisodeID"]

crp_rxdecision[, CentileChange_trun := fcase(
  CentileChange > 30, 30,
  CentileChange >= -30 & CentileChange <= 30, CentileChange,
  CentileChange < -30, -30)]

### Plot: Crude percentage of prescribing decisions after CRP centile change ###

plot_bar_rx_after_centile <-
  ggplot(crp_rxdecision, aes(x = CentileChange_trun, fill = RxDecision, color = RxDecision)) +
  geom_histogram(position = "fill", binwidth = 5, alpha = 0.9, boundary = 0) +
  # geom_vline(xintercept = c(-15, 15), linetype = "dashed", color = "#D8152F", linewidth = 1) +
  scale_x_continuous(breaks = seq(-50, 50, 5)) +
  scale_y_continuous(labels = scales::percent_format(), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Centile change over the past two consecutive days",
    y = "Percentage",
    fill = "Prescribing Decision",
    color = "Prescribing Decision") +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  theme(text = element_text(size = 18),
        legend.position = "right")

plot_freqhist_rx_after_centile <- 
  ggplot(crp_rxdecision, aes(x = CentileChange_trun, fill = RecoveryStatus)) +
  geom_histogram(binwidth = 5, boundary = 0, alpha = 0.7) +
  scale_x_continuous(breaks = seq(-50, 50, 5), limits = c(-30, 30)) +
  scale_fill_manual(
    values = c("Recovering faster than expected" = "#FF9900",
               "Recovering as expected" = "#0073C2",
               "Sub-optimal recovery" = "#D8152F"),
    labels = c("Recovering faster than expected",
               "Recovering as expected",
               "Sub-optimal recovery")) +
  guides(fill = guide_legend(title = "Recovery Status")) +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.title.x=element_blank(),
        legend.position = "top",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

plot_bar_freq_rx_after_centile <- 
  plot_freqhist_rx_after_centile / plot_bar_rx_after_centile +
  plot_layout(heights = c(1, 5))

print(plot_bar_freq_rx_after_centile)

### Modelling and Predictions ###

# Prepare Splines for Modeling
CentileChangeSpline <- ns(crp_rxdecision$CentileChange_trun, knots = c(-15, 0, 15), Boundary.knots = c(-25, 25))

RxDaySplineBoundary <- quantile(crp_rxdecision$RxDay, c(0.05, 0.95))
RxDaySpline <- ns(crp_rxdecision$RxDay, knots = c(3, 5, 7), Boundary.knots = RxDaySplineBoundary)

crp_rxdecision$RxDecision <- factor(crp_rxdecision$RxDecision, levels = c("Unchanged", "Escalation", "De-escalation", "Stop"))

# Multinomial Logistic Regression Model
MLR_RxDecision <- multinom(
    RxDecision ~ 
        CentileChangeSpline * RxDaySpline +
        Age_scaled + Sex + Charlson + Elixhauser + EndRenalDialysis + Diabetes + NEWS_base +
        ImmunoSuppression + PalliativeCare  + CommunityOnset + EpisodeBugGroup + SourceByIndication,
        data = crp_rxdecision, maxit = 1000)

# Generate Prediction Data for CentileChange
RxDecision_Centile <- as.data.table(
  expand.grid(
    CentileChange_trun = seq(-30, 30, 1),
    RxDay = 4,
    Age_scaled = 6.7,
    Sex = "M",
    Charlson = 2,
    Elixhauser = 3,
    EndRenalDialysis = "0",
    Diabetes = "0",
    NEWS_base = 3,
    ImmunoSuppression = "0",
    PalliativeCare = "0",
    CommunityOnset = "1",
    EpisodeBugGroup = "E. coli",
    SourceByIndication = "Urinary"))

# Predict and Plot RxDecision ~ CentileChange
RxDecision_Centile_pred <- predictions(MLR_RxDecision, newdata = RxDecision_Centile, type = "probs")
RxDecision_Centile_pred$group <- factor(RxDecision_Centile_pred$group, levels = c("Escalation", "Unchanged", "De-escalation", "Stop"))

p_rxdecision_centilechange <-
  ggplot(RxDecision_Centile_pred, aes(CentileChange_trun, estimate)) +
  geom_line(aes(color = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .1) +
  scale_x_continuous(breaks = seq(-30, 30, 5)) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(), breaks = seq(0,1,0.2)) +
  labs(x = "Centile change over the past two consecutive days",
       y = "Probability",
       color = "Prescribing Decision",
       fill = "Prescribing Decision") +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.title.y=element_text(size=18))

# Generate Prediction Data for RxDay
RxDecision_RxDay <- as.data.table(
  expand.grid(
    CentileChange_trun = 0,
    RxDay = seq(2, 8.99, 0.1),
    Age_scaled = 6.7,
    Sex = "M",
    Charlson = 2,
    Elixhauser = 3,
    EndRenalDialysis = "0",
    Diabetes = "0",
    NEWS_base = 3,
    ImmunoSuppression = "0",
    PalliativeCare = "0",
    CommunityOnset = "1",
    EpisodeBugGroup = "E. coli",
    SourceByIndication = "Urinary"))

# Predict and Plot RxDecision ~ RxDay
RxDecision_RxDay_pred <- predictions(MLR_RxDecision, newdata = RxDecision_RxDay, type = "probs")
RxDecision_RxDay_pred$group <- factor(RxDecision_RxDay_pred$group, levels = c("Escalation", "Unchanged", "De-escalation", "Stop"))

p_rxdecision_rxday <-
  ggplot(RxDecision_RxDay_pred, aes(RxDay, estimate)) +
  geom_line(aes(color = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .1) +
  scale_x_continuous(breaks = seq(0, 9, 1)) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(), breaks = seq(0,1,0.2)) +
  labs(x = "Time (days after blood culture collection)",
       y = "Probability",
       color = "Prescribing Decision",
       fill = "Prescribing Decision") +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 18),
        axis.title.y=element_text(size=18))

# Plot Outputs
print(p_rxdecision_centilechange)
print(p_rxdecision_rxday)