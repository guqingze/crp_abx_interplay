# Load Necessary Libraries
library(data.table)
library(ggplot2)
library(pROC)

# Load Data
covariate_data <- fread("~/Data/covariate_data.csv")
crp_centile_data <- fread("~/Data/crp_centile_data.csv")

# Preprocess CRP Centile Data
crp_centile_data <- crp_centile_data[Time >= 1 & Time <= 4, .(EpisodeID, Time, CRP, Centile)]
crp_centile_data[, CRP_log := log(CRP)]
crp_centile_data[, N := .N, by = EpisodeID]
crp_centile_data <- crp_centile_data[N >= 2][, N := NULL]

# Calculate Daily Changes using Linear Regression
calculate_daily_change <- function(data, var) {
  coef(lm(var ~ Time, data = data))["Time"]
}

crp_centile_data[, CentileChangePerDay := round(calculate_daily_change(.SD, Centile), 2), by = EpisodeID]
crp_centile_data[, CRPChangePerDay := round(calculate_daily_change(.SD, CRP), 2), by = EpisodeID]
crp_centile_data[, LogCRPChangePerDay := round(calculate_daily_change(.SD, CRP_log), 2), by = EpisodeID]

# Preprocess Covariate Data
covariate_data <- covariate_data[
  , `:=`(
    Age_scaled = Age / 10,
    Mortality30d = fifelse(is.na(DeathTime) | DeathTime > 30, 0L, 1L)
  )
]
covariate_data <- covariate_data[is.na(DeathTime) | DeathTime > 4]

# Merge Data
merged_data <- covariate_data[crp_centile_data, on = "EpisodeID"]
merged_data <- unique(merged_data)

# Split Data into Training and Testing Sets
set.seed(135)
merged_data <- merged_data[sample(.N)]
train_size <- floor(0.8 * nrow(merged_data))
train_set <- merged_data[1:train_size]
test_set <- merged_data[(train_size + 1):.N]

# Fit Logistic Regression Models
fit_logistic_model <- function(formula, data) {
  glm(formula, data = data, family = binomial)
}

models <- list(
  lr_covariate = fit_logistic_model(Mortality30d ~ Age_scaled + Sex + Charlson + Elixhauser + 
                                    EndRenalDialysis + Diabetes + NEWS_base + ImmunoSuppression +
                                    PalliativeCare + CommunityOnset + EpisodeBugGroup + SourceByIndication,
                                    train_set),
  
  lr_centile = fit_logistic_model(Mortality30d ~ ns(CentileChangePerDay, knots = c(-15, -7.5, 0, 7.5, 15), Boundary.knots = c(-25, 25)),
                                  train_set),
  
  lr_covariate_centile = fit_logistic_model(Mortality30d ~ ns(CentileChangePerDay, knots = c(-15, -7.5, 0, 7.5, 15), Boundary.knots = c(-25, 25)) +
                                             Age_scaled + Sex + Charlson + Elixhauser +
                                             EndRenalDialysis + Diabetes + NEWS_base + ImmunoSuppression +
                                             PalliativeCare + CommunityOnset + EpisodeBugGroup + SourceByIndication,
                                             train_set),

  lr_crp = fit_logistic_model(Mortality30d ~ CRPChangePerDay,
                              train_set),
                                             
  lr_covariate_crp = fit_logistic_model(Mortality30d ~ CRPChangePerDay + Age_scaled + Sex + 
                                         Charlson + Elixhauser + EndRenalDialysis + Diabetes +
                                         NEWS_base + ImmunoSuppression + PalliativeCare +
                                         CommunityOnset + EpisodeBugGroup + SourceByIndication,
                                         train_set),

  lr_crplog = fit_logistic_model(Mortality30d ~ LogCRPChangePerDay,
                                 train_set),

  lr_covariate_crplog = fit_logistic_model(Mortality30d ~ LogCRPChangePerDay + Age_scaled + Sex +
                                            Charlson + Elixhauser + EndRenalDialysis + Diabetes +
                                            NEWS_base + ImmunoSuppression + PalliativeCare +
                                            CommunityOnset + EpisodeBugGroup + SourceByIndication,
                                            train_set)
)

# Evaluate Model Performance
calculate_auc <- function(model, data) {
  data[, PredictedProb := predict(model, newdata = data, type = "response")]
  roc_obj <- roc(data$Mortality30d, data$PredictedProb)
  list(auc = auc(roc_obj), ci = ci.auc(roc_obj))
}

auc_results <- lapply(models, calculate_auc, data = test_set)

# Create and Print AUC Summary
auc_summary <- data.table(
  Model = names(auc_results),
  AUC = sapply(auc_results, function(x) x$auc),
  LowerCI = sapply(auc_results, function(x) x$ci[1]),
  UpperCI = sapply(auc_results, function(x) x$ci[3])
)

print(auc_summary)

# Perform ROC Tests
roc_tests <- list(
  roc_test_crp_centile = roc.test(auc_results$lr_crp$roc_obj, auc_results$lr_centile$roc_obj),
  roc_test_crplog_centile = roc.test(auc_results$lr_crplog$roc_obj, auc_results$lr_centile$roc_obj),
  roc_test_covariate_centile = roc.test(auc_results$lr_covariate$roc_obj, auc_results$lr_covariate_centile$roc_obj),
)

lapply(roc_tests, print)