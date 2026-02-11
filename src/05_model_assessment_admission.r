#This script is for assessment of the predictive modelling from the UNTIE COVID data set
# The 04_predictive_model.r script produces a multivariable model from backwards stepwise selection after muiltiple imputation.

# Matt Strammer says ~ forests, calibration and clinical utility.

# Load required packages
library(caret)
library(mice)
library(miceadds)
library(future) # parrallelisation of MICE
library(performance)
library(MASS)
library(pROC)
library(openxlsx)


# source of training and test set
source("src/04_training_predictive_model.r")
source("src/03_imputation_test_data.r")

#construct formula's for the model analysis
formula_cheat <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
  ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN +
  ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
  ICU_RRT_DIAL_YN + ICU_INOTROPES_YN +  RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT +
  RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
  neutrophil_lymphocyte_ratio + ICU_CRP_INT + age_group

model_formula_without_CRP <- update(model_formula, . ~ . - ICU_CRP_INT)

#model contruction
model_admission <- glm(model_formula, data = complete(train_data_imputation), family = binomial(link = "logit"))
model_admission_witout_crp <- glm(model_formula_without_CRP, data = complete(train_data_imputation), family = binomial(link = "logit"))
model_cheating <- glm(formula_cheat, data = complete(train_data_imputation), family = binomial(link = "logit"))
model_simple <- glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CRP_INT*ICU_CORTICO_YN, data = complete(train_data_imputation), family = binomial(link = "logit"))

# produce a prediction on the test set
test_data_completed <- complete(test_data_imputation_scaled)
test_data_completed$predicted_prob_admission <- predict(model_admission, newdata = complete(test_data_imputation_scaled), type = "response")
test_data_completed$predicted_prob_admission_without_crp <- predict(model_admission_witout_crp, newdata = complete(test_data_imputation_scaled), type = "response")
test_data_completed$predicted_prob_cheating <- predict(model_cheating, newdata = complete(test_data_imputation_scaled), type = "response")
test_data_completed$predicted_prob_simple <- predict(model_simple, newdata = complete(test_data_imputation_scaled), type = "response")

# plot the current multivariable model from admission variables
model_admission %>% tbl_regression(exp = TRUE) %>% print()
#exort tbl_regression(exp = TRUE) as excel or csv for figure production
openxlsx::write.xlsx(model_admission %>% tbl_regression(exp = TRUE) %>% as_tibble(), "Figures/model_assessment/model_admission_tbl_regression.xlsx", rowNames = FALSE)

#########################################################
# Construct ROC curves for all models and save outputs
#########################################################

# Outcome vector
y_true <- test_data_completed$OUT_DEAD_DURING_ICU_YN

# Prepare predictions list
pred_list <- list(
  "Admission" = test_data_completed$predicted_prob_admission,
  "Admission_no_CRP" = test_data_completed$predicted_prob_admission_without_crp,
  "Cheating" = test_data_completed$predicted_prob_cheating,
  "Simple_CRPxCortico" = test_data_completed$predicted_prob_simple
)

# Directory for figures
fig_dir <- "Figures/model_assessment"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Compute ROC objects, AUCs, and 95% CIs
roc_objs <- list()
auc_vals <- numeric(length(pred_list))
names(auc_vals) <- names(pred_list)
auc_ci_lower <- numeric(length(pred_list))
auc_ci_upper <- numeric(length(pred_list))
names(auc_ci_lower) <- names(pred_list)
names(auc_ci_upper) <- names(pred_list)

for (nm in names(pred_list)) {
  preds <- pred_list[[nm]]
  roc_obj <- roc(response = y_true, predictor = preds, quiet = TRUE)
  roc_objs[[nm]] <- roc_obj
  auc_vals[nm] <- as.numeric(auc(roc_obj))
  ci_vec <- as.numeric(ci.auc(roc_obj))  # c(lower, auc, upper)
  auc_ci_lower[nm] <- ci_vec[1]
  auc_ci_upper[nm] <- ci_vec[3]

  # Save individual ROC plot
  file_safe <- gsub("[^A-Za-z0-9_]", "_", nm)
  png(file.path(fig_dir, paste0("roc_", file_safe, ".png")), width = 900, height = 700, res = 150)
  plot(roc_obj, col = "#1f77b4", lwd = 2,
       main = sprintf("%s\nAUC = %.3f (95%% CI %.3f–%.3f)", nm, auc_vals[nm], auc_ci_lower[nm], auc_ci_upper[nm]))
  abline(a = 0, b = 1, lty = 2, col = "gray60")
  dev.off()
}

# Combined ROC comparison plot
png(file.path(fig_dir, "roc_comparison_all_models.png"), width = 1000, height = 800, res = 150)
cols <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e")
i <- 1
for (nm in names(roc_objs)) {
  if (i == 1) {
    plot(roc_objs[[nm]], col = cols[i], lwd = 2, main = "ROC Comparison: Test Set")
  } else {
    plot(roc_objs[[nm]], col = cols[i], lwd = 2, add = TRUE)
  }
  i <- i + 1
}
abline(a = 0, b = 1, lty = 2, col = "gray60")
legend("bottomright",
       legend = sprintf("%s (AUC=%.3f; 95%% CI %.3f–%.3f)",
                        names(roc_objs),
                        auc_vals[names(roc_objs)],
                        auc_ci_lower[names(roc_objs)],
                        auc_ci_upper[names(roc_objs)]),
       col = cols[seq_along(roc_objs)], lwd = 2, bty = "n")
dev.off()

# Save AUC summary
auc_df <- data.frame(
  Model = names(auc_vals),
  AUC = round(auc_vals, 4),
  AUC_Lower95 = round(auc_ci_lower[names(auc_vals)], 4),
  AUC_Upper95 = round(auc_ci_upper[names(auc_vals)], 4),
  row.names = NULL
)
auc_df <- auc_df[order(-auc_df$AUC), ]
write.csv(auc_df, file.path(fig_dir, "model_auc_summary.csv"), row.names = FALSE)

cat("ROC analysis completed. Outputs in ", fig_dir, "\n")

#########################################################
# Add calibration metrics table per model
# - Brier score on test set
# - AUC with 95% CI (already computed)
# - Precision, Sensitivity (Recall), F1 at Youden threshold
#########################################################

# Helper to coerce outcome to 0/1 numeric robustly
to_binary01 <- function(v) {
  if (is.factor(v)) {
    lv <- levels(v)
    # Try direct numeric levels first ("0"/"1")
    suppressWarnings(num <- as.numeric(as.character(v)))
    if (!any(is.na(num))) return(as.integer(num))
    # Fallback: map highest level to 1, others to 0
    return(as.integer(v == lv[length(lv)]))
  } else if (is.logical(v)) {
    return(as.integer(v))
  } else {
    return(as.integer(v))
  }
}

y_bin <- to_binary01(y_true)

metrics <- lapply(names(pred_list), function(nm) {
  preds <- pred_list[[nm]]
  roc_obj <- roc_objs[[nm]]
  thr <- as.numeric(coords(roc_obj, "best", ret = "threshold", best.method = "youden"))
  cls <- as.integer(preds >= thr)
  TP <- sum(cls == 1 & y_bin == 1, na.rm = TRUE)
  FP <- sum(cls == 1 & y_bin == 0, na.rm = TRUE)
  FN <- sum(cls == 0 & y_bin == 1, na.rm = TRUE)
  TN <- sum(cls == 0 & y_bin == 0, na.rm = TRUE)
  precision <- if ((TP + FP) == 0) 0 else TP / (TP + FP)
  recall    <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
  f1        <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  # Sensitivity is equivalent to recall (TPR)
  sensitivity <- recall
  # Specificity (TNR)
  specificity <- if ((TN + FP) == 0) 0 else TN / (TN + FP)
  # Brier score on test set
  brier <- mean((preds - y_bin)^2, na.rm = TRUE)
  data.frame(
    Model = nm,
    AUC = round(auc_vals[nm], 4),
    AUC_Lower95 = round(auc_ci_lower[nm], 4),
    AUC_Upper95 = round(auc_ci_upper[nm], 4),
    Youden_Threshold = round(thr, 4),
    Precision = round(precision, 4),
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    Recall = round(recall, 4),
    F1 = round(f1, 4),
    Brier_Score = round(brier, 4),
    row.names = NULL
  )
})

metrics_df <- do.call(rbind, metrics)

# Order by AUC desc for readability
metrics_df <- metrics_df[order(-metrics_df$AUC), ]

# Save and print
metrics_csv <- file.path(fig_dir, "model_metrics_summary.csv")
write.csv(metrics_df, metrics_csv, row.names = FALSE)
cat("Saved model metrics summary to:", metrics_csv, "\n")
print(metrics_df, row.names = FALSE)

#comapre performance between models

print("Model performance comparison:")
compare_performance(
  model_admission,
  model_admission_witout_crp,
  model_cheating,
  model_simple,
  verbose = FALSE
)

