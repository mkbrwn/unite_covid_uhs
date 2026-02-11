# This aims to perform full case predictive modelling on the UNTIE COVID data set

# run the initial data processing script to load and prepare the data for training/test set
source("src/01_raw_to_pre_imputation.R")

# libraries needed for selection and tables
library(MASS)        # stepAIC
library(gtsummary)   # tbl_regression
library(dplyr)

# Scale the continuous variables in the training data set
train_data_scaled = train_data %>%
  mutate(across(where(is.numeric) & !ends_with("_YN"), scale))

# univariable analysis selecting variable P<0.2 

vars <- c(
  "NEW_SITE_ID", "INC_AGE_INT", "INC_LOS_PRIOR_ADM_INT", "INC_INTERVAL_HOSPIT_INT",
  "INC_ADM_SURGE_BED_YN", "INC_CARDIAC_DISEASE_YN", "INC_LIVER_DISEASE_YN", "INC_HBP_YN",
  "INC_NEURO_YN", "INC_PULMO_DISEASE_YN", "INC_ASTHMA_YN", "INC_NEOPLASM_YN",
  "INC_KIDNEY_DISEASE_YN", "INC_IMMUNOSUPPR_YN", "INC_ACE_INHIB_YN", "INC_ANGIO_II_YN",
  "INC_ANTIPLAT_YN", "ICU_RESP_SUPPORT_YN", "ICU_WHITE_CELL_INT", "ICU_NEUTRO_INT",
  "ICU_CRP_INT", "ICU_PLATELETS_INT", "ICU_CARDIAC_THERAPY_YN", "ICU_SEPSIS_YN",
  "ICU_STRESS_MYOC_YN", "ICU_MYOCARDITIS_YN", "ICU_PERICARD_YN", "ICU_PNEUMOTHORAX_YN",
  "ICU_ATELECTASIS_YN", "ICU_DELIRIUM_YN", "ICU_PRESSURE_OTH_YN", "ICU_KIDNEY_INJ_YN",
  "ICU_OBSTRUCTION_YN", "ICU_CORTICO_YN", "ICU_CLIN_TRIAL_YN", "ICU_SEDATION_YN",
  "ICU_SEDAT_DURATION_INT", "ICU_RRT_DIAL_YN", "ICU_BLOOD_PURIF_YN", "ICU_INOTROPES_YN",
  "ICU_TRACHEOS_YN", "RESP_INTUBATED_YN", "RESP_INTUBATED_ICU_STAY_YN", "RESP_NI_VENT_YN",
  "RESP_HFNC_YN", "RESP_INV_VENT_YN", "RESP_DURATION_INV_VENT_INT", "RESP_ECMO_YN",
  "RESP_PRONE_YN", "RESP_NEUROM_BLOC_YN", "RESP_VENT_ROUTINE_YN", "INF_ANTIBIO_YN",
  "INF_ANTIFUNG_YN", "INF_AT_ADMISSION_YN", "INF_DURING_ICU_YN", "NEW_BMI",
  "wave", "ventilation_severity", "INC_DIABETES1_YN", "comorbidity_score",
  "neutrophil_lymphocyte_ratio", "ICU_CRP_CATEGORY", "OUT_DEAD_DURING_ICU_YN"
)
# Keep only the variables listed in vars that exist in the scaled data
vars_present <- intersect(vars, names(train_data_scaled))
if (length(vars_present) == 0) stop("None of the variables in 'vars' are present in train_data_scaled")
missing_vars <- setdiff(vars, vars_present)
if (length(missing_vars) > 0) message("The following vars were not found and will be skipped: ", paste(missing_vars, collapse = ", "))
train_data_scaled <- train_data_scaled %>% dplyr::select(dplyr::all_of(vars_present))

# Ensure outcome is present
if (!"OUT_DEAD_DURING_ICU_YN" %in% names(train_data_scaled)) stop("Outcome OUT_DEAD_DURING_ICU_YN is not present after selection")

# Ensure expected categorical variables are factors (keep YN as numeric 0/1)
factor_vars <- c("NEW_SITE_ID", "wave", "ventilation_severity", "ICU_CRP_CATEGORY")
for (fv in factor_vars) {
  if (fv %in% names(train_data_scaled)) {
    train_data_scaled[[fv]] <- as.factor(train_data_scaled[[fv]])
  }
}

# Drop rows with any missing values (full-case analysis)
train_data_scaled <- train_data_scaled %>% tidyr::drop_na()

# Drop constant predictors (zero variance), excluding the outcome
outcome_var <- "OUT_DEAD_DURING_ICU_YN"
const_preds <- names(train_data_scaled)[sapply(train_data_scaled, function(x) dplyr::n_distinct(x) <= 1)]
const_preds <- setdiff(const_preds, outcome_var)
if (length(const_preds) > 0) {
  message("Removing constant predictors: ", paste(const_preds, collapse = ", "))
  train_data_scaled <- train_data_scaled %>% dplyr::select(-dplyr::all_of(const_preds))
}

# Build full model formula
predictors <- setdiff(names(train_data_scaled), outcome_var)
full_formula <- as.formula(paste(outcome_var, paste(predictors, collapse = " + "), sep = " ~ "))

# Fit full logistic regression model
full_fit <- glm(full_formula, data = train_data_scaled, family = binomial(link = "logit"), control = list(maxit = 100))

# Backward stepwise selection by AIC
final_fit <- suppressWarnings(stepAIC(full_fit, direction = "backward", trace = FALSE))

cat("\nFinal model AIC:", AIC(final_fit), "\n")
cat("Selected terms:\n")
print(attr(terms(final_fit), "term.labels"))

# Summaries
print(summary(final_fit))

# Nicely formatted regression table with ORs
tbl_final <- tbl_regression(final_fit, exponentiate = TRUE)
print(tbl_final)

# Optionally, save selected variables for downstream use
selected_vars <- attr(terms(final_fit), "term.labels")
readr::write_csv(tibble::tibble(variable = selected_vars), "data/processed/selected_vars_backward_aic.csv")

# final model object is `final_fit`, predictors listed in `selected_vars`