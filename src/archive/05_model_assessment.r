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

# source of training and test set 

source("src/01_raw_to_pre_imputation.r")

######################################  univariable and mice ######################################

### columns from univariable analysis 
vars <- c(
  "NEW_SITE_ID", "INC_AGE_INT", "INC_LOS_PRIOR_ADM_INT", "INC_INTERVAL_HOSPIT_INT",
  "INC_ADM_SURGE_BED_YN", "INC_CARDIAC_DISEASE_YN", "INC_LIVER_DISEASE_YN", "INC_HBP_YN",
  "INC_NEURO_YN", "INC_PULMO_DISEASE_YN", "INC_ASTHMA_YN", "INC_NEOPLASM_YN",
  "INC_KIDNEY_DISEASE_YN", "INC_IMMUNOSUPPR_YN", "INC_ACE_INHIB_YN", "INC_ANGIO_II_YN",
  "INC_ANTIPLAT_YN", "ICU_RESP_SUPPORT_YN", "ICU_WHITE_CELL_INT", "ICU_NEUTRO_INT",
  "ICU_CRP_INT", "ICU_PLATENT", "ICU_CARDIAC_THERAPY_YN", "ICU_SEPSIS_YN",
  "ICU_STRESS_MYOC_YN", "ICU_MYOCARDITIS_YN", "ICU_PERICARD_YN", "ICU_PNEUMOTHORAX_YN",
  "ICU_ATELECTASIS_YN", "ICU_DELIRIUM_YN", "ICU_PRESSURE_OTH_YN", "ICU_KIDNEY_INJ_YN",
  "formula_without_CRP_CORTICO_YN", "ICU_CLIN_TRIAL_YN", "ICU_SEDATION_YN",
  "ICU_SEDAT_DURATION_INT", "ICU_RRT_DIAL_YN", "ICU_BLOOD_PURIF_YN", "ICU_INOTROPES_YN",
  "RESP_INTUBATED_YN"model_univariable = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*ICU_CRP_INT, data = complete(train_data_imputed_scaled), family = binomial), "RESP_INTUBATED_ICU_STAY_YN", "RESP_NI_VENT_YN",
  "RESP_HFNC_YN", "RESP_INV_VENT_YN", "RESP_DURATION_INV_VENT_INT", "RESP_ECMO_YN",
  "RESP_PRONE_YN", "RESP_NEUROM_BLOC_YN", "RESP_VENT_ROUTINE_YN", "INF_ANTIBIO_YN",
  "INF_ANTIFUNG_YN", "INF_AT_ADMISSION_YN", "INF_DURING_ICU_YN", "NEW_BMI",
  "wave", "ventilation_severity", "INC_DIABETES1_YN", "comorbidity_score",
  "neutrophil_lymphocyte_ratio", "ICU_CRP_CATEGORY", "OUT_DEAD_DURING_ICU_YN"
)


# Subset the training data to only include columns identified from univariable analysis
Preimputation_train_data_muitltivariable = train_data |> dplyr::select(all_of(vars))

# impute the train set
train_data_imputation = futuremice(Preimputation_train_data_muitltivariable, m = 5, ncore = availableCores()-1, maxit = 10)

# Create age categories and convert back to mids object
# Step 1: Extract each imputed dataset and add age categories
modified_list <- lapply(1:train_data_imputation$m, function(i) {
  complete(train_data_imputation, i) %>%
    mutate(
      # Create age groups
      age_group = case_when(
        INC_AGE_INT < 40 ~ "<40",
        INC_AGE_INT >= 40 & INC_AGE_INT < 50 ~ "40-49",
        INC_AGE_INT >= 50 & INC_AGE_INT < 60 ~ "50-59",
        INC_AGE_INT >= 60 & INC_AGE_INT < 70 ~ "60-69",
        INC_AGE_INT >= 70 & INC_AGE_INT < 80 ~ "70-79",
        INC_AGE_INT >= 80 & INC_AGE_INT < 90 ~ "80-89",
        TRUE ~ "90+"
      )) %>%
      # Scale numeric variables and convert back to mids object
       mutate(across(where(is.numeric), ~ as.vector(scale(.)))) 
    })

# Step 2: Convert back to mids object using datlist2mids
train_data_imputation <- datlist2mids(modified_list)  

#formulat for the final imputed model 

formula <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
             ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
             ICU_RRT_DIAL_YN + ICU_INOTROPES_YN +  RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
             neutrophil_lymphocyte_ratio + ICU_CRP_INT + age_group

train_data_imputed_scaled_muiltivariable = glm(formula, data = complete(train_data_imputation), family = binomial(link = "logit"))

# extract variable from formula object


variables = predictors <- c( "OUT_DEAD_DURING_ICU_YN",
  "NEW_SITE_ID", "INC_LOS_PRIOR_ADM_INT", "INC_CARDIAC_DISEASE_YN", "INC_ASTHMA_YN",
  "ICU_RESP_SUPPORT_YN", "ICU_PLATELETS_INT", "ICU_CARDIAC_THERAPY_YN", "ICU_SEPSIS_YN",
  "ICU_STRESS_MYOC_YN", "ICU_PNEUMOTHORAX_YN", "ICU_ATELECTASIS_YN", "ICU_DELIRIUM_YN",
  "ICU_OBSTRUCTION_YN", "ICU_CORTICO_YN", "ICU_SEDAT_DURATION_INT", "ICU_RRT_DIAL_YN",
  "ICU_INOTROPES_YN", "RESP_INTUBATED_YN", "RESP_HFNC_YN",
  "RESP_INV_VENT_YN", "RESP_DURATION_INV_VENT_INT", "RESP_ECMO_YN", "RESP_PRONE_YN",
  "INF_DURING_ICU_YN", "INC_DIABETES1_YN", "neutrophil_lymphocyte_ratio",
  "ICU_CRP_INT", "INC_AGE_INT"
)
# impute the test set
test_data_preimputation = test_data %>% dplyr::select(all_of(variables))
test_data_imputation = futuremice(test_data_preimputation, m = 5, ncore = availableCores()-1, maxit = 10)


    # recompute age group 
modified_list <- lapply(1:test_data_imputation$m, function(i) {
  complete(test_data_imputation, i) %>%
    mutate(
      # Create age groups
      age_group = case_when(
        INC_AGE_INT < 40 ~ "<40",
        INC_AGE_INT >= 40 & INC_AGE_INT < 50 ~ "40-49",
        INC_AGE_INT >= 50 & INC_AGE_INT < 60 ~ "50-59",
        INC_AGE_INT >= 60 & INC_AGE_INT < 70 ~ "60-69",
        INC_AGE_INT >= 70 & INC_AGE_INT < 80 ~ "70-79",
        INC_AGE_INT >= 80 & INC_AGE_INT < 90 ~ "80-89",
        TRUE ~ "90+"
      )) %>%
      # Scale numeric variables and convert back to mids object
      mutate(across(where(is.numeric), ~ as.vector(scale(.)))) 
    })


# Step 2: Convert back to mids object
test_data_imputed <- datlist2mids(modified_list)


# Create CRP categories from train_data_imputation and keep as mids object
# Step 1: Extract each imputed dataset and add CRP categories
modified_list <- lapply(1:test_data_imputed$m, function(i) {
  complete(test_data_imputed, i) %>%
    mutate(
      # Create CRP categories: low, medium, high
      crp_category = case_when(
        ICU_CRP_INT < 100 ~ "low",
        ICU_CRP_INT >= 100 & ICU_CRP_INT < 200 ~ "medium",
        ICU_CRP_INT >= 200 ~ "high",
        TRUE ~ NA_character_  # for any missing values
      )
    )
})
# Step 2: Convert back to mids object
test_data_imputed <- datlist2mids(modified_list)

# Step 2: Convert back to mids object
test_data_imputed <- datlist2mids(modified_list)

test_data_imputed_scaled_complete = complete(test_data_imputed)

#ROC curve 
library(pROC)

# Generate ROC curve
roc_curve <- roc(test_data_imputed_scaled_complete$OUT_DEAD_DURING_ICU_YN, 
                 predict(train_data_imputed_scaled_muiltivariable, newdata = test_data_imputed_scaled_complete, type = "response"))

# Plot ROC curve and save to file
fig_dir <- "figures/model_assessment"
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
roc_file <- file.path(fig_dir, "roc_curve_test.png")
png(filename = roc_file, width = 800, height = 600, res = 150)
plot(roc_curve, col = "blue", lwd = 2, main = sprintf("ROC curve (AUC = %0.3f)", auc(roc_curve)))
dev.off()
cat("Saved ROC plot to:", roc_file, "\n")

cat(sprintf("AUC = %0.4f\n", auc(roc_curve)))

# optimal cut off and scalibration statistic


pred_prob <- predict(train_data_imputed_scaled_muiltivariable,
                     newdata = test_data_imputed_scaled_complete,
                     type = "response")

y_true <- as.integer(test_data_imputed_scaled_complete$OUT_DEAD_DURING_ICU_YN)


youden_thr <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")  # returns threshold

threshold <- as.numeric(youden_thr)
pred_class <- as.integer(pred_prob >= threshold)

TP <- sum(pred_class == 1 & y_true == 1, na.rm = TRUE)
FP <- sum(pred_class == 1 & y_true == 0, na.rm = TRUE)
FN <- sum(pred_class == 0 & y_true == 1, na.rm = TRUE)

precision <- if ((TP + FP) == 0) 0 else TP / (TP + FP)
recall    <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
f1_score  <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)

cat(sprintf("Youden threshold = %0.3f  Precision = %0.3f  Recall = %0.3f  F1 = %0.3f\n",
            threshold, precision, recall, f1_score))

#Brier score
library(DescTools)
BrierScore(train_data_imputed_scaled_muiltivariable) # 0.1551556


# assess model performance 
model_performance(train_data_imputed_scaled_muiltivariable)

#  AIC    |   AICc |    BIC |    R2 |  RMSE | Sigma
#------------------------------------------------
# 5471.1 | 5471.6 | 5716.0 | 0.308 | 0.394 | 0.395

#Takes a while to complete
#check_collinearity(train_data_imputed_scaled_muiltivariable) # low colinearity 
#check_outliers(train_data_imputed_scaled_muiltivariable) # no outliers 
#check_homogeneity(train_data_imputed_scaled_muiltivariable) # normality assumed
# check_predictions(train_data_imputed_scaled_muiltivariable)

print("Assesment of imputed muiltivariable model completed")

### save train_data_imputed_scaled_muiltivariable regression model 

muiltivariable_model = train_data_imputed_scaled_muiltivariable %>% 
  tbl_regression(exp = TRUE )

#save as csv
write.csv(as_tibble(muiltivariable_model), "Figures/model_assessment/muiltivariable_model_table.csv", row.names = FALSE)


