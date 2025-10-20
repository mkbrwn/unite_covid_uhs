#This script is for assessment of the predictive modelling from the UNTIE COVID data set
# The 04_predictive_model.r script produces a multivariable model from backwards stepwise selection after muiltiple imputation. 

# Matt Strammer says ~ forests, calibration and clinical utility. 

# Load required packages
library(caret)
library(mice)
library(future) # parrallelisation of MICE 
library(performance)

# source of training and test set 

source("src/01_raw_to_pre_imputation.r")


#formulat for the final imputed model 

formula <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
             ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
             + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + ICU_TRACHEOS_YN + RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
             + neutrophil_lymphocyte_ratio + ICU_CRP_INT + age_group

train_data_imputed_scaled_muiltivariable = glm(formula, data = complete(train_data_imputation_scaled),binomial(link = "logit"))

# extract variable from formula object


variables = predictors <- c( "OUT_DEAD_DURING_ICU_YN",
  "NEW_SITE_ID", "INC_LOS_PRIOR_ADM_INT", "INC_CARDIAC_DISEASE_YN", "INC_ASTHMA_YN",
  "ICU_RESP_SUPPORT_YN", "ICU_PLATELETS_INT", "ICU_CARDIAC_THERAPY_YN", "ICU_SEPSIS_YN",
  "ICU_STRESS_MYOC_YN", "ICU_PNEUMOTHORAX_YN", "ICU_ATELECTASIS_YN", "ICU_DELIRIUM_YN",
  "ICU_OBSTRUCTION_YN", "ICU_CORTICO_YN", "ICU_SEDAT_DURATION_INT", "ICU_RRT_DIAL_YN",
  "ICU_INOTROPES_YN", "ICU_TRACHEOS_YN", "RESP_INTUBATED_YN", "RESP_HFNC_YN",
  "RESP_INV_VENT_YN", "RESP_DURATION_INV_VENT_INT", "RESP_ECMO_YN", "RESP_PRONE_YN",
  "INF_DURING_ICU_YN", "INC_DIABETES1_YN", "neutrophil_lymphocyte_ratio",
  "ICU_CRP_INT", "INC_AGE_INT"
)

# impute the test set 
test_data_preimputation = test_data %>% dplyr::select(all_of(variables))#


    # recompute age group 
    test_data_preimputation = test_data_preimputation %>%
     mutate(
      age_group = case_when(
        INC_AGE_INT < 40 ~ "<40",
        INC_AGE_INT >= 40 & INC_AGE_INT < 50 ~ "40-49",
        INC_AGE_INT >= 50 & INC_AGE_INT < 60 ~ "50-59",
        INC_AGE_INT >= 60 & INC_AGE_INT < 70 ~ "60-69",
        INC_AGE_INT >= 70 & INC_AGE_INT < 80 ~ "70-79",
        INC_AGE_INT >= 80 & INC_AGE_INT < 90 ~ "80-89",
        TRUE ~ "90+"
      ))
    #drop INC_AGE_INT
        test_data_preimputation = test_data_preimputation %>%
        select(-INC_AGE_INT)  

# impute the test set
test_data_imputed = futuremice(test_data_preimputation, m = 5, core = 8)
test_data_imputed_complete = complete(test_data_imputed)


#ROC curve 
library(pROC)

# Generate ROC curve
roc_curve <- roc(test_data_imputed_complete$OUT_DEAD_DURING_ICU_YN, 
                 predict(train_data_imputed_scaled_muiltivariable, newdata = test_data_imputed_complete, type = "response"))

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
                     newdata = test_data_imputed_complete,
                     type = "response")

y_true <- as.integer(test_data_imputed_complete$OUT_DEAD_DURING_ICU_YN)


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
#model_performance(train_data_imputed_scaled_muiltivariable)

#  AIC    |   AICc |    BIC |    R2 |  RMSE | Sigma
#------------------------------------------------
# 5471.1 | 5471.6 | 5716.0 | 0.308 | 0.394 | 0.395

#Takes a while to complete
#check_collinearity(train_data_imputed_scaled_muiltivariable) # low colinearity 
#check_outliers(train_data_imputed_scaled_muiltivariable) # no outliers 
#check_homogeneity(train_data_imputed_scaled_muiltivariable) # normality assumed
# check_predictions(train_data_imputed_scaled_muiltivariable)

print("Assesment of imputed muiltivariable model completed")

