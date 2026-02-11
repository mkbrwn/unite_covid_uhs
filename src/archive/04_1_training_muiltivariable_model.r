# muiltivariable model 

# formulat taken from 04_training_predictive_model.r
# data taken from train_data_imputation_scaled which was named in 04_training_predictive_model.r from 03_imputation_training_data.r 


################################################### preopare imputed data set for muiltivariable mode ###################################################

# Load required library
library(miceadds)
library(performance)

# Create CRP categories from train_data_imputation and keep as mids object
# Step 1: Extract each imputed dataset and add CRP categories
modified_list <- lapply(1:train_data_imputation$m, function(i) {
  complete(train_data_imputation, i) %>%
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
train_data_imputation <- datlist2mids(modified_list)

print("CRP categories (low/medium/high) added to train_data_imputation mids object")

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
      ),
    )
})

# Step 2: Convert back to mids object using datlist2mids
library(miceadds)
train_data_imputation <- datlist2mids(modified_list)

# scale variables from train_data_imputation and keep as mids object 

# Scale numeric variables that don't end with "_YN" and keep as mids object
# Step 1: Extract and scale each imputed dataset
scaled_list <- lapply(1:train_data_imputation$m, function(i) {
  data <- complete(train_data_imputation, i)
  
  # Identify numeric columns that don't end with "_YN" and exclude ID/categorical vars
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  cols_to_scale <- numeric_cols[!grepl("_YN$", numeric_cols) & 
                                !grepl("^NEW_SITE_ID$", numeric_cols) &
                                !grepl("^wave$", numeric_cols)]
  
  # Scale only the identified columns and convert to vectors
  #data %>%
  # mutate(across(all_of(cols_to_scale), ~ as.vector(scale(.))))
})

# Step 2: Convert back to mids object
train_data_imputation_scaled <- datlist2mids(scaled_list)




################################################### imputed data ready for muiltivariable model###################################################

formula <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
             ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
             + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + ICU_TRACHEOS_YN + RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
             + neutrophil_lymphocyte_ratio + age_group + ICU_CRP_INT

              # Fit the model
            train_data_imputed_scaled_muiltivariable = glm(formula, data = complete(train_data_imputation_scaled))
            summary(train_data_imputed_scaled_muiltivariable)
            tbl_regression(train_data_imputed_scaled_muiltivariable, exponentiate = TRUE) %>% print()
            
            #save table tbl_regression as csv
            muiltivariable_model_table = tbl_regression(train_data_imputed_scaled_muiltivariable, exponentiate = TRUE)
            muiltivariable_model_table %>%
              as_tibble() %>%
              write.csv("Figures/model_assessment/muiltivariable_model_table.csv", row.names = FALSE)

print("Imputed muiltivariable model complete")

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
