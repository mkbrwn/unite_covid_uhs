# This script trains a predictive model using the imputed training data
# train_data_imputation is the data taken from src/03_imputation_training_data.
#       2. Plan is to scale variables 
#       3. train predictive model with interaction term 

#run previous script (WARNING TAKE about 5 mins)
# source("src/03_imputation_training_data.r")

#Make a crp categories variable

# Load required library
library(miceadds)

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
  data %>%
    mutate(across(all_of(cols_to_scale), ~ as.vector(scale(.))))
})

# Step 2: Convert back to mids object
train_data_imputation_scaled <- datlist2mids(scaled_list)

# Remove INC_AGE_INT and ICU_CRP_INT from mids object 
# Since we have categorical versions (age_group, crp_category), remove continuous versions to avoid multicollinearity

# Extract datasets, remove variables, and convert back to mids
final_list <- lapply(1:train_data_imputation_scaled$m, function(i) {
  complete(train_data_imputation_scaled, i) %>%
    select(-INC_AGE_INT, -ICU_CRP_INT)
})

# Convert back to mids object
train_data_imputation_scaled <- datlist2mids(final_list)

print("Variables scaled and saved as train_data_imputation_scaled mids object")

# fit muiltivariable model 
train_data_imputed_scaled_muiltivariable = glm(OUT_DEAD_DURING_ICU_YN ~ ., data = complete(train_data_imputation_scaled))


x = step(train_data_imputed_scaled_muiltivariable, direction = "backward", trace = FALSE)
# backwards step elimination across all imputed data FROM train_data_imputation_scaled

# Load required libraries for stepwise regression
library(MASS)


# Method 1: Backwards elimination on each imputed dataset separately
cat("Performing backwards stepwise elimination on each imputed dataset...\n")

# Fit full models on each imputed dataset
full_models <- with(train_data_imputation_scaled, {
  glm(OUT_DEAD_DURING_ICU_YN ~ ., family = binomial, data = data)
})

# Perform backwards elimination on each dataset
backwards_models <- lapply(1:train_data_imputation_scaled$m, function(i) {
  cat("Processing imputation", i, "of", train_data_imputation_scaled$m, "\n")
  
  # Get the complete dataset for this imputation
  complete_data <- complete(train_data_imputation_scaled, i)
  
  # Fit full model
  full_model <- glm(OUT_DEAD_DURING_ICU_YN ~ ., family = binomial, data = complete_data)
  
  # Perform backwards elimination using AIC
  step_model <- step(full_model, direction = "backward", trace = FALSE)
  
  return(step_model)
})

# Extract final variables from each model
final_variables_list <- lapply(backwards_models, function(model) {
  # Get variable names from the model terms (not coefficients)
  # This preserves the original variable names without factor level expansions
  vars <- attr(terms(model), "term.labels")
  return(vars)
})

# Find variables that appear in majority of models (>= 3 out of 5)
all_variables <- unique(unlist(final_variables_list))
variable_counts <- sapply(all_variables, function(var) {
  sum(sapply(final_variables_list, function(vars) var %in% vars))
})

# Variables selected in majority of imputations
consensus_variables <- names(variable_counts)[variable_counts >= 3]

cat("\nBackwards elimination results:\n")
cat("Variables selected in each imputation:\n")
for(i in 1:length(final_variables_list)) {
  cat("Imputation", i, ":", length(final_variables_list[[i]]), "variables\n")
  cat("  ", paste(final_variables_list[[i]], collapse = ", "), "\n")
}

cat("\nConsensus variables (selected in ≥3 imputations):\n")
cat(paste(consensus_variables, collapse = "\n"))
cat("\n")

# Fit final consensus model using MICE pooling
if(length(consensus_variables) > 0) {
  # Create formula with consensus variables
  consensus_formula <- as.formula(paste("OUT_DEAD_DURING_ICU_YN ~", paste(consensus_variables, collapse = " + ")))
  
  cat("\nFinal consensus model formula:\n")
  print(consensus_formula)
  
  # Fit consensus model across all imputations
  consensus_models <- with(train_data_imputation_scaled, {
    glm(consensus_formula, family = binomial, data = complete(train_data_imputation_scaled))
  })
  
  # Pool the results
  pooled_consensus_model <- pool(consensus_models)
  
  cat("\nPooled consensus model results:\n")
  print(summary(pooled_consensus_model))
  
} else {
  cat("No consensus variables found. Consider lowering the threshold.\n")
}

print("Backwards stepwise elimination completed")


formula <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_HBP_YN + INC_ASTHMA_YN +
             INC_KIDNEY_DISEASE_YN + ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN + ICU_KIDNEY_INJ_YN + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT 
             + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + ICU_TRACHEOS_YN + RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + RESP_NEUROM_BLOC_YN + INF_DURING_ICU_YN + ventilation_severity + INC_DIABETES1_YN +
              comorbidity_score + neutrophil_lymphocyte_ratio + ICU_CRP_CATEGORY + age_group