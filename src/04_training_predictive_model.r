# This script trains a predictive model using the imputed training data
# train_data_imputation is the data taken from src/03_imputation_training_data.
#       2. Plan is to scale variables
#       3. train predictive model with interaction term

#run previous script (WARNING TAKE about 5 mins)
#source("src/03_imputation_training_data.r")
source("src/01_raw_to_pre_imputation.R")

#library
library(gtsummary)
library(mice) #muiltiple imputation
library(miceadds)
library(future.apply) # parrallelisation of MICE
library(naniar) #visualisation of missingness
library(caret) # For colinearity detection
library(MASS) # backwards step AIC
library(caret)


#
train_data_imputation = futuremice(train_data, m=5)
print("Imputation completed of the training data set")

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
})# Step 2: Convert back to mids object using datlist2mids
library(miceadds)
train_data_imputation <- datlist2mids(modified_list)

print("Age groups created and integrated into train_data_imputation mids object")


#install gtsummary4mice if not already installed
if (!requireNamespace("gtsummary4mice", quietly = TRUE)) {
  remotes::install_github(
    "jrob95/gtsummary4mice",
    build_vignettes = FALSE,
    dependencies    = TRUE
  )
}

# Load required library
library(miceadds)
library(gtsummary4mice) # for tbl_uvregression for mids object 

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
})# Step 2: Convert back to mids object
train_data_imputation_scaled <- datlist2mids(scaled_list)

# Remove INC_AGE_INT and ICU_CRP_INT from mids object 
# Since we have categorical versions (age_group, crp_category), remove continuous versions to avoid multicollinearity

# Extract datasets, remove variables, and convert back to mids
final_list <- lapply(1:train_data_imputation_scaled$m, function(i) {
  complete(train_data_imputation_scaled, i) %>%
    dplyr::select(-INC_AGE_INT, -NEW_BMI)
})

# Convert back to mids object
train_data_imputation_scaled <- datlist2mids(final_list)

print("Variables scaled and saved as train_data_imputation_scaled mids object")

# univariable analysis with P<0.2 
tbl_uv <- tbl_uvregression(
  train_data_imputation_scaled,
  method = glm,
  y = OUT_DEAD_DURING_ICU_YN,
  exp = TRUE
)

# extract varaibles with p<0.2
variables <- tbl_uv$table_body %>%
  filter(p.value < 0.2) %>%
  pull(variable)

#isolate further variable not from admission

variables_to_remove = c(
  #Complications during ICYU stay
  "ICU_PRESSURE_OTH_YN",
  "ICU_KIDNEY_INJ_YN",
  "ICU_OBSTRUCTION_YN",
  "ICU_EXTUB_YN",
  "ICU_ANTIVIRALS_YN",
  "ICU_CARDIAC_THERAPY_YN",
  "ICU_SEPSIS_YN",
  "ICU_STRESS_MYOC_YN",
  "ICU_MYOCARDITIS_YN",
  "ICU_PERICARD_YN",
  "ICU_PNEUMOTHORAX_YN",
  "ICU_ATELECTASIS_YN",
  "ICU_DELIRIUM_YN",
  "ICU_SEIZURE_YN",
  "ICU_PRESSURE_FAC_YN",

  #treatment during ICU stay
  "ICU_INOTROPES_YN",
  "ICU_ANTIMALARIAL_YN",
  "ICU_ANTIVIRALS_YN",
  "ICU_SEDATION_YN",
  "ICU_SEDAT_DURATION_INT",
  "ICU_RRT_DIAL_YN",
  "RESP_INTUBATED_ICU_STAY_YN",
  "RESP_NI_VENT_YN",
  "RESP_HFNC_YN",
  "RESP_DURATION_INV_VENT_INT",
  "RESP_ECMO_YN",
  "RESP_PRONE_YN",
  "RESP_NEUROM_BLOC_YN",
  "RESP_MODE_RAD",
  "INF_ANTIBIO_YN",
  "INF_ANTIFUNG_YN",
  "INF_DURING_ICU_YN",
  "INC_BMI_INT",
  "ventilation_severity",
  "ICU_BLOOD_PURIF_YN",

  #duplicates
  "ICU_CRP_CATEGORY",
  "crp_category",
  "comorbidity_score"

)

variables <- setdiff(variables, variables_to_remove)

# Create restricted imputed data excluding variables_to_remove that actually exist
train_data_imputation_restricted <- lapply(1:train_data_imputation_scaled$m, function(i) {
  dat <- complete(train_data_imputation_scaled, i)
  # Only remove variables that actually exist in the dataset
  vars_to_drop <- intersect(variables_to_remove, names(dat))
  dat %>% dplyr::select(-all_of(vars_to_drop))
}) %>%
  datlist2mids()

# Run univariable regression on restricted dataset, restrict variables which have p<0.2
tbl_uv_restricted <- tbl_uvregression(
  train_data_imputation_restricted,
  method = glm,
  y = OUT_DEAD_DURING_ICU_YN,
  exp = TRUE
)

openxlsx::write.xlsx(tbl_uv_restricted %>% as_tibble(), "Figures/model_assessment/univariable_regression_results.xlsx", rowNames = FALSE)
print("Univariable regression analysis complete.")

#formula for multivariable model using variables with p<0.2 from univariable analysis
formula = as.formula(
  paste("OUT_DEAD_DURING_ICU_YN ~", paste(variables, collapse = " + "))
)

# backwards step elimination for multivariable model across imuted data
cat("Performing backwards stepwise elimination on each imputed dataset...\n")

# Perform backwards elimination on each dataset
backwards_models <- lapply(1:train_data_imputation_scaled$m, function(i) {
  cat("Processing imputation", i, "of", train_data_imputation_scaled$m, "\n")

  # Get the complete dataset for this imputation
  complete_data <- complete(train_data_imputation_scaled, i)

  # Fit full model
  full_model <- glm(formula, family = binomial, data = complete_data)

  # Perform backwards elimination using AIC
  step_model <- step(full_model, direction = "backward", trace = FALSE)

  return(step_model)
})# Extract final variables from each model
final_variables_list <- lapply(backwards_models, function(model) {
  # Get variable names from the model terms (not coefficients)
  # This preserves the original variable names without factor level expansions
  vars <- attr(terms(model), "term.labels")
  return(vars)
})

# resolve final variables by majority vote across imputations
final_variable_counts <- table(unlist(final_variables_list))
final_variables <- names(final_variable_counts[final_variable_counts > (length(backwards_models) / 2)]) 

# In this format for explicit expression of formula
vars <- c(
  "age_group",
  "ICU_CORTICO_YN",
  "ICU_CRP_INT",
  "ICU_NEUTRO_INT",
  "ICU_PLATELETS_INT",
  "ICU_RESP_SUPPORT_YN",
  "ICU_WHITE_CELL_INT",
  "INC_CARDIAC_DISEASE_YN",
  "INC_DIABETES1_YN",
  "INC_IMMUNOSUPPR_YN",
  "INC_KIDNEY_DISEASE_YN",
  "INC_LIVER_DISEASE_YN",
  "INC_LOS_PRIOR_ADM_INT",
  "INC_NEURO_YN",
  "INC_SEX_RAD",
  "INF_AT_ADMISSION_YN",
  "neutrophil_lymphocyte_ratio",
  "RESP_INTUBATED_YN",
  "RESP_INV_VENT_YN",
  "RESP_VENT_ROUTINE_YN"
)

model_formula <- as.formula(
  paste("OUT_DEAD_DURING_ICU_YN ~", paste(vars, collapse = " + "))
)
#final muiltivariable model across imputed data
print(model_formula)