#Data analysis for the the UNTIE COVID study
# Exploration of missingness

# Run 01_raw_to_pre_imputation.R to process the data
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


################################################### IMPUTATION ###################################################

test_data_imputation = futuremice(test_data, m=5)
print("Imputation completed of the training data set")

# Create age categories and convert back to mids object
# Step 1: Extract each imputed dataset and add age categories
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
      ),
    )
})

# Step 2: Convert back to mids object using datlist2mids
test_data_imputation <- datlist2mids(modified_list)

# Scale numeric variables that don't end with "_YN" and keep as mids object
# Step 1: Extract and scale each imputed dataset
scaled_list <- lapply(1:test_data_imputation$m, function(i) {
  data <- complete(test_data_imputation, i)

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
test_data_imputation_scaled <- datlist2mids(scaled_list)

print("Age groups created and scales for test_data_imputation mids object")

