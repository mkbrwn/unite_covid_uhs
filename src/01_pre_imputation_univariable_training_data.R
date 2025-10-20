# Processes raw data files and filters patients who would then be eligable for analsis.

# Load required libraries
library(tidyverse)
library(readr)
library(gtsummary)
library(buildmer)

# Ensure gtsummary4mice package is not loaded (to avoid conflicts)
#detach("package:gtsummary4mice", character.only = TRUE)

#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

source("src/01_raw_to_pre_imputation.R")

univariable_training_output <- train_data %>%
  tbl_uvregression(
    method = glm,
    y = OUT_DEAD_DURING_ICU_YN,
    methods.args = list(family = binomial),
    exponentiate = TRUE
  )
  
# Extract P<0.2
univariable_training_output <- univariable_training_output$table_body %>%
  filter(p.value < 0.2) %>%
  dplyr::select(variable)

#find the significant variables 
sig_vars <- dplyr::pull(univariable_training_output, variable)

# construct full multivariable multilevel model formula
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "))
)

# build the muiltivariable model
multi_model <- glm(multi_formula, data = train_data, family = binomial)

  ##### model building with buildmer #############
print("Building multivariable model starting")

# Ensure subset keeps the outcome and grouping variable 
train_data <- train_data %>%
  dplyr::select(all_of(sig_vars), OUT_DEAD_DURING_ICU_YN, NEW_SITE_ID)

#scale all continuous variables with the dataset
train_data <- train_data %>%
  mutate(across(where(is.numeric), scale))  

# Fit multilevel logistic regression with bobyqa optimizer and increased maxfun 
convergence_model <- glm(OUT_DEAD_DURING_ICU_YN ~ ., data = train_data, family = binomial)

# extract variables from convergence model add OUT_DEAD_DURING_ICU_YN
actual_vars <- attr(terms(convergence_model), "term.labels")

# Create training dataset with only the variables that were actually used in the model
train_data_full <- train_data %>%
  dplyr::select(all_of(actual_vars), OUT_DEAD_DURING_ICU_YN)

multivariable_train_data = step(glm(OUT_DEAD_DURING_ICU_YN ~ ., data = na.omit(train_data_full), family = binomial), )


tbl_regression(multivariable_train_data,
                exponentiate = TRUE
              )

column_names <- colnames(train_data_full)

