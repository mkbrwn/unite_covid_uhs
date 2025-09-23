# Predictive modeling

library(tidyverse)
library(gtsummary)
library(openxlsx)
library(MASS)

# Load the processed data
source("src/01_raw_to_processed.R")

# univariable has to be done in two chunks due to limited ram on uni laptop  

# 1. Extract first 20 variable names
first20_vars <- names(UNITE_2020_corrected_preimputation)[1:20]

# 2. Run tbl_uvregression with those variables
univariable_ICUmort_1 = tbl_uvregression(
  data       = UNITE_2020_corrected_preimputation,
  method     = glm,                       # or lm, coxph, etc.
  y          = OUT_DEAD_DURING_ICU_YN,     # replace with your dependent var
  include    = all_of(first20_vars),      # injects your first 20 names
  method.args = list(family = binomial)   # only for glm; drop or change for other methods
)

# 1. Extract after 20 variable names
last20_vars <- names(UNITE_2020_corrected_preimputation)[21:ncol(UNITE_2020_corrected_preimputation)]

univariable_ICUmort_2 = tbl_uvregression(
  data       = UNITE_2020_corrected_preimputation,
  method     = glm,                       # or lm, coxph, etc.
  y          = OUT_DEAD_DURING_ICU_YN,     # replace with your dependent var
  include    = all_of(last20_vars),      # injects your first 20 names
  method.args = list(family = binomial)   # only for glm; drop or change for other methods
)

# Extract P<0.2
univariable_ICUmort_1_significant <- univariable_ICUmort_1$table_body %>%
  filter(p.value < 0.2) %>%
  select(variable)
sig_vars <- pull(univariable_ICUmort_1_significant, variable)

univariable_ICUmort_2_significant <- univariable_ICUmort_2$table_body %>%
  filter(p.value < 0.2) %>%
  select(variable)

sig_vars = c(sig_vars, pull(univariable_ICUmort_2_significant, variable))

# 2. Build a formula (replace "OUT_DEAD_DURING_ICU_YN" with your outcome)
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "))
)

#multivariable logistic regression
final_model <- stepAIC(glm(multi_formula, data = UNITE_2020_corrected_preimputation, family = "binomial"), direction = "both")

