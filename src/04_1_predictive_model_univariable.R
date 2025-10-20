# GLM Analysis for ICU Mortality - Variable Selection with P < 0.2
# This script performs univariable GLM analysis for all variables in the UNITE_2020_corrected dataset

library(tidyverse)
library(broom)


# Source the data preparation script to get processed data
source("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis/src/01_raw_to_pre_imputation.R")

#test and split data using caret
library(caret)  
set.seed(123)  # For reproducibility
train_indices <- createDataPartition(UNITE_2020_corrected$OUT_DEAD_DURING_ICU_YN, p = 0.7, list = FALSE)
train_data <- UNITE_2020_corrected[train_indices, ]
test_data <- UNITE_2020_corrected[-train_indices, ]

cat("Starting GLM analysis for ICU mortality...\n")

# Create ICU mortality outcome variable
UNITE_2020_corrected <- train_data %>%
  mutate(
    icu_mortality = case_when(
      OUT_DEAD_DURING_ICU_YN == TRUE ~ 1,
      OUT_DEAD_DURING_ICU_YN == FALSE ~ 0,
      TRUE ~ NA_real_
    )
  )

cat("ICU Mortality Distribution:\n")
mortality_table <- table(UNITE_2020_corrected$icu_mortality, useNA = "always")
print(mortality_table)
cat("Proportion of ICU mortality:", 
    round(mean(UNITE_2020_corrected$icu_mortality, na.rm = TRUE), 3), "\n\n")

# Define variables to exclude from analysis
exclude_vars <- c(
  "icu_mortality",
  "OUT_DEAD_DURING_ICU_YN", 
  "OUT_DEATH_DURING_ICU_YN",
  "OUTCOME_LD",
  "OUTCOME_LD_DEATH",
  "OUTCOME_LD_DISCHARGED", 
  "OUTCOME_LD_TRANSFERRED",
  "OUTCOME_LD_OTHER",
  "wave",
  # Exclude identifier variables
  "NEW_CENTRE_ID", "NEW_PATIENT_ID", "NEW_COUNTRY_ID", "NEW_SUBJECT_ID"
)

# Get predictor variables
all_vars <- names(UNITE_2020_corrected)
predictor_vars <- all_vars[!all_vars %in% exclude_vars]

cat("Total predictor variables to analyze:", length(predictor_vars), "\n")

# Function to safely perform GLM
safe_glm <- function(var_name, data) {
  tryCatch({
    # Skip if variable has no variation or too many missing values
    var_data <- data[[var_name]]
    non_missing <- !is.na(var_data) & !is.na(data$icu_mortality)
    
    if (sum(non_missing) < 20) {
      return(data.frame(
        variable = var_name,
        estimate = NA,
        std_error = NA,
        p_value = NA,
        n_obs = sum(non_missing),
        issue = "Insufficient data"
      ))
    }
    
    # Check for variation
    unique_vals <- length(unique(var_data[non_missing]))
    if (unique_vals < 2) {
      return(data.frame(
        variable = var_name,
        estimate = NA,        std_error = NA,
        p_value = NA,
        n_obs = sum(non_missing),
        issue = "No variation"
      ))
    }
    
    # Create analysis dataset
    analysis_data <- data %>%
      select(icu_mortality, all_of(var_name)) %>%
      filter(!is.na(icu_mortality) & !is.na(.data[[var_name]]))
    
    # Fit GLM
    formula_str <- paste("icu_mortality ~", var_name)
    model <- glm(as.formula(formula_str), 
                 data = analysis_data, 
                 family = binomial())
    
    # Extract coefficient for predictor (not intercept)
    model_summary <- tidy(model)
    predictor_coef <- model_summary %>% 
      filter(term != "(Intercept)") %>%
      slice(1)
    
    return(data.frame(
      variable = var_name,
      estimate = predictor_coef$estimate,
      std_error = predictor_coef$std.error,
      p_value = predictor_coef$p.value,
      n_obs = nrow(analysis_data),
      issue = NA
    ))
    
  }, error = function(e) {
    return(data.frame(
      variable = var_name,
      estimate = NA,
      std_error = NA,
      p_value = NA,
      n_obs = NA,
      issue = paste("Error:", e$message)
    ))
  })
}

# Run GLM for each variable
cat("Running univariable GLM analyses...\n")
glm_results <- map_dfr(predictor_vars, ~safe_glm(.x, UNITE_2020_corrected))

cat("Analysis completed!\n\n")

# Summary of results
cat("=== SUMMARY OF RESULTS ===\n")
cat("Total variables analyzed:", nrow(glm_results), "\n")
cat("Variables with valid results:", sum(!is.na(glm_results$p_value)), "\n")
cat("Variables with issues:", sum(!is.na(glm_results$issue)), "\n\n")

# Filter for p < 0.2
significant_vars <- glm_results %>%
  filter(!is.na(p_value) & p_value < 0.2) %>%
  arrange(p_value)

cat("=== VARIABLES WITH P < 0.2 ===\n")
cat("Number of variables with p < 0.2:", nrow(significant_vars), "\n\n")

if (nrow(significant_vars) > 0) {
  # Display results
  results_display <- significant_vars %>%
    select(variable, estimate, std_error, p_value, n_obs) %>%
    mutate(
      estimate = round(estimate, 4),
      std_error = round(std_error, 4),
      p_value = ifelse(p_value < 0.001, "<0.001", round(p_value, 4))
    )
  
  print(results_display)
  
  # Create the variable list
  variables_p_less_than_0_2 <- significant_vars$variable
  
  cat("\n=== VARIABLE LIST FOR MODEL SELECTION ===\n")
  cat("Variables with p < 0.2 (", length(variables_p_less_than_0_2), " variables):\n")
  for (i in 1:length(variables_p_less_than_0_2)) {
    cat(sprintf("%2d. %s\n", i, variables_p_less_than_0_2[i]))
  }
  
  # Save results
  dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
  
  write_csv(glm_results, "data/processed/univariable_glm_results_icu_mortality.csv")
  write_csv(significant_vars, "data/processed/significant_vars_p_less_than_0_2.csv")
  saveRDS(variables_p_less_than_0_2, "data/processed/variables_p_less_than_0_2_list.rds")
  
  cat("\n=== FILES SAVED ===\n")
  cat("✓ data/processed/univariable_glm_results_icu_mortality.csv\n")
  cat("✓ data/processed/significant_vars_p_less_than_0_2.csv\n")
  cat("✓ data/processed/variables_p_less_than_0_2_list.rds\n")
  
  # Also create an R object in global environment
  assign("variables_p_less_than_0_2", variables_p_less_than_0_2, envir = .GlobalEnv)
  
  cat("\n=== R OBJECT CREATED ===\n")
  cat("✓ 'variables_p_less_than_0_2' saved to global environment\n")
  
} else {
  cat("No variables found with p < 0.2\n")
}

# Show any problematic variables
if (sum(!is.na(glm_results$issue)) > 0) {
  cat("\n=== VARIABLES WITH ISSUES ===\n")
  problem_vars <- glm_results %>%
    filter(!is.na(issue)) %>%
    select(variable, issue, n_obs)
  print(problem_vars)
}

cat("\n=== ANALYSIS COMPLETE ===\n")