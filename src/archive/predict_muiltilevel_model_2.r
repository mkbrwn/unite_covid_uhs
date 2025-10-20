# Fit multilevel model across all imputed datasets using with() from micest data was cleaned with src/03_imputation_muiltilevel_publication.r
# This process takes a while and was saved as a RDS file 

#additional packages for this script
library(tidyverse)
library(gtsummary)
library(gtsummary4mice)
library(lme4)
library(broom.mixed)  # Required for tidy() methods with mixed-effects models
library(miceadds)     # For pooling multilevel models

#load RDS file
UNITE_2020_corrected_multiple_imputation <- readRDS("data/processed/multilevel_imputation_object.rds")

# Store the modified complete datasets in a list
modified_datasets <- list()

#categorise age into 10 year groups from 40 years and scale continuous variables
for (i in 1:5) {
    # Extract the complete dataset
    complete_data <- mice::complete(UNITE_2020_corrected_multiple_imputation, action = i)
    
    # Check if INC_AGE_INT exists
    if ("INC_AGE_INT" %in% names(complete_data)) {
        # Create age groups
        complete_data <- complete_data %>%
            mutate(AGE_GROUP = case_when(
                INC_AGE_INT < 40 ~ "<40",
                INC_AGE_INT >= 40 & INC_AGE_INT < 50 ~ "40-49",
                INC_AGE_INT >= 50 & INC_AGE_INT < 60 ~ "50-59",
                INC_AGE_INT >= 60 & INC_AGE_INT < 70 ~ "60-69",
                INC_AGE_INT >= 70 & INC_AGE_INT < 80 ~ "70-79",
                INC_AGE_INT >= 80 ~ "80+",
                TRUE ~ NA_character_
            )) %>%
            mutate(AGE_GROUP = factor(AGE_GROUP, levels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+"))) %>%
            select(-INC_AGE_INT)  # Remove original INC_AGE_INT variable
    }
    
    # Scale and center continuous variables
    complete_data <- complete_data %>%
        mutate(across(where(is.numeric), ~ scale(.) %>% as.vector()))
    
    # Store the modified dataset
    modified_datasets[[i]] <- complete_data
    
    cat("Processed dataset", i, "\n")
}

# Convert the list of modified datasets back to a mids object
# We need to include the original data (.imp = 0) for mice::as.mids() to work
# Extract the original data and apply the same transformations
original_data <- mice::complete(UNITE_2020_corrected_multiple_imputation, action = 0)

# Apply same transformations to original data
if ("INC_AGE_INT" %in% names(original_data)) {
    original_data <- original_data %>%
        mutate(AGE_GROUP = case_when(
            INC_AGE_INT < 40 ~ "<40",
            INC_AGE_INT >= 40 & INC_AGE_INT < 50 ~ "40-49",
            INC_AGE_INT >= 50 & INC_AGE_INT < 60 ~ "50-59",
            INC_AGE_INT >= 60 & INC_AGE_INT < 70 ~ "60-69",
            INC_AGE_INT >= 70 & INC_AGE_INT < 80 ~ "70-79",
            INC_AGE_INT >= 80 ~ "80+",
            TRUE ~ NA_character_
        )) %>%
        mutate(AGE_GROUP = factor(AGE_GROUP, levels = c("<40", "40-49", "50-59", "60-69", "70-79", "80+"))) %>%
        select(-INC_AGE_INT)
}

# Scale continuous variables in original data
original_data <- original_data %>%
    mutate(across(where(is.numeric), ~ scale(.) %>% as.vector()))

# Stack all datasets with original data (.imp = 0) included
all_data <- bind_rows(
    original_data %>% mutate(.imp = 0, .id = row_number()),
    bind_rows(modified_datasets, .id = ".imp") %>% mutate(.imp = as.integer(.imp), .id = row_number())
)

# Create a mids object from the long format
UNITE_2020_corrected_multiple_imputation_modified <- mice::as.mids(all_data)

cat("Successfully created modified mids object.\n")

cat("\nData transformation complete.\n")
cat("First few rows of dataset 1:\n")
head(modified_datasets[[1]]) 

###### procude 


#univariable analysis - This will take 5 mins 
univariable_multilevel_ICU_mort_modified =tbl_uvregression( ,
        method = glmer,
        exponentiate = TRUE,
        method.args = list(family = binomial),
        y = "OUT_DEAD_DURING_ICU_YN",
        formula = "{y} ~ {x}+ (1 | NEW_SITE_ID)", #specifies the random effects component of the mixed effects mode e.g = each site
        )


# Extract P<0.2
univariable_multilevel_ICU_mort_modified <- univariable_multilevel_ICU_mort_modified$table_body %>%
  filter(p.value < 0.2) %>%
  dplyr::select(variable)
sig_vars <- dplyr::pull(univariable_multilevel_ICU_mort_modified, variable)

# construct full multivariable multilevel model formula
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "), " + (1 | NEW_SITE_ID)")
)

print("Initial formula for backwards elimination:")
print(multi_formula)

# ========================================
# LASSO Feature Selection for GLMER
# ========================================

# Install and load required packages for LASSO with mixed models
if (!require(glmmLasso, quietly = TRUE)) {
  install.packages("glmmLasso")
}
library(glmmLasso)

# Alternative: using glmnet with dummy coding for random effects
if (!require(glmnet, quietly = TRUE)) {
  install.packages("glmnet")
}
library(glmnet)

cat("\n========================================\n")
cat("LASSO FEATURE SELECTION\n")
cat("========================================\n")

# Store LASSO results for each imputed dataset
lasso_results <- list()

for (i in 1:5) {
  cat("\n--- Processing imputed dataset:", i, "---\n")
  
  # Extract the i-th complete dataset
  current_data <- modified_datasets[[i]]
  
  # Prepare data: remove outcome and ID variables
  outcome_var <- "OUT_DEAD_DURING_ICU_YN"
  site_var <- "NEW_SITE_ID"
  
  # Check if variables exist
  if (!outcome_var %in% names(current_data)) {
    cat("Warning: Outcome variable", outcome_var, "not found in dataset", i, "\n")
    next
  }
  
  if (!site_var %in% names(current_data)) {
    cat("Warning: Site variable", site_var, "not found in dataset", i, "\n")
    next
  }
  
  # Remove rows with missing outcome
  current_data_clean <- current_data %>%
    filter(!is.na(!!sym(outcome_var)))
  
  # Get all potential predictor variables (from sig_vars if available)
  if (exists("sig_vars") && length(sig_vars) > 0) {
    predictor_vars <- sig_vars
  } else {
    # Use all variables except outcome and site
    predictor_vars <- setdiff(names(current_data_clean), c(outcome_var, site_var))
  }
  
  # Filter to only include variables that exist in the dataset
  predictor_vars <- predictor_vars[predictor_vars %in% names(current_data_clean)]
  
  cat("Number of predictor variables:", length(predictor_vars), "\n")
  
  # Method 1: Using glmmLasso package
  cat("\nMethod 1: glmmLasso\n")
  
  tryCatch({
    # Prepare data for glmmLasso
    # glmmLasso requires a specific format
    lasso_data <- current_data_clean %>%
      select(all_of(c(outcome_var, site_var, predictor_vars)))
    
    # Convert outcome to numeric (0/1) if it's not already
    if (is.factor(lasso_data[[outcome_var]]) || is.character(lasso_data[[outcome_var]])) {
      lasso_data[[outcome_var]] <- as.numeric(as.character(lasso_data[[outcome_var]]))
    }
    
    # Try a range of lambda values
    lambda_seq <- seq(0, 50, by = 5)
    bic_values <- numeric(length(lambda_seq))
    models_list <- list()
    
    for (j in seq_along(lambda_seq)) {
      lambda_val <- lambda_seq[j]
      
      temp_model <- tryCatch({
        glmmLasso(
          fix = as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + "))),
          rnd = list(as.formula(paste("~1|", site_var))),
          data = lasso_data,
          lambda = lambda_val,
          family = binomial(link = "logit"),
          control = list(print.iter = FALSE)
        )
      }, error = function(e) NULL)
      
      if (!is.null(temp_model)) {
        models_list[[j]] <- temp_model
        bic_values[j] <- temp_model$bic
      } else {
        bic_values[j] <- NA
      }
    }
    
    # Select model with minimum BIC
    if (any(!is.na(bic_values))) {
      best_idx <- which.min(bic_values)
      best_lambda <- lambda_seq[best_idx]
      best_model <- models_list[[best_idx]]
      
      # Extract non-zero coefficients
      coefs <- best_model$coefficients
      selected_vars <- names(coefs[coefs != 0 & names(coefs) != "(Intercept)"])
      
      cat("Best lambda:", best_lambda, "\n")
      cat("Number of selected variables:", length(selected_vars), "\n")
      cat("Selected variables:", paste(selected_vars, collapse = ", "), "\n")
      
      lasso_results[[i]] <- list(
        dataset = i,
        method = "glmmLasso",
        lambda = best_lambda,
        selected_vars = selected_vars,
        coefficients = coefs[coefs != 0],
        bic = bic_values[best_idx],
        model = best_model
      )
    } else {
      cat("glmmLasso failed for all lambda values\n")
    }
    
  }, error = function(e) {
    cat("Error in glmmLasso:", e$message, "\n")
  })
  
  # Method 2: Using glmnet with site as a fixed effect (alternative approach)
  cat("\nMethod 2: glmnet with site dummy variables\n")
  
  tryCatch({
    # Create model matrix
    # Include site as dummy variables
    formula_glmnet <- as.formula(paste("~", paste(c(predictor_vars, site_var), collapse = " + ")))
    x_matrix <- model.matrix(formula_glmnet, data = current_data_clean)[, -1]  # Remove intercept
    y_vector <- current_data_clean[[outcome_var]]
    
    # Convert outcome to numeric if needed
    if (is.factor(y_vector) || is.character(y_vector)) {
      y_vector <- as.numeric(as.character(y_vector))
    }
    
    # Fit LASSO with cross-validation
    cv_lasso <- cv.glmnet(
      x = x_matrix,
      y = y_vector,
      family = "binomial",
      alpha = 1,  # LASSO penalty
      nfolds = 10
    )
    
    # Extract coefficients at lambda.min
    coefs_min <- coef(cv_lasso, s = "lambda.min")
    selected_idx <- which(coefs_min != 0)
    selected_names <- rownames(coefs_min)[selected_idx]
    selected_names <- selected_names[selected_names != "(Intercept)"]
    
    # Remove site dummy variable names and keep only predictor variables
    selected_predictors <- selected_names[!grepl(paste0("^", site_var), selected_names)]
    
    cat("Lambda min:", cv_lasso$lambda.min, "\n")
    cat("Lambda 1se:", cv_lasso$lambda.1se, "\n")
    cat("Number of selected features:", length(selected_names), "\n")
    cat("Selected predictors (excluding site dummies):", length(selected_predictors), "\n")
    cat("Selected variables:", paste(selected_predictors, collapse = ", "), "\n")
    
    # Store results
    if (is.null(lasso_results[[i]])) {
      lasso_results[[i]] <- list(
        dataset = i,
        method = "glmnet",
        lambda_min = cv_lasso$lambda.min,
        lambda_1se = cv_lasso$lambda.1se,
        selected_vars = selected_predictors,
        all_selected = selected_names,
        coefficients = as.vector(coefs_min[selected_idx]),
        cv_model = cv_lasso
      )
    } else {
      # Add glmnet results to existing list
      lasso_results[[i]]$glmnet_method <- list(
        lambda_min = cv_lasso$lambda.min,
        lambda_1se = cv_lasso$lambda.1se,
        selected_vars = selected_predictors,
        all_selected = selected_names,
        cv_model = cv_lasso
      )
    }
    
  }, error = function(e) {
    cat("Error in glmnet:", e$message, "\n")
  })
  
  cat("\nDataset", i, "LASSO selection complete.\n")
}

# ========================================
# Summary of LASSO Results Across All Datasets
# ========================================

cat("\n========================================\n")
cat("SUMMARY OF LASSO FEATURE SELECTION\n")
cat("========================================\n")

# Combine selected variables across all datasets
all_selected_vars <- list()

for (i in seq_along(lasso_results)) {
  if (!is.null(lasso_results[[i]])) {
    cat("\nDataset", i, ":\n")
    
    if (!is.null(lasso_results[[i]]$selected_vars)) {
      cat("  Method:", lasso_results[[i]]$method, "\n")
      cat("  Number of selected variables:", length(lasso_results[[i]]$selected_vars), "\n")
      cat("  Variables:", paste(lasso_results[[i]]$selected_vars, collapse = ", "), "\n")
      all_selected_vars[[i]] <- lasso_results[[i]]$selected_vars
    }
    
    if (!is.null(lasso_results[[i]]$glmnet_method)) {
      cat("\n  Alternative (glmnet):\n")
      cat("  Number of selected variables:", length(lasso_results[[i]]$glmnet_method$selected_vars), "\n")
      cat("  Variables:", paste(lasso_results[[i]]$glmnet_method$selected_vars, collapse = ", "), "\n")
    }
  }
}

# Find variables selected in at least 3 out of 5 datasets (consensus)
if (length(all_selected_vars) > 0) {
  all_vars_flat <- unlist(all_selected_vars)
  var_frequency <- table(all_vars_flat)
  consensus_vars <- names(var_frequency[var_frequency >= 3])
  
  cat("\n========================================\n")
  cat("CONSENSUS VARIABLES (selected in ≥3/5 datasets):\n")
  cat("========================================\n")
  cat(paste(consensus_vars, collapse = "\n"), "\n")
  cat("\nNumber of consensus variables:", length(consensus_vars), "\n")
  
  # Create formula with consensus variables
  if (length(consensus_vars) > 0) {
    lasso_formula <- as.formula(
      paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(consensus_vars, collapse = " + "), " + (1 | NEW_SITE_ID)")
    )
    cat("\nRecommended LASSO-selected formula:\n")
    print(lasso_formula)
  }
}

cat("\n========================================\n")
cat("LASSO feature selection complete.\n")
cat("========================================\n")


