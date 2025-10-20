# This script is for multilevel imputation and data cleaning for publication
# First data cleaning and multiple imputation must occur. 
source("src/03_imputation_publication.r")

# Load additional libraries for multilevel imputation
library(mice)
library(miceadds)
library(lme4)        # for multilevel modeling
library(mitml)       # for multilevel multiple imputation
library(pan)         # for multilevel imputation methods
library(furrr)       # for parallel processing with futuremice
library(future)      # for parallel backend

#remove NEW_SITE_ID with observations less than 10
UNITE_2020_corrected_preimputation <- UNITE_2020_corrected_preimputation %>%
  group_by(NEW_SITE_ID) %>%
  filter(n() >= 10) %>%
  ungroup()

# Specify multilevel imputation with NEW_SITE_ID as random intercept
# Convert NEW_SITE_ID to integer for proper clustering (required by mice)
UNITE_2020_corrected_preimputation$NEW_SITE_ID <- as.integer(as.factor(UNITE_2020_corrected_preimputation$NEW_SITE_ID))

# Set up imputation method for multilevel data
# Initialize mice object to get predictor matrix and method
init_mice <- mice(UNITE_2020_corrected_preimputation, maxit = 0)
pred_matrix <- init_mice$predictorMatrix
method <- init_mice$method

# Modify predictor matrix to specify multilevel structure
# Set NEW_SITE_ID as cluster variable (value = -2)
pred_matrix[, "NEW_SITE_ID"] <- -2

# For continuous variables, use multilevel normal imputation (2l.norm)
# For categorical variables, use multilevel logistic/multinomial imputation
continuous_vars <- names(UNITE_2020_corrected_preimputation)[sapply(UNITE_2020_corrected_preimputation, is.numeric)]
categorical_vars <- names(UNITE_2020_corrected_preimputation)[sapply(UNITE_2020_corrected_preimputation, function(x) is.factor(x) | is.character(x))]

# Remove NEW_SITE_ID from variables to be imputed (it's the cluster variable)
continuous_vars <- continuous_vars[continuous_vars != "NEW_SITE_ID"]
categorical_vars <- categorical_vars[categorical_vars != "NEW_SITE_ID"]

# Set imputation methods for multilevel structure
for(var in continuous_vars) {
  if(sum(is.na(UNITE_2020_corrected_preimputation[[var]])) > 0) {
    method[var] <- "2l.norm"  # Multilevel normal imputation
  }
}

for(var in categorical_vars) {
  if(sum(is.na(UNITE_2020_corrected_preimputation[[var]])) > 0) {
    # Check if binary or multinomial
    n_levels <- length(unique(na.omit(UNITE_2020_corrected_preimputation[[var]])))
    if(n_levels == 2) {
      method[var] <- "2l.bin"    # Multilevel binary logistic
    } else {
      method[var] <- "2l.pmm"    # Multilevel predictive mean matching for multinomial
    }
  }
}

# Perform multilevel multiple imputation using futuremice for parallel processing
UNITE_2020_corrected_imputation <- futuremice(
  UNITE_2020_corrected_preimputation, 
  method = method,
  predictorMatrix = pred_matrix,
  m = 5,                    # Number of imputations    
)

print("Multilevel multiple imputation with NEW_SITE_ID as random intercept completed")

# Diagnostic checks for multilevel imputation
print("Generating diagnostic plots for multilevel imputation...")

# Check convergence of the imputation algorithm
pdf("figures/convergence_multilevel_imputation.pdf", width = 12, height = 8)
plot(UNITE_2020_corrected_imputation, layout = c(2, 2))
dev.off()

# Density plots to compare observed and imputed values by site
png("figures/densityplot_multilevel_imputation.png", width = 12, height = 8, units = "in", res = 300, bg = "white")
densityplot(UNITE_2020_corrected_imputation)
dev.off()

# Strip plots to show distribution of imputed values across sites
png("figures/stripplot_multilevel_imputation.png", width = 12, height = 8, units = "in", res = 300, bg = "white")
stripplot(UNITE_2020_corrected_imputation, pch = 20, cex = 1.2)
dev.off()

# Check the between-site and within-site variation
# Extract completed data for analysis
completed_data_list <- complete(UNITE_2020_corrected_imputation, action = "all")

# Summary of imputation by site
print("Summary of observations by site:")
print(table(UNITE_2020_corrected_preimputation$NEW_SITE_ID, useNA = "ifany"))

# Check if imputation preserves site-level relationships
print("Checking preservation of multilevel structure...")

# Example: Compare means of a key variable (e.g., first continuous variable) by site
# before and after imputation
if(length(continuous_vars) > 0) {
  first_cont_var <- continuous_vars[1]
  
  # Original data means by site (excluding missing)
  original_means <- aggregate(
    UNITE_2020_corrected_preimputation[[first_cont_var]], 
    by = list(Site = UNITE_2020_corrected_preimputation$NEW_SITE_ID), 
    FUN = mean, na.rm = TRUE
  )
  
  # Imputed data means by site (average across imputations)
  imputed_means_list <- lapply(completed_data_list, function(data) {
    aggregate(
      data[[first_cont_var]], 
      by = list(Site = data$NEW_SITE_ID), 
      FUN = mean, na.rm = TRUE
    )
  })
  
  # Average across imputations
  imputed_means <- Reduce(function(x, y) {
    merge(x, y, by = "Site", all = TRUE)
  }, imputed_means_list)
  
  # Calculate mean across all imputations
  imputed_means$mean_across_imputations <- rowMeans(imputed_means[, -1], na.rm = TRUE)
  imputed_means <- imputed_means[, c("Site", "mean_across_imputations")]
  
  # Compare
  comparison <- merge(original_means, imputed_means, by = "Site", all = TRUE)
  names(comparison) <- c("Site", "Original_Mean", "Imputed_Mean")
  print("Comparison of site means before and after imputation:")
  print(comparison)
}

print("Multilevel imputation diagnostics completed")

# Save the multilevel imputation object
saveRDS(UNITE_2020_corrected_imputation, "data/processed/multilevel_imputation_object.rds")
print("Multilevel imputation object saved to data/processed/multilevel_imputation_object.rds")


