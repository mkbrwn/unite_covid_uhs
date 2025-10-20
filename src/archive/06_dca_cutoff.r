

train_data_imputation_complete$CRP_bin = ifelse(train_data_imputation_complete$ICU_CRP_INT < 150, 0,1 )

# simple model
simple_imputed_crp_model = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*CRP_bin, data = train_data_imputation_complete, family = binomial)
summary(simple_imputed_crp_model)

# Test multiple CRP thresholds to find where interaction becomes significant
cat("\nTesting multiple CRP thresholds for significant interaction...\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Define range of thresholds to test
thresholds <- seq(50, 300, by = 25)  # Test thresholds from 50 to 300 mg/L

# Initialize results storage
threshold_results <- data.frame(
  threshold = numeric(),
  interaction_p_value = numeric(),
  interaction_coef = numeric(),
  interaction_se = numeric(),
  steroid_main_effect = numeric(),
  steroid_main_p = numeric(),
  steroid_effect_low_crp = numeric(),
  steroid_effect_high_crp = numeric(),
  beneficial_low_crp = logical(),
  beneficial_high_crp = logical(),
  n_low_crp = numeric(),
  n_high_crp = numeric()
)

# Test each threshold
for(threshold in thresholds) {
  # Create binary CRP variable based on current threshold
  train_data_imputation_complete$CRP_bin_temp <- ifelse(train_data_imputation_complete$ICU_CRP_INT < threshold, 0, 1)
  
  # Skip if one group is too small (less than 10 observations)
  crp_table <- table(train_data_imputation_complete$CRP_bin_temp)
  if(min(crp_table) < 10) next
  
  # Fit model with interaction
  model_temp <- glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN * CRP_bin_temp, 
                   data = train_data_imputation_complete, family = binomial)
  
  # Extract interaction term results
  model_summary <- summary(model_temp)
  coef_matrix <- model_summary$coefficients
  
  # Extract all relevant coefficients
  steroid_main_effect <- coef_matrix["ICU_CORTICO_YNTRUE", "Estimate"]
  steroid_main_p <- coef_matrix["ICU_CORTICO_YNTRUE", "Pr(>|z|)"]
  
  # Find interaction term
  interaction_row <- grep("ICU_CORTICO_YNTRUE:CRP_bin_temp", rownames(coef_matrix))
  
  if(length(interaction_row) > 0) {
    interaction_coef <- coef_matrix[interaction_row, "Estimate"]
    interaction_p <- coef_matrix[interaction_row, "Pr(>|z|)"]
    interaction_se <- coef_matrix[interaction_row, "Std. Error"]
    
    # Calculate steroid effect in low CRP group (main effect only)
    steroid_effect_low_crp <- steroid_main_effect
    
    # Calculate steroid effect in high CRP group (main effect + interaction)
    steroid_effect_high_crp <- steroid_main_effect + interaction_coef
    
    # Determine if steroid is beneficial (negative coefficient = reduced mortality)
    beneficial_low_crp <- steroid_effect_low_crp < 0 & steroid_main_p < 0.05
    beneficial_high_crp <- steroid_effect_high_crp < 0
    
    # Store results
    threshold_results <- rbind(threshold_results, data.frame(
      threshold = threshold,
      interaction_p_value = interaction_p,
      interaction_coef = interaction_coef,
      interaction_se = interaction_se,
      steroid_main_effect = steroid_main_effect,
      steroid_main_p = steroid_main_p,
      steroid_effect_low_crp = steroid_effect_low_crp,
      steroid_effect_high_crp = steroid_effect_high_crp,
      beneficial_low_crp = beneficial_low_crp,
      beneficial_high_crp = beneficial_high_crp,
      n_low_crp = crp_table[1],
      n_high_crp = crp_table[2]
    ))
  }
}

# Display results
cat("Threshold analysis results:\n")
print(threshold_results)

# Find thresholds where glucocorticoids are beneficial
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("GLUCOCORTICOID BENEFIT ANALYSIS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Beneficial in low CRP group
beneficial_low <- threshold_results[threshold_results$beneficial_low_crp == TRUE, ]
if(nrow(beneficial_low) > 0) {
  cat("\nThresholds where steroids benefit LOW CRP patients:\n")
  cat("(CRP below threshold, significant main effect, OR < 1)\n")
  print(beneficial_low[, c("threshold", "steroid_effect_low_crp", "steroid_main_p", "n_low_crp")])
}

# Beneficial in high CRP group
beneficial_high <- threshold_results[threshold_results$beneficial_high_crp == TRUE, ]
if(nrow(beneficial_high) > 0) {
  cat("\nThresholds where steroids benefit HIGH CRP patients:\n")
  cat("(CRP above threshold, combined effect OR < 1)\n")
  print(beneficial_high[, c("threshold", "steroid_effect_high_crp", "interaction_p_value", "n_high_crp")])
}

# Find the best thresholds for each group
if(nrow(beneficial_low) > 0) {
  best_low_threshold <- beneficial_low[which.min(beneficial_low$steroid_main_p), ]
  cat("\nBest threshold for LOW CRP benefit:\n")
  cat("Threshold:", best_low_threshold$threshold, "mg/L\n")
  cat("Steroid effect (log OR):", round(best_low_threshold$steroid_effect_low_crp, 4), "\n")
  cat("Odds Ratio:", round(exp(best_low_threshold$steroid_effect_low_crp), 4), "\n")
  cat("P-value:", format(best_low_threshold$steroid_main_p, scientific = TRUE), "\n")
}

if(nrow(beneficial_high) > 0) {
  best_high_threshold <- beneficial_high[which.min(abs(beneficial_high$steroid_effect_high_crp)), ]
  cat("\nBest threshold for HIGH CRP benefit:\n")
  cat("Threshold:", best_high_threshold$threshold, "mg/L\n")
  cat("Steroid effect (log OR):", round(best_high_threshold$steroid_effect_high_crp, 4), "\n")
  cat("Odds Ratio:", round(exp(best_high_threshold$steroid_effect_high_crp), 4), "\n")
}

# Summary of benefit patterns
cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("SUMMARY OF STEROID BENEFIT:\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

if(nrow(beneficial_low) > 0 && nrow(beneficial_high) == 0) {
  cat("Steroids appear beneficial primarily in LOW CRP patients\n")
} else if(nrow(beneficial_low) == 0 && nrow(beneficial_high) > 0) {
  cat("Steroids appear beneficial primarily in HIGH CRP patients\n")
} else if(nrow(beneficial_low) > 0 && nrow(beneficial_high) > 0) {
  cat("Steroids may be beneficial in both groups at different thresholds\n")
} else {
  cat("No clear beneficial effect found at tested thresholds\n")
  cat("Consider: 1) Different threshold ranges, 2) Non-linear relationships, 3) Other confounders\n")
}

# Plot the p-values across thresholds
library(ggplot2)
threshold_plot <- ggplot(threshold_results, aes(x = threshold, y = -log10(interaction_p_value))) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "CRP Threshold Analysis: Interaction Significance",
       x = "CRP Threshold (mg/L)",
       y = "-log10(p-value)",
       subtitle = "Red line indicates p = 0.05") +
  theme_minimal()

print(threshold_plot)

# Save the plot
ggsave("figures/crp_threshold_analysis.png", plot = threshold_plot, 
       width = 10, height = 6, dpi = 300, bg = "white")

cat("\nThreshold analysis plot saved to figures/crp_threshold_analysis.png\n")

