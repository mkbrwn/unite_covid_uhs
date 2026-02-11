#Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
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

# Visualise missingness matrix
vis_miss(train_data)

# resolution of conflicted packages 
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "mice")

# Visualise missingness in the dataset
missingness_plot = gg_miss_var(train_data, show_pct = TRUE) +
  labs(title = "Missingness in train_data dataset") +
  theme_minimal() 

# Save the plot
ggsave("figures/missingness_plot.png", plot = missingness_plot, width = 8, height = 10, dpi = 300,bg = "white")

# compare missingness between ICU_CORTICO_YN variable 
missingness_by_steroids_plot = gg_miss_var(UNITE_2020_corrected, show_pct = TRUE, facet = ICU_CORTICO_YN)
ggsave("figures/missingness_by_steroids_plot.png", plot = missingness_by_steroids_plot, width = 8, height = 10, dpi = 300,bg = "white")

# Other graphs of missingness 
vis_miss(train_data) # 10% missingness 

print("graphical exploration of missingness complete")

#missingness graph 
missingness_plot = gg_miss_var(train_data, show_pct = TRUE) +
  labs(title = "Missingness in train_data dataset") +
  theme_minimal() 

gg_miss_fct(train_data, NEW_SITE_ID) 

# exclusion of variables with >25% missingness
train_data <- as_tibble(train_data) |> 
  dplyr::select(where(~ mean(is.na(.)) < 0.25))

       # give missingness  graph here

missing_preimputation  =  gg_miss_fct(train_data,NEW_SITE_ID) 
ggsave("figures/missing_preimputation.png", plot = missing_preimputation , dpi = 300, bg = "white")

################################################### IMPUTATION ###################################################

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
})

# Step 2: Convert back to mids object using datlist2mids
library(miceadds)
train_data_imputation <- datlist2mids(modified_list)

print("Age groups created and integrated into train_data_imputation mids object")


################################################### POST IMPUTATION ASSSESSMENT ###################################################

### Diagnostic plots of imputed method

# Plot convergence for all imputed variables
pdf("figures/convergence_imputation.pdf", width = 10, height = 8)
plot(train_data_imputation) # convergence looks adequate 
dev.off()

# Create bwplot (box-and-whisker plots) to compare distributions
tryCatch({
  # Set encoding to handle UTF-8 issues
  old_locale <- Sys.getlocale("LC_CTYPE")
  Sys.setlocale("LC_CTYPE", "C")
  
  png("figures/bwplot_imputation.png", width = 10, height = 8, units = "in", res = 300, bg = "white")
  bwplot_imputation <- bwplot(train_data_imputation)
  print(bwplot_imputation)
  dev.off()
  
  # Restore original locale
  Sys.setlocale("LC_CTYPE", old_locale)
  
}, error = function(e) {
  cat("Error creating bwplot:", e$message, "\n")
  cat("Attempting alternative approach...\n")
  
  # Close any open devices
  if(dev.cur() > 1) dev.off()
  
  # Try with different approach - individual variable plots
  tryCatch({
    png("figures/bwplot_imputation_alternative.png", width = 12, height = 10, units = "in", res = 300, bg = "white")
    
    # Get imputed variables
    imputed_vars <- names(train_data_imputation$imp)
    
    if(length(imputed_vars) > 0) {
      # Create layout for multiple plots
      n_vars <- min(6, length(imputed_vars))
      par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
      
      for(i in 1:n_vars) {
        var_name <- imputed_vars[i]
        if(!is.null(train_data_imputation$imp[[var_name]])) {
          # Get completed data for this variable
          complete_data <- complete(train_data_imputation, action = "long")
          
          # Create boxplot comparing imputed vs observed
          if(var_name %in% names(complete_data)) {
            boxplot(complete_data[[var_name]] ~ complete_data$.imp,
                   main = paste("Imputation Check:", var_name),
                   xlab = "Imputation Number (0 = Original)",
                   ylab = var_name,
                   col = rainbow(6))
          }
        }
      }
    } else {
      plot(1, 1, main = "No imputed variables to plot", type = "n")
    }
    
    dev.off()
    cat("Alternative bwplot created successfully\n")
    
  }, error = function(e2) {
    cat("Alternative bwplot also failed:", e2$message, "\n")
    if(dev.cur() > 1) dev.off()
  })
})

# Create densityplot with error handling
tryCatch({
  png("figures/densityplot_imputation.png", width = 10, height = 8, units = "in", res = 300, bg = "white")
  densityplot_imputation <- densityplot(train_data_imputation)
  print(densityplot_imputation)
  dev.off()
}, error = function(e) {
  cat("Error creating densityplot:", e$message, "\n")
  cat("Attempting stripplot instead...\n")
  
  png("figures/stripplot_imputation.png", width = 10, height = 8, units = "in", res = 300, bg = "white")
  stripplot_imputation <- stripplot(train_data_imputation)
  print(stripplot_imputation)
  dev.off()
})

# Create calibration plot to assess quality of imputation
# This compares observed vs imputed data distributions
tryCatch({
  png("figures/calibration_plot_imputation.png", width = 12, height = 10, units = "in", res = 300, bg = "white")
  
  # Use VIM package for calibration plots if available
  if(require(VIM, quietly = TRUE)) {
    # Margin plot for key variables
    marginplot(train_data_imputation$imp$ICU_CRP_INT, 
               main = "Calibration Plot: ICU CRP")
  } else {
    # Alternative using mice built-in functions
    # Create comparison plots for key continuous variables
    library(lattice)
    
    # Get list of imputed variables
    imputed_vars <- names(train_data_imputation$imp)
    
    if(length(imputed_vars) > 0) {
      # Create a multi-panel plot for first few imputed variables
      n_vars <- min(4, length(imputed_vars))
      par(mfrow = c(2, 2))
      
      for(i in 1:n_vars) {
        var_name <- imputed_vars[i]
        if(!is.null(train_data_imputation$imp[[var_name]])) {
          # Plot observed vs imputed distributions
          hist(complete(train_data_imputation)[[var_name]], 
               main = paste("Calibration:", var_name),
               xlab = var_name,
               col = "lightblue",
               alpha = 0.7)
        }
      }
    } else {
      plot(1, 1, main = "No imputed variables found", type = "n")
    }
  }
  
  dev.off()
}, error = function(e) {
  cat("Error creating calibration plot:", e$message, "\n")
})

print("Diagnostic plots of multiple imputation completed")


