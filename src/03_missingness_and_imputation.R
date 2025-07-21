# Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
# Exploration of missingness

#library 
library(tidyverse)
library(gtsummary)
library(mice) #muiltiple imputation 
library(naniar) #visualisation of missingness

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.r") 

# Visualise missingness in the dataset
gg_miss_var(UNITE_2020_corrected, show_pct = TRUE) +
  labs(title = "Missingness in UNITE_2020_corrected dataset") +
  theme_minimal()       