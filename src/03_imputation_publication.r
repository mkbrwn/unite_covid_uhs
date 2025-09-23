#Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
# Exploration of missingness

#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

# Run 01_raw_to_processed.r to process the data

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "mice")

source("src/01_raw_to_processed.r") 

#library 
library(gtsummary)
library(mice) #muiltiple imputation
library(miceadds) 
library(future.apply) # parrallelisation of MICE 
library(naniar) #visualisation of missingness
library(caret) # For colinearity detection
library(mixgb) # mi with xgboost
library(MASS) # backwards step AIC 

# make sure select is used with dplyr by default 


# Visualise missingness matrix
vis_miss(UNITE_2020_corrected)

# Visualise missingness in the dataset
missingness_plot = gg_miss_var(UNITE_2020_corrected, show_pct = TRUE) +
  labs(title = "Missingness in UNITE_2020_corrected dataset") +
  theme_minimal() 

# Save the plot
ggsave("figures/missingness_plot.png", plot = missingness_plot, width = 8, height = 6, dpi = 300,bg = "white")

# compare missingness between ICU_CORTICO_YN variable 
missingness_by_steroids_plot = gg_miss_var(UNITE_2020_corrected, show_pct = TRUE, facet = ICU_CORTICO_YN)
ggsave("figures/missingness_by_steroids_plot.png", plot = missingness_by_steroids_plot, width = 8, height = 6, dpi = 300,bg = "white")

# Other graphs of missingness 
vis_miss(UNITE_2020_corrected) # 10% missingness 
gg_miss_fct(UNITE_2020_corrected, NEW_SITE_ID) 
gg_miss_fct(UNITE_2020_corrected, NEW_SITE_ID) 

print("graphical exploration of missingness complete")

# remove unessasary variables 


vars_to_remove <- c(
  "ICU_CORTICO_INDICATION_SHOCK",
  "ICU_CORTICO_INDICATION_HYPERINFLAMMATION",
  "ICU_CORTICO_INDICATION_PNEUMONITIS",
  "ICU_CORTICO_INDICATION_PRE_EXISTING_CONDITION",
  "ICU_CORTICO_INDICATION_OTHER",
  "ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN",
  "ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN",
  "ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN",
  "cardiac_d", "liver_d", "neuro_d", "diabetes_d", "kidney_d", "htn_d",
  "asthma_d", "resp_d", "immunosup_d", "hiv_d",
  "NEW_SUBJECT_ID", "NEW_CENTRE_ID", "NEW_PATIENT_ID", "NEW_COUNTRY_ID",
  "ICU_SUP_WHITE_CELL_INT", "ICU_SUP_NEUTRO_INT", "ICU_SUP_LYMPH", "ICU_SUP_FERRITINE_DEC",
  "COAG_DIMERS_INT", "COAG_PLATELETS_INT",
  "RESP_PRONE_YN", "RESP_ECMO_YN", "RESP_NEUROM_BLOC_YN", "RESP_INTUBATED_YN",
  "RESP_INTUB_DAYS_AFT_ADM_INT", "RESP_DURATION_INV_VENT_INT", "RESP_INV_VENT_YN", "RESP_NI_VENT_YN",
  "OUTCOME_LD_OTHER", "OUTCOME_LD_TRANSFERRED", "OUTCOME_LD_DISCHARGED", "OUT_HOSP_DURATION_OVERALL_INT",
  "ICU_ADM_DIAG_RAD", "ICU_CORTICO_INDICATION_RAD", "ICU_ANTIVIRALS_RAD", "OUTCOME_LD",
  "ICU_OTHER_ANTIVIRALS_RAD", "INC_HIV_YN", "ICU_SUPP_TYPE_RAD", "INC_DIABETES_YN",
  "OUTCOME_LD_DEATH", "OUT_HOSP_DURATION_INT", "OUT_ICU_DURATION_INT",
  "ICU_CRP_RAD50", "ICU_CRP_RAD"
)

UNITE_2020_corrected_preimputation <-
  UNITE_2020_corrected[ , 
    ! names(UNITE_2020_corrected) %in% vars_to_remove
  ]


  print("Preimputation data set cleaned")

#missingness graph 

missingness_plot = gg_miss_var(UNITE_2020_corrected_preimputation, show_pct = TRUE) +
  labs(title = "Missingness in UNITE_2020_corrected dataset") +
  theme_minimal() 

gg_miss_fct(UNITE_2020_corrected_preimputation,NEW_SITE_ID) 

# exclusion of variables with >25% missingness
UNITE_2020_corrected_preimputation <- as_tibble(UNITE_2020_corrected_preimputation) |> 
  dplyr::select(where(~ mean(is.na(.)) < 0.25))

#Remove colinear variables 
# 1. Select only numeric cols
num_df <- UNITE_2020_corrected_preimputation %>%
  dplyr::select(where(is.numeric))

# 2. Compute correlation matrix
cor_mat <- cor(num_df, use = "pairwise.complete.obs")

# 3. Identify variables to remove (e.g. cutoff = 0.8)
to_remove_idx <- findCorrelation(cor_mat, cutoff = 0.8)

# 4. Get their names
to_remove_names <- colnames(num_df)[to_remove_idx]

# 5. Drop them from the original data frame
df_uncorrelated <- UNITE_2020_corrected_preimputation %>%
  dplyr::select(-all_of(to_remove_names))

       # give missingness  graph here

missing_preimputation  =  gg_miss_fct(UNITE_2020_corrected_preimputation,NEW_SITE_ID) 
ggsave("figures/missing_preimputation.png", plot = missing_preimputation , width = 8, height = 6, dpi = 300, bg = "white")


##############################################imputation##############################################

UNITE_2020_corrected_imputation = futuremice(UNITE_2020_corrected_preimputation, m=5, n.cores=8 ) 
print("multiple imputation completed")

# synthesis of vriables post imputation

# neutrophil_lymphocyte_ratio, comorbidity_score,
# INC_DIABETES1_YN,INC_DIABETES2_YN, ventilation_severity






### Diagnositcs of imputed method
 
# Plot convergence for all imputed variables
pdf("figures/convergence_imputation.pdf", width = 10, height = 8)
plot(UNITE_2020_corrected_imputation) # convergence looks adequate 
dev.off()

#bwplot and save
png("figures/bwplot_imputation.png", width = 8, height = 6, units = "in", res = 300, bg = "white")
bwplot(UNITE_2020_corrected_imputation)
dev.off()

pdf("figures/densityplot_imputation.pdf", width = 8, height = 6, units = "in", res = 300, bg = "white")
densityplot_imputation = densityplot(UNITE_2020_corrected_imputation)
dev.off()

print("Diagnostic plots of multiple imputation completed")