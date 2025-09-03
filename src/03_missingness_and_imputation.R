# Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
# Exploration of missingness

#library 
library(tidyverse)
library(gtsummary)
library(mice) #muiltiple imputation 
library(future.apply) # parrallelisation of MICE 
library(naniar) #visualisation of missingness


#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.r") 

print("Exploring missingness in UNITE_2020_corrected dataset")  

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

# exploration of missingness patterns

# work in progress  - should identify important variables - I suspect that missingness will be related center.
# It is important to note some variables have high proportion of missingess (>25%). E.g
  # ICU_PRO_CALCIT_DEC
  # ICU_PRO_FERRITINE_INT
  # ICU_PRO_DIMERS_INT
  # COAG_DIMERS_INT
  # ICU_SUP_FerrITINE_INT 


#removal of corticosteroid related variables
#UNITE_2020_corrected_mice = UNITE_2020_corrected |> 
#         select(-c(ICU_CORTICO_DURATION_INT, ICU_CORTICO_INTERV_INT, ICU_CORTICO_INDICATION_RAD, ICU_CORTICO_INDICATION_SHOCK,
#           ICU_CORTICO_INDICATION_HYPERINFLAMMATION, ICU_CORTICO_INDICATION_PNEUMONITIS, ICU_CORTICO_INDICATION_PRE_EXISTING_CONDITION, ICU_CORTICO_INDICATION_OTHER))



######################################################### Regression with steroids with adjustment for CRP ##############################################

# Complete case analysis

  # univariable logstic regression 
complete_case_glm = glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN, data = UNITE_2020_corrected,family = binomial)
summary(complete_case_glm)
complete_case_glm = tbl_regression(complete_case_glm, exponentiate = TRUE)
complete_case_glm <- tibble::as_tibble(as.data.frame(complete_case_glm))


  # logistic regression with adjustment for CRP
complete_case_glm_crp = glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ICU_CRP_INT +comorbidity_score , data = UNITE_2020_corrected,family = binomial)
complete_case_glm_crp = tbl_regression(complete_case_glm_crp, exponentiate = TRUE)
complete_case_glm_crp <- tibble::as_tibble(as.data.frame(complete_case_glm_crp))

print("complete case analysis")
####################################################################### imputation ######################################################################

#Minimal set used for imputation 
UNITE_2020_corrected_mice_minimum = UNITE_2020_corrected |> 
          select( ICU_CORTICO_YN, ICU_CRP_INT, ICU_CRP_RAD, INC_BMI_INT, NEW_PATIENT_ID, NEW_SITE_ID, comorbidity_score,  OUTCOME_LD_DEATH)

#muiltiple imputation with parallelisation
UNITE_2020_corrected_mice_minimum = futuremice(UNITE_2020_corrected_mice_minimum,  m=10, n.cores=4)
print("Imputation complete")

#pool models and results across imputated data sets
model_imputed = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN,  family = binomial))
imputation_glm = tbl_regression(model_imputed, exponentiate = TRUE)
imputation_glm <- tibble::as_tibble(as.data.frame(imputation_glm))

model_imputed = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ICU_CRP_INT,  family = binomial))
imputation_glm_crp = tbl_regression(model_imputed, exponentiate = TRUE)
imputation_glm_crp <- tibble::as_tibble(as.data.frame(imputation_glm_crp))

model_imputed = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ICU_CRP_RAD,  family = binomial))
imputation_glm_crp_cat = tbl_regression(model_imputed, exponentiate = TRUE)
imputation_glm_crp_cat <- tibble::as_tibble(as.data.frame(imputation_glm_crp_cat))



# imputation diagnostics
densityplot(UNITE_2020_corrected_mice_minimum)

#export regression analysis to excel
library(openxlsx)

# Create a new workbook
wb_regression <- createWorkbook()

# Add a worksheet for each imputation model
addWorksheet(wb_regression, "Complete case Model ")
addWorksheet(wb_regression, "Complete case Model CRP")
addWorksheet(wb_regression, "Imputation Model")
addWorksheet(wb_regression, "Imputation Model CRP")
addWorksheet(wb_regression, "Imputation Model CRP Category")

writeData(wb_regression, "Complete case Model ", complete_case_glm)
writeData(wb_regression, "Complete case Model CRP", complete_case_glm_crp)
writeData(wb_regression, "Imputation Model", imputation_glm)
writeData(wb_regression, "Imputation Model CRP", imputation_glm_crp)
writeData(wb_regression, "Imputation Model CRP Category", imputation_glm_crp_cat)


# Save the workbook
saveWorkbook(wb_regression, "data/processed/regression_models.xlsx", overwrite = TRUE)
print("regression models saved to excel")