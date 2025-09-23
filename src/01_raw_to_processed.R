# Processes raw data files and filters patients who would then be eligable for analsis.

# Load required libraries
library(tidyverse)
library(readr)

#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

# Loading data for analysis 
UNITE_2020_corrected = read.csv("data/original/UNITE_2020_corrected.csv", header=TRUE)

######################################################## Filter dataset: PRISMA diagram #############################################################################
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  filter( !is.na(ICU_CORTICO_YN)) |> # Filter out patients without corticosteroid information, N = 218
  filter( !is.na(OUTCOME_LD)) |> # Filter out patients without outcome information, N = 111
  filter( !is.na(INC_PREGNANT_YN)) # Filter out patient who are pregnant, N = 0 

#Variables for synethesis 
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  mutate(INC_BMI_INT = INC_WEIGHT_INT / (INC_HEIGHT_INT/100)^2) |> # BMI
  mutate(NEW_CENTRE_ID = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID)) |> #ID for centre
  mutate(NEW_PATIENT_ID = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID, "_", NEW_SUBJECT_ID)) |>  #ID for patient
  mutate(OUT_HOSP_DURATION_OVERALL_INT = ifelse( is.na(OUT_HOSP_DURATION_INT), 
                                                 OUT_ICU_DURATION_INT + INC_LOS_PRIOR_ADM_INT, 
                                                 OUT_HOSP_DURATION_INT)) |> #calculating the overall days in hospital (PreICU +ICU) if a patient died in ICU
  mutate(COAG_THROMBO_COMPLICATION = ifelse( COAG_THROMBO_NONE_CB == FALSE, TRUE, FALSE)) # Any thromboembolic complication


# Severity of respiratory failure denoted by the requirement for respiratory support.

    # replacing NA with FALSE for ventilation data
    UNITE_2020_corrected <- UNITE_2020_corrected %>% replace_na(list(RESP_PRONE_YN = FALSE, RESP_ECMO_YN = FALSE, RESP_NEUROM_BLOC_YN = FALSE, ICU_ANTIMALARIAL_YN=FALSE,
                                                ICU_ANTIVIRALS_YN = FALSE))


    #ventilation_severity - score the ventilation severity using different levels of ventilatory support
    UNITE_2020_corrected <- UNITE_2020_corrected %>% mutate(ventilation_severity = case_when(RESP_ECMO_YN == TRUE ~3, RESP_PRONE_YN==TRUE ~3,
                                                                        RESP_NEUROM_BLOC_YN ==TRUE & RESP_PRONE_YN ==FALSE ~2,
                                                                        RESP_PRONE_YN==FALSE & RESP_NEUROM_BLOC_YN == FALSE ~1))



######################################################### restrict to the variables of interest ###############################################################
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  select(
  #ID 
  "NEW_CENTRE_ID",
  "NEW_PATIENT_ID",
  "NEW_COUNTRY_ID",
  "NEW_SITE_ID",              
  "NEW_SUBJECT_ID",
  # Demographics
  "INC_AGE_INT",
  "INC_SEX_RAD",
  "INC_BMI_INT",
  # Grouping
  "ICU_CORTICO_YN",
  # Comorbidities (range)
  "INC_CARDIAC_DISEASE_YN",
  "INC_LIVER_DISEASE_YN",
  "INC_HBP_YN",
  "INC_NEURO_YN",
  "INC_PULMO_DISEASE_YN",
  "INC_DIABETES_YN",
  "INC_ASTHMA_YN",
  "INC_NEOPLASM_YN",
  "INC_KIDNEY_DISEASE_YN",
  "INC_IMMUNOSUPPR_YN",
  "INC_HIV_YN",
  # Pre-ICU admission
  "ICU_ADM_DIAG_RAD",
  "INC_LOS_PRIOR_ADM_INT",
  "ICU_RESP_SUPPORT_YN",
  "ICU_SUPP_TYPE_RAD",
  "ICU_WHITE_CELL_INT", 
  "ICU_NEUTRO_INT",
  "ICU_FERRITINE_INT", # CRF -Ferritin recorded as mg/L 
  "ICU_DIMERS_INT",
  "ICU_CRP_INT",
  "ICU_PRO_CALCIT_DEC",
  "ICU_PLATELETS_INT",
  "ICU_LYMPH_DEC",
  # ICU medications
  "ICU_ANTIVIRALS_YN",
  "INF_ANTIBIO_YN",
  # Corticosteroids
  "ICU_CORTICO_DURATION_INT",
  "ICU_CORTICO_INTERV_INT",
  "ICU_CORTICO_INDICATION_RAD",
  "INF_AT_ADMISSION_YN",
  "INF_DURING_ICU_YN",
       # "INF_SEVERITY", CRF indicates this is only if developed infection during ICU 
  # During ICU
  "ICU_SUP_WHITE_CELL_INT",
  "ICU_SUP_NEUTRO_INT",
  "ICU_SUP_LYMPH",
  "ICU_SUP_FERRITINE_DEC",
  "COAG_DIMERS_INT",
  "COAG_PLATELETS_INT",
  "ICU_ANTIVIRALS_YN",	
  "ICU_ANTIVIRALS_RAD", # contains Lopinavir/Ritonavir /Remdesivir etc. 
  "ICU_OTHER_ANTIVIRALS_RAD", # contains Convalescent plasma/Tocilizumab / Anakinra / Interferon alpha / Interferon beta
  "ICU_ANTIMALARIAL_YN",
  "ICU_CLIN_TRIAL_YN",
  # Respiratory support
  "RESP_PRONE_YN",
  "RESP_ECMO_YN",
  "RESP_NEUROM_BLOC_YN",
  "RESP_INTUBATED_YN",
  "RESP_INTUB_DAYS_AFT_ADM_INT", # DAYS FROM ICU ADMISSION TO INTUBATION
  "RESP_DURATION_INV_VENT_INT", # DURATION OF VENTILATION
  # Outcomes
  "OUTCOME_LD",
  "OUT_ICU_DURATION_INT",
  "OUT_HOSP_DURATION_INT",
  "OUT_HOSP_DURATION_OVERALL_INT",
  "OUT_DEAD_DURING_ICU_YN",
  # Secondary outcomes
  "ICU_RRT_DIAL_YN",
  "ICU_INOTROPES_YN",
  "RESP_INV_VENT_YN",
  "RESP_NI_VENT_YN",
  "COAG_THROMBO_COMPLICATION",
  "ventilation_severity"
)

######################################################## Ensure format of data is correct ########################################################


              ###########glucocortioid metrics ###########
# Convert categorical variables of the indication of steriods to factors according to the CRF
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(ICU_CORTICO_INDICATION_RAD = case_when(
    ICU_CORTICO_INDICATION_RAD == 1 ~ "Shock",
    ICU_CORTICO_INDICATION_RAD == 2 ~ "Hyperinflammation",
    ICU_CORTICO_INDICATION_RAD == 3 ~ "Pneumonitis",
    ICU_CORTICO_INDICATION_RAD == 4 ~ "Pre-existing condition",
    ICU_CORTICO_INDICATION_RAD == 5 ~ "Other",
    TRUE ~ as.character(ICU_CORTICO_INDICATION_RAD)
  ))

# Convert ICU_CORTICO_INDICATION_RAD to muiltile binary variables
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(
    ICU_CORTICO_INDICATION_SHOCK = ICU_CORTICO_INDICATION_RAD == "Shock",
    ICU_CORTICO_INDICATION_HYPERINFLAMMATION = ICU_CORTICO_INDICATION_RAD == "Hyperinflammation",
    ICU_CORTICO_INDICATION_PNEUMONITIS = ICU_CORTICO_INDICATION_RAD == "Pneumonitis",
    ICU_CORTICO_INDICATION_PRE_EXISTING_CONDITION = ICU_CORTICO_INDICATION_RAD == "Pre-existing condition",
    ICU_CORTICO_INDICATION_OTHER = ICU_CORTICO_INDICATION_RAD == "Other"
  )

# Indicator of those who started steroids prior to ICU admission
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate( 
    ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN = ifelse(ICU_CORTICO_INTERV_INT == 0, 1, 0), #Coritcoids started at admission
    ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN = ifelse(ICU_CORTICO_INTERV_INT <= INC_LOS_PRIOR_ADM_INT, 1, 0), # corticoids started prior to icu
    ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN = ifelse( INC_LOS_PRIOR_ADM_INT <= ICU_CORTICO_INTERV_INT, 1, 0) # corticoids started during icu
  )

# how many do not have admission to ICU and admission to corticoid information and receive steroids
#table(ifelse(is.na(UNITE_2020_corrected$ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN) &
#     is.na(UNITE_2020_corrected$ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_ICU_YN) &
#      UNITE_2020_corrected$ICU_CORTICO_YN == TRUE, TRUE, FALSE))

##################################################################
# Convert OUTCOME_LD to multiple indicator (dummy) variables
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(
    OUTCOME_LD_DEATH = OUTCOME_LD == "Death",
    OUTCOME_LD_DISCHARGED = OUTCOME_LD == "Discharged",
    OUTCOME_LD_TRANSFERRED = OUTCOME_LD == "Transferred",
    OUTCOME_LD_OTHER = OUTCOME_LD == "Other"
  )

# replacing NA with FALSE for ventilation data - again assuming TRUE would hav ebeen entered if it had been provided
UNITE_2020_corrected<- UNITE_2020_corrected %>% replace_na(list(RESP_PRONE_YN = FALSE, RESP_ECMO_YN = FALSE, RESP_NEUROM_BLOC_YN = FALSE, ICU_ANTIMALARIAL_YN=FALSE,
                                                ICU_ANTIVIRALS_YN = FALSE))


# converting character variables to factors
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(
    INC_SEX_RAD = as.factor(INC_SEX_RAD),
    INC_DIABETES1_YN = ifelse(INC_DIABETES_YN =="Type II diabetes", TRUE, FALSE),
    INC_DIABETES2_YN = ifelse(INC_DIABETES_YN =="Type I diabetes", TRUE, FALSE)

  )

################################################################### Variables which are NA not allowed to be #################################################################

UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(ICU_OTHER_ANTIVIRALS_RAD = ifelse(is.na(ICU_OTHER_ANTIVIRALS_RAD), "None", ICU_OTHER_ANTIVIRALS_RAD)) |>
  mutate(ICU_ANTIVIRALS_RAD = ifelse(is.na(ICU_ANTIVIRALS_RAD), "None", ICU_ANTIVIRALS_RAD))

# Save the processed data
write.csv(UNITE_2020_corrected, "data/processed/UNITE_2020_corrected_processed.csv", row.names = FALSE)

# missingess in ICU_SUPP_TYPE_RAD
 UNITE_2020_corrected <- UNITE_2020_corrected |> 
  mutate(ICU_SUPP_TYPE_RAD = ifelse( ICU_RESP_SUPPORT_YN  == FALSE, "No respiratory support prior to ICU", ICU_SUPP_TYPE_RAD))

# Organ support outcome - assumed to be FALSE if NA
 UNITE_2020_corrected <- UNITE_2020_corrected |> 
  mutate(
    RESP_NI_VENT_YN = ifelse(is.na(RESP_NI_VENT_YN), FALSE, RESP_NI_VENT_YN),
    RESP_INTUBATED_YN  = ifelse(is.na(RESP_INTUBATED_YN), FALSE, RESP_INTUBATED_YN),
    ICU_INOTROPES_YN = ifelse(is.na(ICU_INOTROPES_YN), FALSE, ICU_INOTROPES_YN),
    ICU_RRT_DIAL_YN = ifelse(is.na(ICU_RRT_DIAL_YN), FALSE, ICU_RRT_DIAL_YN)
  )


####################################### Creation of new variables #######################################

# Stratification based on admission CRP  by 100
 UNITE_2020_corrected =  UNITE_2020_corrected |> 
    mutate(ICU_CRP_RAD = cut(ICU_CRP_INT,
        breaks = c( 0 ,100, 200, 300, Inf ),
        labels = c( "0-100", "101-200", "201-300", "≥300" )
    ))

# Stratification based on admission CRP  by 50
 UNITE_2020_corrected =  UNITE_2020_corrected |> 
    mutate(ICU_CRP_RAD50 = cut(ICU_CRP_INT,
        breaks = c( 0 ,50, 100, 150, 200, 250, 300, Inf ),
        labels = c( "0-50", "51-100", "101-150", "151-200", "201-250", "251-300", "≥301" )
    ))

# comorbidity score - replace NA with FALSE to allow for scoring (Taken from cambridge)
UNITE_2020_corrected<- UNITE_2020_corrected %>% replace_na(list(INC_CARDIAC_DISEASE_YN = FALSE, INC_PULMO_DISEASE_YN = FALSE,
                                                INC_ASTHMA_YN = FALSE, INC_DIABETES_YN = "No diabetes", INC_NEURO_YN = FALSE,
                                                INC_HBP_YN = FALSE, INC_NEOPLASM_YN = FALSE, INC_LIVER_DISEASE_YN = FALSE,
                                                INC_IMMUNOSUPPR_YN = FALSE, INC_HIV_YN = "No HIV", INC_KIDNEY_DISEASE_YN= FALSE))

UNITE_2020_corrected<- UNITE_2020_corrected %>% mutate(cardiac_d = if_else(INC_CARDIAC_DISEASE_YN ==TRUE, 1,0),
                                       liver_d = if_else(INC_LIVER_DISEASE_YN ==TRUE, 1,0),
                                       neuro_d = if_else(INC_NEURO_YN ==TRUE, 1,0),
                                       diabetes_d = if_else(INC_DIABETES_YN !="No diabetes", 1,0),
                                       kidney_d = if_else(INC_KIDNEY_DISEASE_YN ==TRUE, 1,0),
                                       htn_d = if_else(INC_HBP_YN ==TRUE, 1,0),
                                       asthma_d = if_else(INC_ASTHMA_YN ==TRUE, 1,0),
                                       resp_d = if_else(INC_CARDIAC_DISEASE_YN ==TRUE, 1,0),
                                       immunosup_d = if_else(INC_IMMUNOSUPPR_YN ==TRUE, 1,0),
                                       hiv_d = if_else(INC_HIV_YN !="No HIV", 1,0))

# compute an additive comorbidity score where each comorbidity has equal weight:
UNITE_2020_corrected<- UNITE_2020_corrected %>% mutate(comorbidity_score = cardiac_d+liver_d+neuro_d+diabetes_d+kidney_d+htn_d+asthma_d+resp_d+immunosup_d+hiv_d)

# neutrophil to lymphocyte ratio

UNITE_2020_corrected <- UNITE_2020_corrected %>%
  mutate(neutrophil_lymphocyte_ratio = ifelse(ICU_LYMPH_DEC == 0, NA, ICU_NEUTRO_INT / ICU_LYMPH_DEC))
