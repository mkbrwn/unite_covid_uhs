# Processes raw data files and filters patients who would then be eligable for analsis.

# Load required libraries
library(tidyverse)
library(readr)
library(gtsummary)

# Loading data for analysis 
UNITE_2020_corrected = read.csv("data/original/UNITE_2020_corrected.csv", header=TRUE)
UNITE_2020_corrected = UNITE_2020_corrected %>% mutate(wave = 1)
max = max(UNITE_2020_corrected$NEW_SITE_ID) # used for unique ID assignment

UNITE_2021_corrected = read.csv("data/original/UNITE_2021_corrected.csv", header=TRUE)
UNITE_2021_corrected = UNITE_2021_corrected %>% mutate(wave = 2)

# Combine datasets
UNITE_2020_corrected = bind_rows(UNITE_2020_corrected, UNITE_2021_corrected)

# give unique SITE ID between first and second wave
UNITE_2020_corrected = UNITE_2020_corrected %>% 
      mutate(NEW_SITE_ID = ifelse(wave ==2,
      NEW_SITE_ID + max, 
      NEW_SITE_ID)
      )

# to restrict to wave 2 only
#UNITE_2020_corrected = UNITE_2020_corrected %>% filter(wave == 2)


######################################################## Filter dataset: PRISMA diagram #############################################################################
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  filter( !is.na(ICU_CORTICO_YN)) |> # Filter out patients without corticosteroid information, N = 218
  filter( !is.na(OUTCOME_LD))  # Filter out patients without outcome information, N = 111

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

UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(ICU_TRACHEOS_YN = ifelse(ICU_TRACHEOS_YN == "Before this admission", "Yes", "No"))


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

# convert ICU_ATELECTASIS_YN, ICU_PRESSURE_FAC_YN , ICU_PRESSURE_OTH_YN,  ICU_OBSTRUCTION_YN to boolean from strings
UNITE_2020_corrected <- UNITE_2020_corrected %>%
  mutate(
    ICU_ATELECTASIS_YN = case_when(
      tolower(ICU_ATELECTASIS_YN) == "true" ~ TRUE,
      tolower(ICU_ATELECTASIS_YN) == "false" ~ FALSE,
      TRUE ~ NA
    ),
    ICU_PRESSURE_FAC_YN = case_when(
      tolower(ICU_PRESSURE_FAC_YN) == "true" ~ TRUE,
      tolower(ICU_PRESSURE_FAC_YN) == "false" ~ FALSE,
      TRUE ~ NA
    ),
    ICU_PRESSURE_OTH_YN = case_when(
      tolower(ICU_PRESSURE_OTH_YN) == "true" ~ TRUE,
      tolower(ICU_PRESSURE_OTH_YN) == "false" ~ FALSE,
      TRUE ~ NA
    ),
    ICU_OBSTRUCTION_YN = case_when(
      tolower(ICU_OBSTRUCTION_YN) == "true" ~ TRUE,
      tolower(ICU_OBSTRUCTION_YN) == "false" ~ FALSE,
      TRUE ~ NA
    )
  )


# remove unessasary variables 

vars_to_remove <- c(
  #Related to the glucocorticoids 
  "ICU_CORTICO_INDICATION_SHOCK",
  "ICU_CORTICO_INDICATION_HYPERINFLAMMATION",
  "ICU_CORTICO_INDICATION_PNEUMONITIS",
  "ICU_CORTICO_INDICATION_PRE_EXISTING_CONDITION",
  "ICU_CORTICO_INDICATION_OTHER",
  "ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN",
  "ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN",
  "ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN",
  "ICU_CORTICO_DURATION_INT",
  "ICU_CORTICO_INTERV_INT",

  # comorbidity dummy variables
  "cardiac_d", "liver_d", "neuro_d", "diabetes_d", "kidney_d", "htn_d",
  "asthma_d", "resp_d", "immunosup_d", "hiv_d",

  #ID variables 
  "NEW_SUBJECT_ID", "NEW_CENTRE_ID", "NEW_PATIENT_ID", "NEW_COUNTRY_ID", "Iteration",
  "INC_CRIT1_YN", "INC_CRIT2_YN", "INC_CRIT3_YN", "INC_CRIT4_YN", "DURATION_INC_CRIT4_YN",
  "INC_PREGNANT_YN","INC_HELTH_WORK_YN", "DURATION_LOSINC_HELTH_WORK_YN",

  #superceded by other variables
  "INC_HEIGHT_INT", "INC_WEIGHT_INT",

  #Labs during ICU stay 
  "ICU_SUP_WHITE_CELL_INT", "ICU_SUP_NEUTRO_INT", "ICU_SUP_LYMPH", "ICU_SUP_FERRITINE_DEC",
  "COAG_DIMERS_INT", "COAG_PLATELETS_INT", "COAG_ANTIPLAT_DOSE_TXT",

  #Coagulation variables to exclude 
  "COAG_DVT_DOSE_TXT", "COAG_LIFETHREAT_HEMO_YN", "COAG_TRANSF_NB_INT", "NEW_COAG_DVT_DAILY_DOSE", "INF_WITHOUT_ANTIMIC_INT",

  # variables for ventilation parameters
   "RESP_TIDAL_INT",
  "RESP_PEEP_INT",
  "RESP_FIO2_INT",
  "RESP_PF_RATIO_INT",
  "RESP_PACO2_INT",
  "RESP_DRIV_PRESS_INT",
  "RESP_PRONE_INTUB_INT",
  "RESP_PRONE_NOTINTUB_INT",


  #Outcomes for exclusion 
  "OUTCOME_LD_OTHER", "OUTCOME_LD_TRANSFERRED", "OUTCOME_LD_DISCHARGED", "OUT_HOSP_DURATION_OVERALL_INT",
  "ICU_ADM_DIAG_RAD", "ICU_ADM_DIAG_YN", "ICU_ANTIVIRALS_RAD", "OUTCOME_LD",
  "ICU_OTHER_ANTIVIRALS_RAD", "INC_HIV_YN", "ICU_SUPP_TYPE_RAD", "INC_DIABETES_YN",
  "OUTCOME_LD_DEATH", "OUT_HOSP_DURATION_INT", "OUT_ICU_DURATION_INT",
  "ICU_CRP_RAD50", "ICU_CRP_RAD",
  "INC_HEIGHT", "ICU_PRESSURE_FAC_YN", "ICU_TEMPERATURE_DEC", "RESP_INTUBATED_YNTRUE",

  # not suitable for modeling
  "ICU_TRACHEO_METHOD_RAD",
  "INF_ANTIFUNG4_INT",
  "INF_ANTIBIO1_INT",
  "INF_ANTIBIO2_INT",

  # infection parameters
  "INF_PULMO_YN",
  "INF_RESPIR_YN",
  "INF_ABDO_YN",
  "INF_BACTEREMIA_YN",
  "INF_URINARY_YN",
  "INF_CNS_YN",
  "INF_OTHER_YN",
  "INF_CLABSI_YN",
  "INF_MDR_PATHO_YN",

  # irrelevant
  "REH_MOBIL_72H_YN",
  "REH_MOBIL_VENT_72H_YN"
)

UNITE_2020_corrected <- UNITE_2020_corrected[ , 
    ! names(UNITE_2020_corrected) %in% vars_to_remove
  ]

# exclude any variable with > 25% missingness
  missingness <- colMeans(is.na(UNITE_2020_corrected)) > 0.25
  vars_to_remove <- c(vars_to_remove, names(missingness[missingness]))

UNITE_2020_corrected <-
  UNITE_2020_corrected[ , 
    ! names(UNITE_2020_corrected) %in% vars_to_remove
  ]

#exclude any variable with "CRIT"
vars_to_remove <- c(vars_to_remove, grep("CRIT", names(UNITE_2020_corrected), value = TRUE))

#exclude variables with "COAG_"
vars_to_remove <- c(vars_to_remove, grep("COAG_", names(UNITE_2020_corrected), value = TRUE))

#exclude variables with "COAG_BLEED_"
vars_to_remove <- c(vars_to_remove, grep("COAG_BLEED_", names(UNITE_2020_corrected), value = TRUE))

#remove variables which end in "_UNK"
vars_to_remove <- c(vars_to_remove, grep("_UNK$", names(UNITE_2020_corrected), value = TRUE))


vars_to_remove <- c(vars_to_remove, grep("ICU_THROMBO", names(UNITE_2020_corrected), value = TRUE))

# exclude any variable with "RAD"
#rad_vars <- grep("RAD", names(UNITE_2020_corrected), value = TRUE)
#vars_to_remove <- c(vars_to_remove, rad_vars)


UNITE_2020_corrected <-
  UNITE_2020_corrected[ , 
    ! names(UNITE_2020_corrected) %in% vars_to_remove
  ]


print("data cleaned")



#####################################  PRODUCE categories of CRP #####################################

UNITE_2020_corrected <- UNITE_2020_corrected %>%
  mutate(
    ICU_CRP_CATEGORY = case_when(
      ICU_CRP_INT < 100 ~ "Low",
      ICU_CRP_INT >= 100 & ICU_CRP_INT < 200 ~ "Moderate",
      ICU_CRP_INT >= 200 & ICU_CRP_INT < 300 ~ "High",
      ICU_CRP_INT >= 300 ~ "Very High",
      TRUE ~ "Unknown"
    )
  )


#test and split data using caret
library(caret)  
set.seed(123)  # For reproducibility
train_indices <- createDataPartition(UNITE_2020_corrected$OUT_DEAD_DURING_ICU_YN, p = 0.8, list = FALSE)
train_data <- UNITE_2020_corrected[train_indices, ]
test_data <- UNITE_2020_corrected[-train_indices, ]
