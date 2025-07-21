# Processes raw data files and filters patients who would then be eligable for analsis.

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)

# Loading data for analysis 
UNITE_2020_corrected = read.csv("data/original/UNITE_2020_corrected.csv", header=TRUE)

# PRISMA diagram - filter 
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

# filter data set to the variables of interest
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
  "INC_HEIGHT_INT",
  "INC_WEIGHT_INT",
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
  "ICU_RESP_SUPPORT_YN",
  "ICU_SUPP_TYPE_RAD",
  "ICU_WHITE_CELL_INT", 
  "ICU_NEUTRO_INT",
  "ICU_FERRITINE_INT",
  "ICU_DIMERS_INT",
  "ICU_CRP_INT",
  "ICU_PRO_CALCIT_DEC",
  "ICU_PLATELETS_INT",
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
  # Outcomes
  "OUTCOME_LD",
  "OUT_ICU_DURATION_INT",
  "OUT_HOSP_DURATION_OVERALL_INT",
  # Secondary outcomes
  "ICU_RRT_DIAL_YN",
  "ICU_INOTROPES_YN",
  "RESP_INTUBATED_YN",
  "RESP_NI_VENT_YN",
  "COAG_THROMBO_COMPLICATION"
)

# Ensure format of data is correct

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

# Convert OUTCOME_LD to multiple indicator (dummy) variables
UNITE_2020_corrected <- UNITE_2020_corrected |>
  mutate(
    OUTCOME_LD_DEATH = OUTCOME_LD == "Death",
    OUTCOME_LD_DISCHARGED = OUTCOME_LD == "Discharged",
    OUTCOME_LD_TRANSFERRED = OUTCOME_LD == "Transferred",
    OUTCOME_LD_OTHER = OUTCOME_LD == "Other"
  )
# Save the processed data
write.csv(UNITE_2020_corrected, "data/processed/UNITE_2020_corrected_processed.csv", row.names = FALSE)