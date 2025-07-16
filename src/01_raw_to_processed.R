# Processes raw data files and filters patients who would then be eligable for analsis.

# Load required libraries
library(dplyr)
library(tidyr)
library(readr)

# Loading data for analysis 
UNITE_2020_corrected = read.csv("data/original/UNITE_2020_corrected.csv", header=TRUE)
UNITE_2021_corrected = read.csv("data/original/UNITE_2021_corrected.csv", header=TRUE)

# PRISMA diagram - filter 
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  filter( !is.na(ICU_CORTICO_YN)) |> # Filter out patients without corticosteroid information, N = 218
  filter( !is.na(OUTCOME_LD))  # Filter out patients without outcome information, N = 111

#Variables for synethesis 
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  mutate(INC_BMI_INT = INC_WEIGHT_INT / (INC_HEIGHT_INT/100)^2) |> # BMI
  mutate(centre = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID)) |> #ID for centre
  mutate(patient = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID, "_", NEW_SUBJECT_ID)) |>  #ID for patient
  mutate(OUT_HOSP_DURATION_OVERALL_INT = ifelse( is.na(OUT_HOSP_DURATION_INT), 
                                                 OUT_ICU_DURATION_INT + INC_LOS_PRIOR_ADM_INT, 
                                                 OUT_HOSP_DURATION_INT)) |> #calculating the overall days in hospital (PreICU +ICU) if a patient died in ICU
  mutate(COAG_THROMBO_COMPLICATION = ifelse( COAG_THROMBO_NONE_CB == FALSE, TRUE, FALSE)) # Any thromboembolic complication
    

