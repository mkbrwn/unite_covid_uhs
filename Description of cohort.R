# Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
# Description of cohort characteristics 

#library 
library(tidyverse)
library(gtsummary)
library(openxlsx)

# Loading data from onedrive -Analysis restricted to 2020 as the 2021 had a 90% prevlance of glucocorticoid use. 
UNITE_2020_corrected = read.csv("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Projects/Dushi/UNITE COVID data analysis/Data/UNITE_2020_corrected.csv", header= TRUE)
UNITE_2021_corrected = read.csv("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Projects/Dushi/UNITE COVID data analysis/Data/UNITE_2021_corrected.csv", header= TRUE)

# PRISMA diagram - filter 
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  filter( !is.na(ICU_CORTICO_YN))

#Variables for synethesis 
UNITE_2020_corrected  = UNITE_2020_corrected |> 
  mutate(INC_BMI_INT = INC_WEIGHT_INT / (INC_HEIGHT_INT/100)^2) |> # BMI
  mutate(centre = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID)) |> #ID for centre
  mutate(patient = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID, "_", NEW_SUBJECT_ID)) |>  #ID for patient
  mutate(OUT_HOSP_DURATION_OVERALL_INT = ifelse( is.na(OUT_HOSP_DURATION_INT), 
                                                 OUT_ICU_DURATION_INT + INC_LOS_PRIOR_ADM_INT, 
                                                 OUT_HOSP_DURATION_INT)) |> #calculating the overall days in hospital (PreICU +ICU) if a patient died in ICU
  mutate(COAG_THROMBO_COMPLICATION = ifelse( COAG_THROMBO_NONE_CB == FALSE, TRUE, FALSE)) # Any thromboembolic complication
        
# summary tables using gtsummary based on the Dushi's protocol 

table_characteristics <- UNITE_2020_corrected |>
  tbl_summary( include = c( INC_AGE_INT, INC_SEX_RAD, INC_HEIGHT_INT, INC_WEIGHT_INT, INC_BMI_INT), 
               by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

table_comorbidity <- UNITE_2020_corrected |>
  tbl_summary( include = c( INC_CARDIAC_DISEASE_YN:INC_HIV_YN
    ), 
               by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

table_preICUaddmission <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_RESP_SUPPORT_YN, ICU_SUPP_TYPE_RAD, ICU_WHITE_CELL_INT:ICU_FERRITINE_INT, ICU_DIMERS_INT, ICU_PLATELETS_INT
  ), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

#Missing medications during ICU 
table_icu_meds <- UNITE_2020_corrected |>
  tbl_summary( include = c(ICU_ANTIVIRALS_YN, INF_ANTIBIO_YN
  ), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

table_corticosteroids <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_CORTICO_YN:ICU_CORTICO_INDICATION_RAD, INF_AT_ADMISSION_YN, INF_DURING_ICU_YN, INF_SEVERITY
  ), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

# PCT during ICU Admission was not collected 
table_duringICUaddmission <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_SUP_WHITE_CELL_INT:ICU_SUP_FERRITINE_DEC, COAG_DIMERS_INT, COAG_PLATELETS_INT ), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

table_primary_outcome <- UNITE_2020_corrected |>
  tbl_summary( include = c( OUTCOME_LD, OUT_ICU_DURATION_INT, OUT_HOSP_DURATION_OVERALL_INT ), 
               by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()
 
table_secondary_outcome <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_RRT_DIAL_YN, ICU_INOTROPES_YN, RESP_INTUBATED_YN, RESP_NI_VENT_YN, COAG_THROMBO_COMPLICATION ), 
               by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

# Save tables as an .csv

  # Convert to data frame 
  table_characteristics <- as_tibble(table_characteristics, col_labels = TRUE)
  table_comorbidity <- as_tibble(table_comorbidity, col_labels = TRUE)
  table_preICUaddmission <- as_tibble(table_preICUaddmission, col_labels = TRUE)
  table_icu_meds <- as_tibble(table_icu_meds, col_labels = TRUE)
  table_corticosteroids <- as_tibble(table_corticosteroids, col_labels = TRUE)
  table_duringICUaddmission <- as_tibble(table_duringICUaddmission, col_labels = TRUE)
  table_primary_outcome <- as_tibble(table_primary_outcome, col_labels = TRUE)
  table_secondary_outcome <- as_tibble(table_secondary_outcome, col_labels = TRUE)
  
  # Create a workbook and add a sheet
  wb <- createWorkbook()
  addWorksheet(wb, "Characteristics")
  addWorksheet(wb, "Comorbidities")
  addWorksheet(wb, "Preadmission")
  addWorksheet(wb, "ICU medications")
  addWorksheet(wb, "Corticosteroids")
  addWorksheet(wb, "During ICU")
  addWorksheet(wb, "Primary outcome")
  addWorksheet(wb, "Secondary outcome")
  
  
  writeData(wb, "Characteristics", table_characteristics)
  writeData(wb, "Comorbidities", table_comorbidity)
  writeData(wb, "Preadmission", table_preICUaddmission)
  writeData(wb, "ICU medications", table_icu_meds)
  writeData(wb, "Corticosteroids", table_corticosteroids)
  writeData(wb, "During ICU", table_duringICUaddmission)
  writeData(wb, "Primary outcome", table_primary_outcome)
  writeData(wb, "Secondary outcome", table_secondary_outcome)
  
  
  # Save the workbook
  saveWorkbook(wb, "summary_tables.xlsx", overwrite = TRUE)


