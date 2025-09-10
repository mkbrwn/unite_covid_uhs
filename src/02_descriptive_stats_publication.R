# Data analysis for the the UNTIE COVID study publication -descriptive statistics
#library 
library(tidyverse)
library(gtsummary)
library(openxlsx)

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.R") 

############### Glucocorticoid treatment therapy description ###############

table_corticosteroid_treatment <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_CORTICO_DURATION_INT:ICU_CORTICO_INDICATION_RAD, ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN,
   ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN, ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

############### Demographics table by steroid treatment ###############

table_comorbidity <- UNITE_2020_corrected |>  tbl_summary(
    include = c(INC_AGE_INT, INC_SEX_RAD, INC_BMI_INT, ICU_ADM_DIAG_RAD,  INC_CARDIAC_DISEASE_YN:INC_HIV_YN),
    by  = ICU_CORTICO_YN,
        statistic = list(
      all_continuous() ~ "{median} [{p25}, {p75}]"
    ),
    digits = list(
      all_continuous() ~ 1
    )
  ) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

# This extra summary for comorbidity score which didn't summarise within the tbl_summary
UNITE_2020_corrected %>%
  mutate(comorbidity_score = as.numeric(comorbidity_score)) %>%
  group_by(ICU_CORTICO_YN) %>%
  summarise(
    median = median(comorbidity_score, na.rm = TRUE),
    p25 = quantile(comorbidity_score, 0.25, na.rm = TRUE),
    p75 = quantile(comorbidity_score, 0.75, na.rm = TRUE),
    n = sum(!is.na(comorbidity_score))
  )
# Kruskal-Wallis test for comorbidity_score by ICU_CORTICO_YN
kruskal_test_result <- kruskal.test(comorbidity_score ~ ICU_CORTICO_YN, data = UNITE_2020_corrected)
print(kruskal_test_result)

# admission labratory values table


table_admission_labs    <- UNITE_2020_corrected |>
  tbl_summary( include = c( ICU_CRP_INT, ICU_WHITE_CELL_INT:ICU_FERRITINE_INT, ICU_DIMERS_INT, ICU_PLATELETS_INT, ICU_LYMPH_DEC, neutrophil_lymphocyte_ratio
  ), 
  by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()

table_outcome <- UNITE_2020_corrected |>
  tbl_summary( include = c( OUTCOME_LD, OUT_DEAD_DURING_ICU_YN, OUT_ICU_DURATION_INT, OUT_HOSP_DURATION_OVERALL_INT, ICU_RRT_DIAL_YN, ICU_INOTROPES_YN, RESP_INV_VENT_YN, RESP_NI_VENT_YN, ventilation_severity), 
               by  = ICU_CORTICO_YN) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Corticoids Received**") |>
  add_overall() |>
  modify_header(label = "**Variable**") |>
  bold_labels()
 
# This extra summary for ventilation severity which didn't summarise within the tbl_summary
UNITE_2020_corrected %>%
  mutate(ventilation_severity = as.numeric(ventilation_severity)) %>%
  summarise(
    median = median(ventilation_severity, na.rm = TRUE),
    p25 = quantile(ventilation_severity, 0.25, na.rm = TRUE),
    p75 = quantile(ventilation_severity, 0.75, na.rm = TRUE),
    n = sum(!is.na(ventilation_severity))
  )
# Kruskal-Wallis test for comorbidity_score by ICU_CORTICO_YN
kruskal_test_result <- kruskal.test(ventilation_severity ~ ICU_CORTICO_YN, data = UNITE_2020_corrected)
print(kruskal_test_result)
