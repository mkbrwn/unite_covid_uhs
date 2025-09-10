# Data analysis for the the UNTIE COVID study - inflammatory phenotypes stratification and corticosteroid response
# Dushi analysis 
    # stratification by CRP / exploration of inflammatory markers / Does BMI influence steroid responce 

#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

#library 
library(tidyverse)
library(gtsummary)
library(openxlsx)
library(visreg)

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.r") 

#ggplot of CRP - I would like to smooth this out 
ggplot(UNITE_2020_corrected, aes(x = ICU_CRP_INT)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "white") +
  labs(title = "Distribution of CRP", x = "CRP (mg/L)", y = "Count") +
  theme_minimal()

ggplot(UNITE_2020_corrected, aes(y = ICU_CRP_INT, x = ICU_CRP_RAD)) +
  geom_violin(scale = "width") +
  labs(title = "Distribution of CRP", x = "CRP (mg/L)", y = "Count") +
  theme_minimal()

# table of inflammtory marks/ ICU support/ outcomes by CRP categories 
table_crp_cat_inf_markers <- UNITE_2020_corrected |>
  tbl_summary( include = c("ICU_WHITE_CELL_INT", 
  "ICU_NEUTRO_INT",
  "ICU_SUP_LYMPH",
  "ICU_FERRITINE_INT",
  "ICU_DIMERS_INT",
  "ICU_PRO_CALCIT_DEC",
  "ICU_PLATELETS_INT",
  ), 
  by  = ICU_CRP_RAD) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2), test = all_categorical() ~ 'fisher.test') |> 
  add_overall() |>
  modify_header(label = "**Variable**") |> 
  modify_spanning_header(all_stat_cols() ~ "**C-reactive protein**") |>
  bold_labels()

table_crp_cat_primary_outcome <- UNITE_2020_corrected |>
  tbl_summary( include =  c( OUTCOME_LD, OUT_ICU_DURATION_INT, OUT_HOSP_DURATION_OVERALL_INT ), 
  by  = ICU_CRP_RAD) |> 
  add_n() |> 
  add_p(pvalue_fun = ~style_sigfig(.,digits =2),
   test = list(OUTCOME_LD ~ "fisher.test"),
  test.args = list(OUTCOME_LD ~ list(simulate.p.value = TRUE)) # simulating p value allows computation of p value, alternative would be chi2
) |> 
  add_overall() |>
  modify_header(label = "**Variable**") |>
  modify_spanning_header(all_stat_cols() ~ "**C-reactive protein**") |>
  bold_labels()

table_crp_cat_secondary_outcome <- UNITE_2020_corrected |>
  tbl_summary( include = c(ICU_RRT_DIAL_YN, ICU_INOTROPES_YN, RESP_INTUBATED_YN, RESP_NI_VENT_YN, COAG_THROMBO_COMPLICATION
  ), 
  by  = ICU_CRP_RAD) |>
  add_p(pvalue_fun = ~style_sigfig(.,digits =2)) |> 
  add_overall() |>
  modify_header(label = "**Variable**") |>
  modify_spanning_header(all_stat_cols() ~ "**C-reactive protein**") 

#Save above tables

# Convert to data frame 
  table_crp_cat_inf_markers <- as_tibble(table_crp_cat_inf_markers, col_labels = TRUE)
  table_crp_cat_primary_outcome <- as_tibble(table_crp_cat_primary_outcome, col_labels = TRUE)
  table_crp_cat_secondary_outcome <- as_tibble(table_crp_cat_secondary_outcome, col_labels = TRUE)

  # Create a workbook and add a sheet
  wb_crp <- createWorkbook()
  addWorksheet(wb_crp, "Inflammatory markers")
  addWorksheet(wb_crp, "Primary outcomes")
  addWorksheet(wb_crp, "Secondary outcomes")
  
  writeData(wb_crp, "Inflammatory markers", table_crp_cat_inf_markers)
  writeData(wb_crp, "Primary outcomes", table_crp_cat_primary_outcome)
  writeData(wb_crp, "Secondary outcomes", table_crp_cat_secondary_outcome)
  
  # Save the workbook
  saveWorkbook(wb_crp, "data/processed/statified_crp.xlsx", overwrite = TRUE)

########################################### logistic regression ###########################################

# univariable logistic regression 

    # CRP as continuous 
 simple_regression_CRP = glm(OUTCOME_LD_DEATH ~ ICU_CRP_INT, data =  UNITE_2020_corrected,  family = binomial)
 summary(simple_regression_CRP)

     # plot simple regression model 
    plot(visreg(simple_regression_CRP,"ICU_CRP_INT",scale   = "response",gg = TRUE) +
    labs(x = "ICU CRP (mg/L)", y = "Predicted Probability of Death") +
    theme_minimal())

    #CRP as categorical 
 simple_regression_CRP_cat = glm(OUTCOME_LD_DEATH ~ ICU_CRP_RAD, data =  UNITE_2020_corrected,  family = binomial)
 summary(simple_regression_CRP_cat)

    # plot simple regression model 
    plot(visreg(simple_regression_CRP_cat,"ICU_CRP_RAD",scale   = "response",gg = TRUE) +
    labs(x = "ICU CRP (mg/L)", y = "Predicted Probability of Death") +
    theme_minimal())

######  With adjustment for steroid use.
regression_CRP_adjsteroid = glm(OUTCOME_LD_DEATH ~ ICU_CRP_INT + ICU_CORTICO_YN, data =  UNITE_2020_corrected,  family = binomial)
 summary(regression_CRP_adjsteroid)

regression_CRP_cat_adjsteroid = glm(OUTCOME_LD_DEATH ~ ICU_CRP_RAD + ICU_CORTICO_YN, data =  UNITE_2020_corrected,  family = binomial)
 summary(regression_CRP_cat_adjsteroid)

#polnomial finding best fit 
 poly_regression_CRP = glm(OUTCOME_LD_DEATH ~ poly(ICU_CRP_INT,2) , data =  UNITE_2020_corrected |> filter( !is.na(ICU_CRP_INT)))
 summary(poly_regression_CRP)

 plot(visreg(poly_regression_CRP,"ICU_CRP_INT",scale   = "response",gg = TRUE) +
    labs(x = "ICU CRP (mg/L)", y = "Predicted Probability of Death") +
    theme_minimal())


# further stratification 

#save plots and tables 
