# run DCA with different cut off for CRP 
#source("src/01_raw_to_pre_imputation.R")
source("src/archive/Preprocessed_marginal.r")

# libraries needed for DCA
library(rmda)
library(dplyr)
library(dcurves)
library(purrr)
library(broom)
library(ggplot2)
library(marginaleffects)


library(mgcv) # for GAM models


# mice for data imputation
library(mice)
library(future)

imputed_data = complete(futuremice(UNITE_2020_corrected, meth = "rf", m = 5))

#ventilated patients only
#imputed_data = imputed_data %>% filter( RESP_INTUBATED_YN == 1 )

imputed_data = imputed_data %>%
  mutate(crp_category100 = case_when(
    ICU_CRP_INT < 100 ~ "<100",
    ICU_CRP_INT >= 100 & ICU_CRP_INT <= 200 ~ "100-200",
    ICU_CRP_INT >= 200 & ICU_CRP_INT <= 300 ~ "200-300",
    ICU_CRP_INT >= 300 & ICU_CRP_INT <= 400 ~ "300-400",
    ICU_CRP_INT > 400 ~ ">400"
  ),
    crp_category50 = case_when(
    ICU_CRP_INT < 50 ~ "<50",
    ICU_CRP_INT >= 50 & ICU_CRP_INT < 100 ~ "50-100",
    ICU_CRP_INT >= 100 & ICU_CRP_INT < 150 ~ "100-150",
    ICU_CRP_INT >= 150 & ICU_CRP_INT < 200 ~ "150-200",
    ICU_CRP_INT >= 200 & ICU_CRP_INT < 250 ~ "200-250",
    ICU_CRP_INT >= 250 & ICU_CRP_INT <= 300 ~ "250-300",
    ICU_CRP_INT > 300 ~ ">300",
  ))

# make baseline categories and reference levels for crp_category50 and crp_category100
imputed_data$crp_category50 <- factor(
  imputed_data$crp_category50,
  levels = c("<50", "50-100", "100-150", "150-200", "200-250", "250-300", ">300"))
imputed_data$crp_category50 <- relevel(imputed_data$crp_category50, ref = "<50")    

imputed_data$crp_category100 <- factor(
  imputed_data$crp_category100,
  levels = c("<100", "100-200", "200-300", "300-400", ">400"))
imputed_data$crp_category100 <- relevel(imputed_data$crp_category100, ref = "<100")


# Ensure factors and reference levels
imputed_data$crp_category100 <- factor(
  imputed_data$crp_category100,
    levels = c("<100", "100-200", "200-300", "300-400", ">400"))

imputed_data$crp_category100 <- stats::relevel(imputed_data$crp_category100, ref = "<100")
if (!is.factor(imputed_data$ICU_CORTICO_YN)) {
  imputed_data$ICU_CORTICO_YN <- factor(imputed_data$ICU_CORTICO_YN)
}

# construction model 
glm_model_cat = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*crp_category100 + INC_AGE_INT + INC_SEX_RAD + comorbidity_score, data = imputed_data ,binomial(link = "logit"))
glm_model_cont = glm_model_cont = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*ICU_CRP_INT + INC_AGE_INT + INC_SEX_RAD + comorbidity_score +  ventilation_severity + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN, data = imputed_data ,binomial(link = "logit"))

tbl_regression(glm_model_cat, exp = TRUE) %>% print()

# Comparison of recept of corticosteroids at different CRP levels
avg_comparisons(
  glm_model_cat,
  variable = "ICU_CORTICO_YN",
   by =  "crp_category100",
   comparison = "difference"
  )

avg_comparisons(
  glm_model_cat,
  hypothesis = "reference",
  variable = "ICU_CORTICO_YN",
  by =  "crp_category100")

# Plot the predicted probabilities for different CRP categories and corticosteroid use

# Use explicit namespace, assign, print, and save to file to ensure rendering
p <- marginaleffects::plot_comparisons(
  glm_model_cat,
  variables = "ICU_CORTICO_YN",
  condition = "crp_category100",
  type = "response",
  comparison = "difference"
) + ggplot2::labs(y = "Conditional risk difference") + ggplot2::theme_minimal( )

print(p)

# Persist to disk in case IDE/terminal doesn't auto-show plots

ggplot2::ggsave(
  filename = "Figures/marginal_models/steroid_effect_by_crp.png",
  plot = p,
  width = 7, height = 5, dpi = 300,
  bg = "white"
)
################################################################## GAMM MODEL #######################################################
  
  gamm_model <- gamm(
  OUT_DEAD_DURING_ICU_YN ~ s(ICU_CRP_INT, by = ICU_CORTICO_YN)+ INC_AGE_INT + INC_SEX_RAD + comorbidity_score +  ventilation_severity + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN,
    # random intercept for patients
  family = "binomial",           # change to binomial if OutcomeVar is binary
  data = imputed_data
)
summary(gamm_model$gam)

# Plot the smooth terms
plot(gamm_model$gam, 
     pages = 1,        # all plots on one page
     rug = TRUE,       # adds rug marks for data points
     seWithMean = TRUE, # include standard error
     shade = TRUE)     # shaded confidence intervals