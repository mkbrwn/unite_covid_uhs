#set working directory 
setwd("C:/Users/brownmq/OneDrive - University Hospital Southampton NHS Foundation Trust/Documents/R/UNITE COVID data analysis/UNITE COVID data analysis")

#library 
library(tidyverse)
library(gtsummary)
library(openxlsx)
library(visreg)
library(lme4)

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.r") 

### model variables 
#ICU_CORTICO_YN - eposure - reciept of steroids 
#ICU_CRP_INT - CRP as a contineous varaible
#ICU_CRP_RAD - CRP as a ordinal variable - < 100 mg/L, 100-200 mg/L, > 200 mg/L
#RESP_INV_VENT_YN - invasive ventilation
#OUT_DEAD_DURING_ICU_YN - outcome - death during ICU stay
#NEW_CENTRE_ID - random effect - hospital site

####################################### specify a mixed effects glm model #######################################

#simple model 
model <- glmer(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN + (1|NEW_CENTRE_ID), data = UNITE_2020_corrected, family = binomial)
print(summary(model))

tbl <- tbl_regression(model,    
  conf.int     = TRUE,
  exponentiate = TRUE
)
tbl

#model with interaction of CRP 

model_interaction_crplinear <- glmer(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN * ICU_CRP_INT + (1|NEW_CENTRE_ID), data = UNITE_2020_corrected, family = binomial)
print(summary(model_interaction_crplinear))

tbl <- tbl_regression(model_interaction_crplinear,    
  conf.int     = TRUE,
  exponentiate = TRUE
)
tbl


# It's definitly not linear 

### model diagnostics

# Check for overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rdev <- sum(residuals(model, type = "pearson")^2)
  rdev / rdf
}
overdisp_fun(model_interaction_crpcat)
