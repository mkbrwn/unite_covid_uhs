#Ths script is for predictive modelling from the UNTIE COVID data set

# Fit multilevel model across all imputed datasets using with() from micest data was cleaned with src/03_imputation_muiltilevel_publication.r
# This process takes a while and was saved as a RDS file 

#additional packages for this script
library(tidyverse)
library(gtsummary)
library(gtsummary4mice)
library(lme4)
library(broom.mixed)  # Required for tidy() methods with mixed-effects models
library(miceadds)     # For pooling multilevel models

#load RDS file
UNITE_2020_corrected_multiple_imputation <- readRDS("data/processed/multilevel_imputation_object.rds")

 remotes::install_github(
    "jrob95/gtsummary4mice",
    build_vignettes = FALSE,
    dependencies    = TRUE
    )

# Install broom.mixed if not already installed
if (!require(broom.mixed, quietly = TRUE)) {
  install.packages("broom.mixed")
}

library(gtsummary4mice)
library(broom.mixed)

#univariable analysis - This will take 5 mins 
univariable_multilevel_ICU_mort =tbl_uvregression(UNITE_2020_corrected_multiple_imputation,
        method = glmer,
        exponentiate = TRUE,
        method.args = list(family = binomial),
        y = "OUT_DEAD_DURING_ICU_YN",
        formula = "{y} ~ {x}+ (1 | NEW_SITE_ID)", #specifies the random effects component of the mixed effects mode e.g = each site
        )

print("Univariate table with mice produced")

# Extract P<0.2
univariable_multilevel_ICU_mort <- univariable_multilevel_ICU_mort$table_body %>%
  filter(p.value < 0.2) %>%
  dplyr::select(variable)
sig_vars <- dplyr::pull(univariable_multilevel_ICU_mort, variable)

# construct full multivariable multilevel model formula
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "), " + (1 | NEW_SITE_ID)")
)

print("Initial formula for backwards elimination:")
print(multi_formula)


######### model building with buildmer #############

#load buildmer 
library(buildmer)


convergence_model = buildmer::buildmer(multi_formula, 
                            data = complete(UNITE_2020_corrected_multiple_imputation), family = binomial,
                            args = list(control = list(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))))


x = glmer( OUT_DEAD_DURING_ICU_YN ~ + ( 1| NEW_SITE_ID), data = complete(UNITE_2020_corrected_multiple_imputation), family = binomial)
summary(x)

plot(x)