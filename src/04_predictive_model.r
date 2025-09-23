#Ths script is for predictive modelling from the UNTIE COVID data set

# First data cleaning and muiltiple imputation must occur. 
source("src/03_imputation_publication.r")

# The methodology for predictive modeling involves several steps
# 1 - univariable analysis 
# 2 - feature selection with backwards step and AIC 
# 3 - final muiltivariable model


## univariable regression

# This is a really helpful package which takes the tbl_uvregression in the gtsummary packages and allows it with mice. 

# remotes::install_github(
#     "jrob95/gtsummary4mice",
#     build_vignettes = FALSE,
#     dependencies    = TRUE
# )

library(gtsummary4mice)

univariable_ICU_mort = tbl_uvregression(UNITE_2020_corrected_imputation,
        method = glm,
        exponentiate = TRUE,
        y = "OUT_DEAD_DURING_ICU_YN"
        )

print("Univariate table with mice produced")

# Extract P<0.2
univariable_ICU_mort <- tbl$table_body %>%
  filter(p.value < 0.2) %>%
  dplyr::select(variable)
sig_vars <- pull(univariable_ICU_mort, variable)

# construct full multivariable formula 
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "))
)

   data = complete(UNITE_2020_corrected_imputation, 1) 

# loop over imputed data sets    

multivariable_backwardsstep = lapply(seq_len(1), function(i) {

    #create a dataset
   data = complete(UNITE_2020_corrected_imputation, i) 

    #run model
   fit <- glm(multi_formula, data = data, family = "binomial")

   #run stepwise selection
   x = stepAIC(fit, direction = "backward", trace  = FALSE )
   print(x)
   print(i)
})

# use mice for final model
glm(formula = OUT_DEAD_DURING_ICU_YN ~ INC_AGE_INT + ICU_CORTICO_YN + 
    INC_HBP_YN + INC_NEURO_YN + INC_ASTHMA_YN + INC_NEOPLASM_YN +
    ICU_RESP_SUPPORT_YN + ICU_WHITE_CELL_INT + ICU_NEUTRO_INT +
    ICU_PLATELETS_INT + ICU_ANTIVIRALS_YN + INF_DURING_ICU_YN +
    ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + ventilation_severity +
    comorbidity_score + neutrophil_lymphocyte_ratio, family = "binomial",
    data = data)