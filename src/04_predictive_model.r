#Ths script is for predictive modelling from the UNTIE COVID data set

# First data cleaning and muiltiple imputation must occur. 
source("src/03_imputation_publication.r")

# The methodology for predictive modeling involves several steps
# 1 - univariable analysis 
# 2 - feature selection with backwards step and AIC 
# 3 - final muiltivariable model


## univariable regression

# This is a really helpful package which takes the tbl_uvregression in the gtsummary packages and allows it with mice. 

 remotes::install_github(
    "jrob95/gtsummary4mice",
    build_vignettes = FALSE,
    dependencies    = TRUE,
    force = TRUE
)

library(gtsummary4mice)

univariable_ICU_mort =tbl_uvregression(train_data_imputation,
        method = glm,
        exponentiate = TRUE,
        y = "OUT_DEAD_DURING_ICU_YN"
        )
print("Univariate table with mice produced")

# Extract P<0.2
univariable_ICU_mort <- univariable_ICU_mort$table_body %>%
  filter(p.value < 0.2) %>%
  dplyr::select(variable)
sig_vars <- dplyr::pull(univariable_ICU_mort, variable)

# construct full multivariable formula 
multi_formula <- as.formula(
  paste0("OUT_DEAD_DURING_ICU_YN ~ ", paste(sig_vars, collapse = " + "))
)

   data = complete(train_data_imputation, 1) 

# loop over imputed data sets    

multivariable_backwardsstep = lapply(seq_len(5), function(i) {

    #create a dataset
   data = complete(train_data_imputation, i) 

    #run model
   fit <- glm(multi_formula, data = data, family = "binomial")

   #run stepwise selection
   x = step(fit, direction = "backward", trace  = FALSE )
   print(x)
   print(i)
})



# complete data 
data = within(data, {
  CRP100 = cut(ICU_CRP_INT, breaks = c(-Inf, 100, 200, Inf), labels = c("Low", "Medium", "High"))
})

data = complete(train_data_imputation) 

# Define the final multivariable model formula
final_formula <- OUT_DEAD_DURING_ICU_YN ~ INC_AGE_INT + ICU_CORTICO_YN* + 
    INC_CARDIAC_DISEASE_YN + INC_LIVER_DISEASE_YN + INC_HBP_YN +
    INC_ASTHMA_YN + INC_KIDNEY_DISEASE_YN + INC_LOS_PRIOR_ADM_INT +
    ICU_RESP_SUPPORT_YN + ICU_WHITE_CELL_INT + ICU_NEUTRO_INT +
    ICU_PLATELETS_INT + INF_DURING_ICU_YN + ICU_CLIN_TRIAL_YN +
    ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + ventilation_severity +
    INC_DIABETES1_YN + comorbidity_score + neutrophil_lymphocyte_ratio
    
# multivariable model

muiltivariable_model <- glm(
  formula = final_formula,
  family = "binomial",
  data = data
) 

muiltivariable_model_gtsummary = muiltivariable_model |> 
  tbl_regression(
    exponentiate = TRUE
  )