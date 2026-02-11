# This script will run the decision curves analysis for the multivariable model

# load libraries
library(dcurves)

test_data_imputed_scaled = as.data.frame(test_data_imputed_scaled_complete)

# formula 
formula <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
             ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
             + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN + RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
             + neutrophil_lymphocyte_ratio + ICU_CRP_INT + age_group

formula1 <- OUT_DEAD_DURING_ICU_YN ~ NEW_SITE_ID + INC_LOS_PRIOR_ADM_INT + INC_CARDIAC_DISEASE_YN + INC_ASTHMA_YN +
             ICU_RESP_SUPPORT_YN + ICU_PLATELETS_INT + ICU_CARDIAC_THERAPY_YN + ICU_SEPSIS_YN + ICU_STRESS_MYOC_YN + 
             ICU_PNEUMOTHORAX_YN + ICU_ATELECTASIS_YN + ICU_DELIRIUM_YN  + ICU_OBSTRUCTION_YN + ICU_CORTICO_YN + ICU_SEDAT_DURATION_INT +
             + ICU_RRT_DIAL_YN + ICU_INOTROPES_YN +  + RESP_INTUBATED_YN + RESP_HFNC_YN + RESP_INV_VENT_YN + RESP_DURATION_INV_VENT_INT + 
             RESP_ECMO_YN + RESP_PRONE_YN + INF_DURING_ICU_YN + INC_DIABETES1_YN +
             + neutrophil_lymphocyte_ratio + crp_category + age_group

train_data_imputed_scaled_muiltivariable = glm(formula, data = complete(train_data_imputation_scaled),binomial(link = "logit"))
train_data_imputed_scaled_muiltivariable_catCRP = glm(formula1, data = complete(train_data_imputation_scaled),binomial(link = "logit"))

#name muiltivariable model 
muiltivariable_model = train_data_imputed_scaled_muiltivariable
muiltivariable_model_catcrp = train_data_imputed_scaled_muiltivariable_catCRP
    # extract predicted values
    test_data_imputed_scaled$pred_multivariable <- predict(muiltivariable_model, newdata = test_data_imputed_scaled, type = "response")
    test_data_imputed_scaled$pred_multivariablecatcrp <- predict(muiltivariable_model_catcrp, newdata = test_data_imputed_scaled, type = "response")

    #run decision curve analysis
    muivariable_dca = dca(OUT_DEAD_DURING_ICU_YN ~ pred_multivariable + pred_multivariablecatcrp +  ICU_CORTICO_YN, data = test_data_imputed_scaled,
    label = list(pred_multivariable = "Prediction Model",
                pred_multivariablecatcrp = "Prediction Model (CRP)",
                ICU_CORTICO_YN = "Glucocorticoid therapy"
                )
        )  %>%
    plot(smooth = TRUE)

    # save plot
    ggsave("figures/decision_curve_analysis/muiltivariable_dca_plot.png", plot = muivariable_dca, width = 10, height = 6, dpi = 300, bg = "white")

    muivariable_dca_treatment_avoided = dca(OUT_DEAD_DURING_ICU_YN ~ pred_multivariable + pred_multivariablecatcrp + ICU_CORTICO_YN, data = test_data_imputed_scaled,
    label = list(pred_multivariable = "Prediction Model",
                  pred_multivariablecatcrp = "Prediction Model (CRP)",
                  ICU_CORTICO_YN = "Glucocorticoid therapy"
                )
        )  %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE)


    # save plot
    ggsave("figures/decision_curve_analysis/muiltivariable_dca_treatment_avoided_plot.png", plot = muivariable_dca_treatment_avoided, width = 10, height = 6, dpi = 300, bg = "white")


test_data_imputed_scaled

