#Comparison of imputed models 

#model 1 
model_imputed = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN,  family = binomial))


 # model 2 with the addition of CRP 
model_imputed_crp = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ICU_CRP_INT,  family = binomial))


model_imputed_crp_rad = with(data=UNITE_2020_corrected_mice_minimum,
 expr=glm(OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ICU_CRP_RAD,  family = binomial))



 # comparison of models 

 #Walds test 
 print(summary(D1( model_imputed_crp, model_imputed)))

#likelihood test
 print(summary(D3( model_imputed_crp, model_imputed_crp_rad)))



##################################################################  visualising spline ##################################################################

library(mice)
library(splines)

model_spline <- with(
  data = UNITE_2020_corrected_mice_minimum,
  expr = glm(
    OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ns(ICU_CRP_INT, df = 3),
    family = binomial
  )
)

# Define a fine grid over observed CRP
crp_grid <- seq(
  min(UNITE_2020_corrected_mice_minimum$ICU_CRP_INT, na.rm = TRUE),
  max(UNITE_2020_corrected_mice_minimum$ICU_CRP_INT, na.rm = TRUE),
  length.out = 100
)

# Expand over steroid use
newdata <- expand.grid(
  ICU_CRP_INT      = crp_grid,
  ICU_CORTICO_YN   = c(FALSE, TRUE)
)

#### 
library(mice)

imputed_list <- complete(
  UNITE_2020_corrected_mice_minimum,
  action = "all"
)

# pull out CRP from each imputed dataset
all_crp_values <- unlist(
  lapply(imputed_list, function(df) df$ICU_CRP_INT)
)

# compute finite range
crp_range <- range(all_crp_values, na.rm = TRUE)
crp_range

crp_grid <- seq(
  from       = crp_range[1],
  to         = crp_range[2],
  length.out = 100
)

newdata <- expand.grid(
  ICU_CRP_INT    = crp_grid,
  ICU_CORTICO_YN = c(FALSE, TRUE)
)


# Extract all completed datasets
imputed_list <- complete(UNITE_2020_corrected_mice_minimum, action = "all")

# Generate predicted probabilities per imputation
pred_mat <- sapply(imputed_list, function(d){
  fit  <- glm(
    OUTCOME_LD_DEATH ~ ICU_CORTICO_YN + ns(ICU_CRP_INT, df = 3),
    data = d,
    family = binomial
  )
  predict(fit, newdata = newdata, type = "response")
})

# Pool: mean and 95% “imputation” intervals
newdata$pred_mean <- rowMeans(pred_mat)
newdata$pred_lwr  <- apply(pred_mat, 1, quantile, probs = 0.025)
newdata$pred_upr  <- apply(pred_mat, 1, quantile, probs = 0.975)

library(ggplot2)

ggplot(newdata, aes(x = ICU_CRP_INT, y = pred_mean, color = ICU_CORTICO_YN)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = pred_lwr, ymax = pred_upr, fill = ICU_CORTICO_YN),
    alpha = 0.2, color = NA
  ) +
  scale_color_manual(
    values = c("FALSE" = "#2C3E50", "TRUE" = "#E74C3C"),
    labels = c("No Steroids", "Steroids")
  ) +
  scale_fill_manual(
    values = c("FALSE" = "#2C3E50", "TRUE" = "#E74C3C"),
    labels = c("No Steroids", "Steroids")
  ) +
  labs(
    x     = "CRP (mg/L)",
    y     = "Predicted Probability of Death",
    color = "ICU Corticosteroids",
    fill  = "ICU Corticosteroids"
  ) +
  theme_minimal(base_size = 14)