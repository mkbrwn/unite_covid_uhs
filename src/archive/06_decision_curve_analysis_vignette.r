# decision curve analysis using the dca package

#plan
# 1. use complete case 
# 2. uuse imputed data sets 
# 3. figure our the muiltivariable model. 


#run previous script 
source("01_UNITE COVID data analysis/src/01_pre_imputation_univariable_training_data.R")

#load required libraries
library(dcurves)
library(mice)   

# the subsequent analysis follows the vignette from the dcurves package - https://cran.r-project.org/web//packages/dcurves/vignettes/dca.html

#Create bins of CRP
train_data$CRP_bin = cut(train_data$ICU_CRP_INT, breaks = c(-Inf, 100, 200, 300, Inf), labels = c("Low", "Medium", "High", "Very High"))

#complete case
train_data_complete <- train_data[complete.cases(train_data[c("CRP_bin", "OUT_DEAD_DURING_ICU_YN", "ICU_CORTICO_YN")]), ]

# simple model
simple_model = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN , data = train_data_complete, family = binomial(link = "logit"))
simple_model_table = tbl_regression(simple_model, exp = TRUE)

# Create decision curve and ensure it displays in VS Code
dca_result <- dca(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN, data = train_data_complete)
dca_plot <- plot(dca_result, smooth = TRUE)

# Force display in VS Code
print(dca_plot)

# Also save to file so you can view it
ggsave("figures/decision_curve_analysis/dca_simple_corticosteroid_model.png", plot = dca_plot, width = 10, height = 6, dpi = 300, bg = "white")

######## Conditional analysis based on CRP > 100 as this is where mortality increases from baseline 
simple_model_CRP = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN  + CRP_bin, data = train_data_complete, family = binomial(link = "logit"))
simple_model_CRP_table = tbl_regression(simple_model_CRP, exp = TRUE)

#produce a variable indlcate non low risk
train_data_complete = train_data_complete %>%
    mutate(CRP_highrisk = ifelse(CRP_bin != "Low", 1, 0)) #keep as binary

# Create decision curve and ensure it displays in VS Code
dca_CRP_result <- dca(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN + CRP_highrisk, data = train_data_complete)
dca_CRP_plot <- plot(dca_CRP_result, smooth = TRUE)

ggsave("figures/decision_curve_analysis/dca_simple_corticosteroid_crp__model.png", plot = dca_CRP_plot, width = 10, height = 6, dpi = 300, bg = "white")


############################################## With imputation dataset ##############################################


#Create bins of CRP
train_data_imputation_complete = complete(train_data_imputation)
train_data_imputation_complete$CRP_bin = cut(train_data_imputation_complete$ICU_CRP_INT, breaks = c(-Inf, 100, 200, 300, Inf), labels = c("Low", "Medium", "High", "Very High"))


# simple model
simple_imputed_model = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN , data = train_data_imputation_complete, family = binomial)
simple_imputed_model_table = tbl_regression(simple_imputed_model, exp = TRUE)

# Create decision curve and ensure it displays in VS Code
dca_result <- dca(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN, data = train_data_complete)
dca_plot <- plot(dca_result, smooth = TRUE)

# Force display in VS Code
print(dca_plot)

# Also save to file so you can view it
ggsave("figures/decision_curve_analysis/dca_simple_imputed_corticosteroid_model.png", plot = dca_plot, width = 10, height = 6, dpi = 300, bg = "white")

# simple model
simple_imputed_crp_model = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN + CRP_bin, data = train_data_imputation_complete, family = binomial)
simple_imputed_crp_model_table = tbl_regression(simple_imputed_crp_model, exp = TRUE)

train_data_imputation_complete = train_data_imputation_complete %>%
    mutate(CRP_highrisk = ifelse(ICU_CRP_INT >300, 1, 0)) #keep as binary

# Create decision curve and ensure it displays in VS Code
dca_result <- dca(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN + CRP_highrisk, data = train_data_complete)
dca_plot <- plot(dca_result, smooth = TRUE)

# Force display in VS Code
print(dca_plot)

# Also save to file so you can view it
ggsave("figures/decision_curve_analysis/dca_simple_imputed_crp_corticosteroid_model.png", plot = dca_plot, width = 10, height = 6, dpi = 300, bg = "white")



train_data_imputation_complete$CRP_bin = ifelse(train_data_imputation_complete$ICU_CRP_INT < 150, 0,1 )

# simple model
simple_imputed_crp_model = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*CRP_bin, data = train_data_imputation_complete, family = binomial)
simple_imputed_crp_model_table = tbl_regression(simple_imputed_crp_model, exp = TRUE)
simple_imputed_crp_model_table

