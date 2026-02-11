# prospensity score modeling from imputed data 

#load script 
source("src/archive/Preprocessed_propensity_score.r")

#prodce age categories
UNITE_2020_corrected = UNITE_2020_corrected %>%
  mutate(INC_AGE_INT= cut(
    INC_AGE_INT,
    breaks = c(-Inf, 40, 50, 60, 70, 80, Inf),
    labels = c("<40", "40-50", "50-60", "60-70", "70-80", ">80")
  ))

#remove variables unsuitable for imputation
  UNITE_2020_corrected = UNITE_2020_corrected %>%
   dplyr::select(-ICU_CRP_CATEGORY, -RESP_MODE_RAD)



# to install gtsummary4mice if not already installed
remotes::install_github(
    "jrob95/gtsummary4mice",
    build_vignettes = FALSE,
    dependencies    = TRUE
    )
#load packages 
library(MatchIt)
library(mice)
library(gtsummary4mice)
library(MASS)
library(cobalt)
library(sandwich)
library(lmtest)
library(marginaleffects)

#contruct prodictive model from imputed data 
imputed_data = complete(futuremice(UNITE_2020_corrected, method = "rf", m = 5))

# Check for columns with NA in imputed_data
na_cols <- names(imputed_data)[sapply(imputed_data, function(x) any(is.na(x)))]
if(length(na_cols) > 0) {
  cat("Columns with NA in imputed_data:\n")
  print(na_cols)
} else {
  cat("No columns with NA in imputed_data\n")
}

#remove na_cols from imputed_data
imputed_data = imputed_data %>%
    dplyr::select( -ICU_TRACHEOS_YN, -NEW_BMI)

# produce univariable summary table for corticoids
uni_estimate_corticoids = tbl_uvregression(
    imputed_data,
    method = glm,
    y = ICU_CORTICO_YN,
    exp = TRUE
) 

#extract varaibles with p<0.2
uni_estimate_corticoids_0.2 <- uni_estimate_corticoids$table_body %>%
  filter(p.value < 0.2) %>%
  pull(variable)

#remove OUT_DEAD_DURING_ICU_YN from uni_estimate_corticoids_0.2 if present
uni_estimate_corticoids_0.2 <- uni_estimate_corticoids_0.2[uni_estimate_corticoids_0.2 != "OUT_DEAD_DURING_ICU_YN"] 

# Fit full model with selected variables
full_formula <- as.formula(paste("ICU_CORTICO_YN ~", paste(uni_estimate_corticoids_0.2, collapse = " + ")))

full_model <- glm(full_formula, data = imputed_data, family = binomial())

# Backwards step elimination using AIC
step_model <- stepAIC(full_model, direction = "backward", trace = FALSE)

# Use step_model for propensity score modeling
ps_match <- matchit(formula(step_model), 
                   data = imputed_data, 
                   method = "nearest", 
                   distance = "glm",
                   caliper = 0.2)

# Summary of matching
summary(ps_match)

################################# Propensity score diagnostics ############################################
# Love plot (standardized mean differences)
png("Figures/propensity score/love_plot.png", width = 800, height = 600)
love.plot(ps_match, 
          threshold = 0.1,
          title = "Covariate Balance Before and After Matching")
dev.off()

# Assess balance across all variables used in ps_match
cat("\n=== VARIABLES USED IN PROPENSITY SCORE MODEL ===\n")

bal.tab(ps_match)

# Comprehensive balance plots for all matched variables
pdf("Figures/propensity score/balance_assessment_all_variables.pdf", width = 16, height = 12)

################## assess ICU_CORTICO_YN on OUT_DEAD_DURING_ICU_YN from propensity score ########################

matched_data <- match.data(ps_match)

model <- glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN, data = matched_data, family = binomial)
tbl_regression(model, exp = TRUE) %>% print()

##### form categories of CRP 
matched_data <-matched_data %>%  mutate(
    crp_category100 = cut(
      ICU_CRP_INT,
      breaks = c(-Inf, 100, 200, 300, 400, Inf),
      labels = c("<100", "100-200", "200-300", "300-400", ">400")
    )) 
matched_data$crp_category100 <- relevel(matched_data$crp_category100, ref = "<100")

model_cat = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*crp_category100, data = matched_data, family = binomial)


imputed_data$predicted_corticoids = predict(glm(formula(step_model), data = imputed_data, family = binomial),
                                       type = "response")
model_propensity_cortico = glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*crp_category100 + predicted_corticoids,
                                data = imputed_data,
                                family = binomial)
summary(model_propensity_cortico)

model_propensity_cortico %>% tbl_regression(exp = TRUE) %>% print()

# Cluster by matched pair ID (from MatchIt object)
rlibrary(sandwich)
library(lmtest)
library(broom)
library(gtsummary)

# Fit model
model_cat <- glm(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN * crp_category100,
                 data = matched_data,
                 family = binomial)

# Cluster-robust SE
robust_se <- vcovCL(model_cat, cluster = matched_data$subclass)

# Tidy function that accepts extra arguments
tidy_cluster <- function(x, ...) {
  broom::tidy(x, conf.int = TRUE, exponentiate = TRUE, vcov = robust_se)
}

# gtsummary table
tbl <- tbl_regression(
  model_cat,
  exponentiate = TRUE,
  tidy_fun = tidy_cluster
)

tbl
#Get marginal effects or risk ratios by CRP group
avg_comparisons(model_cat,
                variables = "ICU_CORTICO_YN",
                by = "crp_category100")



plot_cap(model_cat, condition = "crp_category100",
         variables = "ICU_CORTICO_YN") +
  labs(y = "Predicted probability of ICU death",
       x = "Glucocorticoid use")



avg_comparisons(
  model_cat,
  hypothesis = "reference",
  variable = "ICU_CORTICO_YN",
  by =  "crp_category100")
