
# Load required packages
library(caret)
library(mice)
library(miceadds)
library(future) # parrallelisation of MICE 
library(performance)
library(MASS)

# source of training and test set 

source("src/01_raw_to_pre_imputation.r")

##### lasso model 
library(glmnet)

# any variable ending with _YN is logicin train_data
train_data_clean <- train_data %>%
  dplyr::select(- ICU_TRACHEOS_YN, - NEW_BMI, -    RESP_MODE_RAD  ) %>%
  mutate(across(ends_with("_YN"), ~ case_when(
    is.logical(.) ~ .,
    is.numeric(.) ~ . == 1,
    is.character(.) ~ tolower(.) %in% c("y","yes","true","1"),
    TRUE ~ as.logical(.)
  )))


train_data_imputation = futuremice(train_data_clean, m = 5, ncore = availableCores()-1, maxit = 10)

data = complete(train_data_imputation)

#lasso for feature slection of logistic regression model on complete(train_data_imputation
x <- model.matrix(OUT_DEAD_DURING_ICU_YN ~ ., data )[, -1]  # Remove intercept
y <- data$OUT_DEAD_DURING_ICU_YN

# Fit LASSO model
    cv_model <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance")

coef_min <- coef(cv_model, s = "lambda.1se")
nonzero_vars <- rownames(coef_min)[which(coef_min != 0)]
cat("Non-zero variables:\n")
print(nonzero_vars)
