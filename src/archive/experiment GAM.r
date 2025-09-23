# Install and load mgcv if you haven’t already
library(mgcv)

# Assume your data frame is called df, with:
# – outcome: 0/1 (e.g. death or survival)
# – ICU_CORTICO_YN: factor or 0/1 indicating corticosteroid use
# – ICU_CRP: numeric CRP value (continuous) or converted from categories

# Fit the GAM
UNITE_2020_corrected_gam = UNITE_2020_corrected |> 
  mutate(
    ICU_CORTICO_YN = as.numeric(ICU_CORTICO_YN == TRUE)
  )

gam_model <- gam(
  OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN
         + s(ICU_CRP_INT)                          # baseline CRP effect
         + s(ICU_CRP_INT, by = ICU_CORTICO_YN),    # different smooth if cortico = TRUE
  family = binomial(link = "logit"),
  data   = UNITE_2020_corrected_gam
)

# Summarize the model
summary(gam_model)

# Base plot of smooth terms
plot(gam_model, pages = 1, shade = TRUE, rug = FALSE)

# If you want predicted probabilities over a grid:

# Example: drop rows where CRP is NA
df_clean <- subset(UNITE_2020_corrected_gam, !is.na(ICU_CRP_INT))

# Or coerce factor→numeric (if that’s the issue)
df_clean$ICU_CRP_INT <- as.numeric(as.character(df_clean$ICU_CRP_INT))


newdat <- with(df_clean, expand.grid(
  ICU_CRP_INT        = seq(min(ICU_CRP_INT), max(ICU_CRP_INT), length = 200),
  ICU_CORTICO_YN = c(0, 1)
))
newdat$predicted <- predict(gam_model, newdata = newdat, type = "response")

#plot
ggplot(newdat, aes(ICU_CRP_INT, predicted, color = factor(ICU_CORTICO_YN))) +
  geom_line(linewidth = 1.2) +
  labs(color = "Steroid\nUse", y = "Predicted Probability", x = "CRP")
