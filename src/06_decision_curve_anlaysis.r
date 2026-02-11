# This script will run the decision curves analysis for the multivariable model

# load libraries
library(dcurves)

#run previous scipts for models - run time < 5 mins
source("src/05_model_assessment_admission.r")

# place to save figures
fig_dir <- "Figures/decision_curve_analysis"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)



#decision curve analysis for muiltivariable model

muivariable_dca <- dca(OUT_DEAD_DURING_ICU_YN ~ predicted_prob_admission + predicted_prob_admission_without_crp + predicted_prob_cheating + predicted_prob_simple, data = test_data_completed,
  label = list(predicted_prob_admission = "Model from admission variables",
               predicted_prob_admission_without_crp = "Model from admission variables without CRP",
               predicted_prob_cheating = "Cheating model - interventions and complications",
               predicted_prob_simple = "Simple logistic - CRPxCortico"
  )
) %>%
  plot(smooth = TRUE)

ggsave(filename = file.path(fig_dir, "DCA_models.png"), plot = muivariable_dca, width = 12, height = 10, dpi = 150)


muivariable_dca_treatment_avoided <- dca(OUT_DEAD_DURING_ICU_YN ~ predicted_prob_admission + predicted_prob_admission_without_crp + predicted_prob_cheating + predicted_prob_simple, data = test_data_completed,
  label = list(predicted_prob_admission = "Admission Model",
               predicted_prob_admission_without_crp = "Prediction Model without CRP",
               predicted_prob_cheating = "Cheating model",
               predicted_prob_simple = "Simple CRPxCortico model"
  )
) %>%
  net_intervention_avoided() %>%
  plot(smooth = TRUE)

ggsave(filename = file.path(fig_dir, "DCA_models_treatment_avoided.png"), plot = muivariable_dca_treatment_avoided, width = 12, height = 10, dpi = 150)

# Histogram/density of predicted probabilities from models
library(ggplot2)
library(tidyr)
library(dplyr)

# Ensure required columns are present
prob_cols <- c(
  "predicted_prob_admission",
  "predicted_prob_admission_without_crp",
  "predicted_prob_cheating",
  "predicted_prob_simple"
)

missing_cols <- setdiff(prob_cols, names(test_data_completed))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing prediction columns in test_data_completed: %s. Ensure 05_model_assessment_admission.r created them.",
    paste(missing_cols, collapse = ", ")
  ))
}

# Long format for plotting
prob_long <- test_data_completed %>%
  dplyr::select(all_of(c("OUT_DEAD_DURING_ICU_YN", prob_cols))) %>%
  tidyr::pivot_longer(cols = all_of(prob_cols), names_to = "Model", values_to = "Probability") %>%
  mutate(
    Model = recode(Model,
                   predicted_prob_admission = "Admission Model",
                   predicted_prob_admission_without_crp = "Admission (no CRP)",
                   predicted_prob_cheating = "Cheating model",
                   predicted_prob_simple = "Simple CRPxCortico"),
    Outcome = ifelse(as.integer(OUT_DEAD_DURING_ICU_YN) == 1, "Event", "No event")
  )

# Overlayed histogram
p_hist <- ggplot(prob_long, aes(x = Probability, fill = Model)) +
  geom_histogram(alpha = 0.35, position = "identity", bins = 30, color = "white") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Predicted probability distributions (all models)", x = "Predicted probability", y = "Count") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        legend.background = element_rect(fill = "white", colour = NA))

ggsave(filename = file.path(fig_dir, "probability_histograms_overlay.png"), plot = p_hist, width = 10, height = 7, dpi = 150)

# Faceted density by model and colored by outcome
p_dens <- ggplot(prob_long, aes(x = Probability, color = Outcome, fill = Outcome)) +
  geom_density(alpha = 0.25) +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(~ Model, ncol = 2) +
  labs(title = "Predicted probability density by model and outcome", x = "Predicted probability", y = "Density") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA)
  )

ggsave(filename = file.path(fig_dir, "probability_densities_by_model.png"), plot = p_dens, width = 12, height = 10, dpi = 150)

message("Saved probability distribution plots to ", fig_dir)


UNITE_2020_corrected <- UNITE_2020_corrected %>%
  mutate(CRP_cat100 = ifelse(
    ICU_CRP_INT < 100, 1, 0),
    crp_200 = ifelse(
      ICU_CRP_INT < 200, 1, 0)
  )


x <- dca(OUT_DEAD_DURING_ICU_YN ~ ICU_CORTICO_YN*CRP_cat100 + ICU_CORTICO_YN*crp_200,
  data = UNITE_2020_corrected,
  prevalence = 0.66
) %>%
  plot(smooth = TRUE)
