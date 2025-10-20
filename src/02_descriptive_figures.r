# Data analysis for the the UNTIE COVID study publication -descriptive statistics
#library 
library(tidyverse)
library(gtsummary)
library(openxlsx)

# Run 01_raw_to_processed.r to process the data
source("src/01_raw_to_processed.R") 

############### Glucocorticoid treatment therapy initiation during admission to hospital ###############

treatment_admission_figure  = UNITE_2020_corrected |> 
  select(ICU_CORTICO_INTERV_INT, ICU_CORTICO_DURATION_INT, ICU_CORTICO_YN, ICU_CORTICO_INDICATION_RAD, NEW_SUBJECT_ID,
  ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN, ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN, ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN,
  INC_LOS_PRIOR_ADM_INT) |> 
  filter(ICU_CORTICO_YN == 1) |> 
  mutate(cortico_start  = ifelse( ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN == 1, ICU_CORTICO_INTERV_INT,
                                ifelse(ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN == 1, -(INC_LOS_PRIOR_ADM_INT + (INC_LOS_PRIOR_ADM_INT - ICU_CORTICO_INTERV_INT)),
                                0)),
        cortico_end  = ifelse( ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN == 1, ICU_CORTICO_INTERV_INT + ICU_CORTICO_DURATION_INT,
                                ifelse(ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN == 1, -INC_LOS_PRIOR_ADM_INT +ICU_CORTICO_INTERV_INT + ICU_CORTICO_DURATION_INT, 
                                NA)   ) 
        )


################################# Correlation matrix #################################

library(ggcorrplot)

UNITE_2020_labs = UNITE_2020_corrected |>
select(
  #ICU_CORTICO_YN,NEW_PATIENT_ID,
  ICU_CRP_INT,
ICU_WHITE_CELL_INT, 
ICU_NEUTRO_INT,
ICU_LYMPH_DEC,
neutrophil_lymphocyte_ratio,
ICU_FERRITINE_INT,
ICU_DIMERS_INT
) |> 
rename(
  `C-reactive protein` = ICU_CRP_INT, 
    `White cell count`= ICU_WHITE_CELL_INT,
  `Neutrophil count` = ICU_NEUTRO_INT,
  `Lymphocyte count` = ICU_LYMPH_DEC,
  `Neutrophil to Lymphocyte ratio` = neutrophil_lymphocyte_ratio,
  `Ferritin` = ICU_FERRITINE_INT,
  `D-dimer` = ICU_DIMERS_INT
)

# correlation matrix of lab values as graph
 ggcorrplot(cor(UNITE_2020_labs, use = "pairwise.complete.obs"), lab = TRUE, type = "lower")
ggsave( "figures/correlation_matrix.png",dpi = 900, bg = "white")

############## Now with multiple imputation ##############
library(mice)
library(miceadds)
library(ggmice)
UNITE_2020_labs_mice = mice(UNITE_2020_labs,  m=10, n.cores=4)
print("completed imputation")

#correlation analysis after imputation
correlations= miceadds::micombine.cor(UNITE_2020_labs_mice)
attr(correlations, "r_matrix")

#as graph
plot_corr(UNITE_2020_labs)


############## stacked frequency polygons describing admission in relation to steroid therapy ##############

#get data first 
UNITE_2020_steroid_timing  = UNITE_2020_corrected |> 
  mutate( Day_discharge = OUT_HOSP_DURATION_INT - INC_LOS_PRIOR_ADM_INT) |>
  select(ICU_CORTICO_INTERV_INT, ICU_CORTICO_DURATION_INT, ICU_CORTICO_YN, 
  INC_LOS_PRIOR_ADM_INT, OUT_ICU_DURATION_INT, RESP_INTUB_DAYS_AFT_ADM_INT, RESP_DURATION_INV_VENT_INT, RESP_INTUB_DAYS_AFT_ADM_INT,OUT_HOSP_DURATION_OVERALL_INT,Day_discharge
   ) |>
  filter(ICU_CORTICO_YN == 1) |>
  mutate(
    `Hospital stay prior to ICU`  = -`INC_LOS_PRIOR_ADM_INT`, 
    `Initiation of glucocorticoid treatment` = ICU_CORTICO_INTERV_INT, 
    `Discharge from ICU` = OUT_ICU_DURATION_INT,
    `Hospital discharge` = Day_discharge,
    `Invasive ventilation started` = RESP_INTUB_DAYS_AFT_ADM_INT,
    `Extubation from invasive ventilation` = RESP_INTUB_DAYS_AFT_ADM_INT + RESP_DURATION_INV_VENT_INT
    )|> 
    select(
      `Hospital stay prior to ICU`,
      `Initiation of glucocorticoid treatment`,
      `Discharge from ICU`,
      `Hospital discharge`,
      `Invasive ventilation started`,
      `Extubation from invasive ventilation`
    ) |>
 mutate( `Hospital discharge` = ifelse(`Hospital discharge` < 0, -`Hospital discharge`, `Hospital discharge`))  # negative values for hospital discharge
 
 # stacked frequency polygon with ggplot
UNITE_2020_steroid_timing_long = UNITE_2020_steroid_timing |>
  pivot_longer(cols = everything(), names_to = "event", values_to = "days") |>
  filter(!is.na(days))

# Order the ridgeline plot

UNITE_2020_steroid_timing_long$event <- factor(
  UNITE_2020_steroid_timing_long$event,
  levels = rev(c("Hospital stay prior to ICU", "Invasive ventilation started",
  "Initiation of glucocorticoid treatment", "Extubation from invasive ventilation",
   "Discharge from ICU", "Hospital discharge"))
)

#################### Annotate the ridgeline plot with median and IQR ####################

library(dplyr)
library(ggplot2)
library(ggridges)

# 1. Compute summary stats
summary_stats <- UNITE_2020_steroid_timing_long %>%
  group_by(event) %>%
  summarise(
    n        = n(),
    median   = median(days),
    q1       = quantile(days, .25),
    q3       = quantile(days, .75),
    iqr      = q3 - q1,
    interquartiles = paste0( q1,"-", q3)
  ) %>%
  # keep the factor ordering
  mutate(event = factor(event, levels = levels(UNITE_2020_steroid_timing_long$event)))

# 2. Build the ridgeline plot + annotate
x = ggplot(UNITE_2020_steroid_timing_long, aes(x = days, y = event, fill = event,)) +
  geom_density_ridges(
    scale            = 1.2,
    rel_min_height   = 0.01,
    alpha            = 0.8,
    color            = "white"
  ) + 
  # median line + IQR crossbar
  stat_summary(
    aes(y = as.numeric(event)),
    fun        = median,
    geom       = "point",
    color      = "black",
    size       = 2
  ) +
  stat_summary(
    aes(y = as.numeric(event)),
    fun.data   = function(x) {
      d <- quantile(x, c(0.25, 0.75))
      data.frame(y = median(x), ymin = d[1], ymax = d[2])
    },
    geom       = "crossbar",
    width      = 0.1,
    color      = "black"
  ) +
  # text labels pulled from summary_stats
    geom_text(
      data = summary_stats,
      aes(
        x     = q3 + 5, # nudge the label to the right of the ridge
        y     = as.numeric(event) - 0.1, # position just below the ridge
        label = sprintf("%d, %.0f days (%s)", n, median, interquartiles)
      ),
      hjust = 0,
      size  = 6
    ) +
  labs(
    x        = "Days",
    y        = NULL,
    fill     = "Event"
  ) +
  scale_x_continuous(
    limits = c(-25, 100),
    breaks = c(-25, 0, 25, 50, 75),
    labels = c("-25", "ICU admission", "25", "50", "75")
  ) +
  theme_ridges(font_size = 20, grid = TRUE) +
  theme(
    plot.subtitle  = element_text(size = 20, hjust = 0.5),
    axis.title.x   = element_text(size = 20, hjust = 0.5),
    axis.text.y    = element_text(size = 20),
    axis.text.x    = element_text(size = 20),
    legend.position = "none"
  ) +
  scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8)

# Save the plot
ggsave("figures/UNITE_2020_steroid_timing.png", width = 16, height = 10, dpi = 900, bg = "white")
print("figure saved")