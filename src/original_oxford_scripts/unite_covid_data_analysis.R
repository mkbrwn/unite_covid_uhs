library(tidyverse)
library(dplyr)
library(lubridate)
library(MatchIt)
library(MASS)
library(stats)

## propensity matching for covid unite

#load data
covid_data<-read.csv("UNITE_2020_corrected.csv")


# Select only the columns needed for the analysis

reduced_c_df <- covid_data %>% dplyr::select(NEW_COUNTRY_ID,NEW_SUBJECT_ID, NEW_SITE_ID, INC_AGE_INT, INC_SEX_RAD,
                                             INC_CARDIAC_DISEASE_YN, INC_LIVER_DISEASE_YN, INC_HBP_YN, INC_NEURO_YN,
                                             INC_NEOPLASM_YN, INC_PULMO_DISEASE_YN, INC_DIABETES_YN, INC_ASTHMA_YN,
                                             INC_KIDNEY_DISEASE_YN, INC_HIV_YN,INF_SEVERITY,
                                             INC_IMMUNOSUPPR_YN, NEW_BMI, ICU_CORTICO_INTERV_INT,
                                             ICU_KIDNEY_INJ_YN, ICU_RRT_DIAL_YN,ICU_INOTROPES_YN,
                                             RESP_INTUBATED_YN, RESP_INTUBATED_ICU_STAY_YN, RESP_NI_VENT_YN, RESP_HFNC_YN,
                                             RESP_INV_VENT_YN, RESP_ECMO_YN,RESP_PRONE_YN, RESP_NEUROM_BLOC_YN,
                                             INF_ANTIBIO_YN, INF_AT_ADMISSION_YN, INF_DURING_ICU_YN,INF_SEVERITY,INF_MDR_PATHO_YN,
                                             OUTCOME_LD, OUT_ICU_DURATION_INT,RESP_DURATION_INV_VENT_INT,
                                             ICU_CORTICO_YN,ICU_CORTICO_DURATION_INT,ICU_CORTICO_INDICATION_RAD, ICU_CORTICO_INTERV_INT, 
                                             ICU_ANTIMALARIAL_YN, ICU_ANTIVIRALS_YN, ICU_ANTIVIRALS_RAD, ICU_OTHER_ANTIVIRALS_RAD, 
                                             OUT_ICU_DURATION_INT, OUT_HOSP_DURATION_INT, OUT_DEATH_DURING_ICU_YN, OUT_DEAD_DURING_ICU_YN,OUTCOME_LD)

# create a column that combines site and country code to have single identifier for a site
reduced_c_df<- reduced_c_df %>% mutate(centre = paste(NEW_COUNTRY_ID, "_", NEW_SITE_ID))



# Different parameters needed for the analysis:
# age, sex, comorbidity, BMI
# steroids y/n, steroid duration
# antivirals
# ventilated y/n, ventilation duration, ventilation severity
# icu LOS, hospital LOS, death yn
# infection on ICU adm, infection in ICU, antibiotics yn, infection severity, infection MDR

# Data cleaning to correct NAs and other challenges
# OUTCOME
## Correcting NA's in DEAD_IN_ICU to "FALSE" based on OUTCOME_LD not Death/NA
reduced_c_df$OUT_DEAD_DURING_ICU_YN[reduced_c_df$OUTCOME_LD == "Still in ICU"] <- "FALSE"
reduced_c_df$OUT_DEAD_DURING_ICU_YN[reduced_c_df$OUTCOME_LD == "Hospitalized"] <- "FALSE"
reduced_c_df$OUT_DEAD_DURING_ICU_YN[reduced_c_df$OUTCOME_LD == "Transfer to other facility"] <- "FALSE"
reduced_c_df$OUT_DEAD_DURING_ICU_YN[reduced_c_df$OUTCOME_LD == "Discharged alive"] <- "FALSE"
reduced_c_df$OUT_DEAD_DURING_ICU_YN[reduced_c_df$OUTCOME_LD == "Palliative discharge"] <- "FALSE"

## Antivirals - correcting for NAs assuming that NA is FALSE as they would have entered TRUE if they had received them
reduced_c_df$ICU_ANTIVIRALS_YN[!is.na(reduced_c_df$ICU_ANTIVIRALS_RAD)] <- TRUE
reduced_c_df$ICU_ANTIVIRALS_YN[!is.na(reduced_c_df$ICU_OTHER_ANTIVIRALS_RAD)] <- TRUE


# replacing NA with FALSE for ventilation data - again assuming TRUE would hav ebeen entered if it had been provided
reduced_c_df<- reduced_c_df %>% replace_na(list(RESP_PRONE_YN = FALSE, RESP_ECMO_YN = FALSE, RESP_NEUROM_BLOC_YN = FALSE, ICU_ANTIMALARIAL_YN=FALSE,
                                                ICU_ANTIVIRALS_YN = FALSE))


#ventilation_severity - score the ventilation severity using different levels of ventilatory support
reduced_c_df <-reduced_c_df %>% mutate(ventilation_severity = case_when(RESP_ECMO_YN == TRUE ~3, RESP_PRONE_YN==TRUE ~3,
                                                                        RESP_NEUROM_BLOC_YN ==TRUE & RESP_PRONE_YN ==FALSE ~2,
                                                                        RESP_PRONE_YN==FALSE & RESP_NEUROM_BLOC_YN == FALSE ~1))

# comorbidity score - replace NA with FALSE to allow for scoring
reduced_c_df<- reduced_c_df %>% replace_na(list(INC_CARDIAC_DISEASE_YN = FALSE, INC_PULMO_DISEASE_YN = FALSE,
                                                INC_ASTHMA_YN = FALSE, INC_DIABETES_YN = FALSE, INC_NEURO_YN = FALSE,
                                                INC_HBP_YN = FALSE, INC_NEOPLASM_YN = FALSE, INC_LIVER_DISEASE_YN = FALSE,
                                                INC_IMMUNOSUPPR_YN = FALSE, INC_HIV_YN = FALSE, INC_KIDNEY_DISEASE_YN= FALSE))

reduced_c_df<- reduced_c_df %>% mutate(cardiac_d = if_else(INC_CARDIAC_DISEASE_YN ==TRUE, 1,0),
                                       liver_d = if_else(INC_LIVER_DISEASE_YN ==TRUE, 1,0),
                                       neuro_d = if_else(INC_NEURO_YN ==TRUE, 1,0),
                                       diabetes_d = if_else(INC_DIABETES_YN ==TRUE, 1,0),
                                       kidney_d = if_else(INC_KIDNEY_DISEASE_YN ==TRUE, 1,0),
                                       htn_d = if_else(INC_HBP_YN ==TRUE, 1,0),
                                       asthma_d = if_else(INC_ASTHMA_YN ==TRUE, 1,0),
                                       resp_d = if_else(INC_CARDIAC_DISEASE_YN ==TRUE, 1,0),
                                       immunosup_d = if_else(INC_IMMUNOSUPPR_YN ==TRUE, 1,0),
                                       hiv_d = if_else(INC_HIV_YN ==TRUE, 1,0))

# compute an additive comorbitdity score where each comorbidity has equal weight:
reduced_c_df<- reduced_c_df %>% mutate(comorb = cardiac_d+liver_d+neuro_d+diabetes_d+kidney_d+htn_d+asthma_d+resp_d+immunosup_d+hiv_d)

# remove extraneous columns
reduced_c_df<- reduced_c_df %>% dplyr::select(-INC_CARDIAC_DISEASE_YN, - INC_LIVER_DISEASE_YN, - INC_NEURO_YN,
                                       - INC_DIABETES_YN, -INC_KIDNEY_DISEASE_YN, - INC_HBP_YN,
                                       - INC_ASTHMA_YN, - INC_CARDIAC_DISEASE_YN, - INC_IMMUNOSUPPR_YN,
                                       - INC_HIV_YN, - cardiac_d, - liver_d, - neuro_d, - diabetes_d,
                                       - kidney_d, -htn_d, -asthma_d, -resp_d, - immunosup_d,-hiv_d)
reduced_c_df<- reduced_c_df %>% dplyr::select(-RESP_ECMO_YN, - RESP_NEUROM_BLOC_YN, -RESP_PRONE_YN,
                                              -NEW_COUNTRY_ID, -NEW_SITE_ID)


# slightly redundant but ensuring the correct columns are in the dataset only
reduced_c_df<- reduced_c_df %>% dplyr::select(centre, ventilation_severity, comorb, NEW_SUBJECT_ID,
                                              INC_AGE_INT, INC_SEX_RAD, INF_SEVERITY, NEW_BMI,
                                              RESP_INTUBATED_YN, RESP_INTUBATED_ICU_STAY_YN,
                                              RESP_INV_VENT_YN, INF_ANTIBIO_YN, INF_AT_ADMISSION_YN,
                                              INF_DURING_ICU_YN,ICU_CORTICO_INTERV_INT,ICU_RRT_DIAL_YN,ICU_INOTROPES_YN,
                                              INF_MDR_PATHO_YN, OUTCOME_LD, OUT_ICU_DURATION_INT,
                                              RESP_DURATION_INV_VENT_INT, ICU_CORTICO_YN, ICU_CORTICO_DURATION_INT,
                                              ICU_CORTICO_INDICATION_RAD, ICU_ANTIMALARIAL_YN, ICU_ANTIVIRALS_YN, ICU_ANTIVIRALS_RAD,
                                              ICU_OTHER_ANTIVIRALS_RAD,OUT_HOSP_DURATION_INT, OUT_DEATH_DURING_ICU_YN, OUT_DEAD_DURING_ICU_YN)
                                              
saveRDS(reduced_c_df, "reduced_0406.rds")
reduced_c_df<- readRDS("reduced_0406.rds")


# consider replacing NA in steroids or antivirals with false as well - assuming that anyone who had steroids would have been marked
#*

## select only the ventilated patients
# ventilated with no infection on admission
filter_vent_no_inf <- reduced_c_df %>% filter(RESP_INV_VENT_YN == TRUE) %>% filter(INF_AT_ADMISSION_YN ==FALSE)

# ventilated with potential infection at admission
filter_vent_allow_infection<-reduced_c_df %>% filter(RESP_INV_VENT_YN == TRUE)

# ventilated with definitively postitive for additional infection on admission
filter_vent_definiteinf<- reduced_c_df %>% filter(RESP_INV_VENT_YN==TRUE) %>% filter(INF_AT_ADMISSION_YN== TRUE)
saveRDS(filter_vent_no_inf, "filter_vent_no_inf.rds")
saveRDS(filter_vent_allow_infection, "filter_vent_allow_infection.rds")




#missing data - most important variables are
# ventilation severity, ventilation duration, age, sex, length of stay hospital and legnth of stay ICU, death
# just a check for missing values
missing_df<- filter_vent_no_inf %>% group_by(centre) %>% 
  summarise(na_severity = sum(is.na(ventilation_severity)) , na_duration = sum(is.na(RESP_DURATION_INV_VENT_INT)),
            na_age = sum(is.na(INC_AGE_INT)), na_loshosp = sum(is.na(OUT_HOSP_DURATION_INT)),
            na_losicu = sum(is.na(OUT_ICU_DURATION_INT)), n=max(row_number()))


missing_addition<- filter_vent_no_inf %>% filter(is.na(RESP_DURATION_INV_VENT_INT)) %>% group_by(centre) %>% 
  summarise(na_mean_severity = mean(ventilation_severity), na_mean_comorbidity = mean(comorb),
            na_mean_age= mean(INC_AGE_INT), na_mean_losicu= mean(OUT_ICU_DURATION_INT),
            na_mean_loshosp= mean(OUT_HOSP_DURATION_INT))

missing_add2<- filter_vent_no_inf%>% group_by(centre) %>% 
  summarise(mean_severity = mean(ventilation_severity), mean_comorbidity = mean(comorb),
            mean_age= mean(INC_AGE_INT), mean_losicu= mean(OUT_ICU_DURATION_INT, na.rm=TRUE),
            mean_loshosp= mean(OUT_HOSP_DURATION_INT, na.rm=TRUE))

missing_df<- missing_df %>% left_join(missing_addition, by = "centre")
missing_df<- missing_df %>% left_join(missing_add2, by = "centre")

missing_df<- missing_df %>% mutate(percent_missing = na_duration/n)
lots_missing<- missing_df %>% filter(percent_missing>0.25)

missingtests.lm <- lm(mean_losicu ~ na_mean_losicu, data = lots_missing)
summary(missingtests.lm)



#look at missing things overall not per centre
overall_missing<- filter_vent_no_inf %>%
  summarise(na_severity = sum(is.na(ventilation_severity)) , na_duration = sum(is.na(RESP_DURATION_INV_VENT_INT)),
            na_age = sum(is.na(INC_AGE_INT)), na_loshosp = sum(is.na(OUT_HOSP_DURATION_INT)),
            na_losicu = sum(is.na(OUT_ICU_DURATION_INT)), n=max(row_number()))
overall_miss_add<-filter_vent_no_inf %>% filter(is.na(RESP_DURATION_INV_VENT_INT)) %>% 
  summarise(na_mean_severity = mean(ventilation_severity), na_mean_comorbidity = mean(comorb),
            na_mean_age= mean(INC_AGE_INT))
overall_add2<- filter_vent_no_inf%>% 
  summarise(mean_severity = mean(ventilation_severity), mean_comorbidity = mean(comorb),
            mean_age= mean(INC_AGE_INT, na.rm = TRUE))
overall_missing<- overall_missing %>% cbind(overall_miss_add)
overall_missing<- overall_missing %>% cbind(overall_add2)

# 
# missingness is centre dependent
# the ones where the duration of ventilation is missing is centre dependent and they are at random more or less sick.




# check which centres have significant missing data
library(naniar)
gg_miss_var(filter_vent_no_inf, facet = centre)

#
##################################################################
# Analysis proceeding with ventilated and definitively no infection at admission

ventilation_ix<- filter_vent_no_inf
ventilation_ix<- ventilation_ix%>% mutate(cort_steroids = if_else(ICU_CORTICO_YN ==TRUE, 1,0), 
                                          inf_during_icu = if_else(INF_DURING_ICU_YN ==TRUE, 1,0)) 
ventilation_ix<- ventilation_ix %>% dplyr::select(NEW_SUBJECT_ID, ventilation_severity, 
                                                  comorb, INF_DURING_ICU_YN, INF_AT_ADMISSION_YN, INC_AGE_INT, INC_SEX_RAD)%>% drop_na()

ventilation_ix %>% group_by(ventilation_severity) %>% summarise(mean_inf_status = median(INF_DURING_ICU_YN), sd_inf = sd(INF_DURING_ICU_YN))
ventilation_ix$INF_DURING_ICU_YN<-factor(ventilation_ix$INF_DURING_ICU_YN, c(FALSE,TRUE), labels = c("no infection", "infection"))
ventilation_ix$ventilation_severity <-factor(ventilation_ix$ventilation_severity, c(1,2,3), labels = c("mild", "moderate", "severe"))

cross<- table(ventilation_ix$INF_DURING_ICU_YN, ventilation_ix$ventilation_severity)
addmargins(cross)
round(100*prop.table(cross,2),digits = 0)
chisq.test(table(ventilation_ix$INF_DURING_ICU_YN, ventilation_ix$ventilation_severity))

plot(ventilation_ix$ventilation_severity, ventilation_ix$INF_DURING_ICU_YN)

## the code below is for analysis only - sparsely commented
##################################################################
#corticosteroid analysis - no infection on admission
#exclude anyone receiving anti malarials or antivirals and then mark the ones with steroids

# pick the dataset needed for the specific analysis depending on if patients with infection are excluded or not
filter_vent<-readRDS("filter_vent_no_inf.rds")
filter_vent<- readRDS("filter_vent_allow_infection.rds")

filtered_cortico_only<- filter_vent %>% filter(ICU_ANTIMALARIAL_YN == "FALSE") %>% filter(ICU_ANTIVIRALS_YN =="FALSE")
filtered_cortico_only<- filtered_cortico_only %>% mutate(cort_steroids = if_else(ICU_CORTICO_YN ==TRUE, 1,0), 
                                                         inf_during_icu = if_else(INF_DURING_ICU_YN ==TRUE, 1,0))

filtered_cortico_only<- filtered_cortico_only %>% dplyr::select(NEW_SUBJECT_ID, centre,inf_during_icu, comorb, ventilation_severity, 
                                                                cort_steroids, NEW_BMI,INF_MDR_PATHO_YN, ICU_CORTICO_DURATION_INT, 
                                                                INC_SEX_RAD, INC_AGE_INT, RESP_DURATION_INV_VENT_INT,ICU_RRT_DIAL_YN, ICU_INOTROPES_YN,
                                                                ICU_CORTICO_INTERV_INT, ICU_CORTICO_INDICATION_RAD, OUT_DEAD_DURING_ICU_YN, OUT_ICU_DURATION_INT, 
                                                                OUT_HOSP_DURATION_INT)



# dealing with NA values
# cortico steroids given, infection during icu and steroid indication 
# and MDR assumed unknown and False if not entered and infection during ICU
filtered_cortico_only<- filtered_cortico_only %>% 
  replace_na(list(INF_MDR_PATHO_YN = FALSE, ICU_CORTICO_INDICATION_RAD = "unknown", 
                  cort_steroids = 0, inf_during_icu=0))
#if steroids is 0 then steroid indication and interval and length should be zero

filtered_cortico_only<- filtered_cortico_only %>% mutate(ICU_CORTICO_DURATION_INT = if_else(cort_steroids ==0, as.numeric(0), as.numeric(ICU_CORTICO_DURATION_INT)),
                                                         ICU_CORTICO_INTERV_INT = if_else(cort_steroids ==0, as.numeric(0), as.numeric(ICU_CORTICO_INTERV_INT)))

# remove the rows where cort duration and interval are still NA
filtered_cortico_only<- filtered_cortico_only %>% drop_na(ICU_CORTICO_DURATION_INT)
filtered_cortico_only<- filtered_cortico_only %>% drop_na(ICU_CORTICO_INTERV_INT)
filtered_cortico_only<- filtered_cortico_only %>% rename(comorbidity = comorb,
                                                         age = INC_AGE_INT,
                                                         sex = INC_SEX_RAD,
                                                         ventilation_duration = RESP_DURATION_INV_VENT_INT)
filtered_cortico_only<- filtered_cortico_only %>% rename(inotropes = ICU_INOTROPES_YN,
                                                         rrt = ICU_RRT_DIAL_YN)
filtered_cortico_only<- filtered_cortico_only %>% drop_na(age)
save_filtered_cortico_only<- filtered_cortico_only

# impute the duration of ventilation

imputation_vent_duration<- filtered_cortico_only %>% dplyr::select(ventilation_duration, age, comorbidity, sex, cort_steroids, ventilation_severity)

#install.packages("mice")                        # Install mice package
library("mice") 
##### Impute data via predictive mean matching (single imputation)#####


imp_multi <- mice(imputation_vent_duration, m = 5, method = c("pmm", "", "", "", "", ""))  # Impute missing values multiple times
data_imp_multi_all <- complete(imp_multi,       # Store multiply imputed data
                               "repeated",
                               include = TRUE)

data_imp_multi <- data.frame(                   # Combine imputed Y and X1-X4 (for convenience)
  data_imp_multi_all[ , 1:6], imputation_vent_duration[, 2:6])
data_imp_multi<- data_imp_multi %>% mutate(ventilation_duration = (+ventilation_duration.1+ventilation_duration.2+ventilation_duration.3+ventilation_duration.4+ventilation_duration.5)/5)

imputed_values<- data_imp_multi %>% dplyr::select(age, sex, comorbidity, ventilation_duration, ventilation_severity, cort_steroids)

# add the other things back on

filtered_cortico_only<- filtered_cortico_only %>% dplyr::select(-age,-sex,-comorbidity,-ventilation_severity, -ventilation_duration, -cort_steroids)
filtered_cortico_only<- filtered_cortico_only %>% cbind(imputed_values)

saveRDS(filtered_cortico_only, "filtered_cortico_line260.rds")
filtered_cortico_only<-readRDS("filtered_cortico_line260.rds")



#############################filtered cortico only analysis:
###### steroid analysis propensity matching
#match - using "full" matching to improve the match
library(optmatch)
m_steroid_no_infection_full.out<-matchit(cort_steroids ~ comorbidity +ventilation_severity+sex+age+ventilation_duration+rrt+inotropes , data = filtered_cortico_only, method = "full", ratio = 1)
summary(m_steroid_no_infection_full.out)
plot(m_steroid_no_infection_full.out, type = "jitter", interactive = FALSE)
m_full.sum<-summary(m_steroid_no_infection_full.out)
plot(m_full.sum, var.order = "matched", threshold = c(0.1,0.05), abs = FALSE)


###compare estimated ps with true values
library(cobalt)
v<- data.frame(old = c("age", "rrt", "inotropes", "ventilation_severity", "ventilation_duration", "comorbidity", "sex_Male"), new = c("Age", "Renal replacement therapy","Inotropes/vasopressors", "Ventilation severity", "Ventilation duration", "Comorbidity", "Sex (Male)"))


tiff("Love plot steroids.tiff", width = 6, height = 5, units = 'in', res = 300)
love.plot(bal.tab(m_steroid_no_infection_full.out, m.threshold = 0.1,stat="mean.diffs", grid = TRUE, stars="raw", abs=FALSE),var.names = v, var.order = "unadjusted", title = "Covariate Balance Love Plot",sample.names = c("Unmatched", "Matched"),abs=F, stars = "raw" )
dev.off()

plot(m_steroid_no_infection_full.out, type = "qq", interactive = FALSE, which.xs = c("age","comorbidity","ventilation_severity","rrt", "inotropes", "ventilation_duration"))
m.data2<- match.data(m_steroid_no_infection_full.out)
library("lmtest")
library("sandwich")
library("survey")
fit2<-glm(inf_during_icu ~cort_steroids +ventilation_severity+comorbidity+sex+age+ventilation_duration, data = m.data2, weights = weights, family = "quasibinomial")
coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)

library("rms")
fit_rms<-lrm(inf_during_icu ~cort_steroids +ventilation_severity+comorbidity+sex+age+ventilation_duration+rrt+inotropes, data = m.data2, method = "lrm.fit",model = TRUE, na.action = na.delete,
             weights = weights, normwt = TRUE)

### density plots
p<-ggplot(match_data, aes(x=weight))+geom_density()+facet_grid(sex ~ .)


filtered_cortico_only<-filtered_cortico_only %>% rename(ventilatioventilation_duration)
m_plot.out<-matchit(cort_steroids ~ comorbidity +ventilation_severity+sex+age+ventilation_duration+rrt+inotropes , data = filtered_cortico_only, method = "full", ratio = 1)

library(ggpubr)
p1<-cobalt::bal.plot(m_plot.out, var.name = "age", which = "both")+ labs(title = "Distribution for age")+xlab("Age")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p2<-cobalt::bal.plot(m_plot.out, var.name = "sex", which = "both")+ labs(title = "Distribution for sex (Male)")+xlab("Sex (Male)")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p3<-cobalt::bal.plot(m_plot.out, var.name = "comorbidity", which = "both")+ labs(title = "Distribution for comorbidity score")+xlab("Comorbidity score")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p4<-cobalt::bal.plot(m_plot.out, var.name = "ventilation_severity", which = "both")+ labs(title = "Distribution for ventilation severity score")+xlab("Ventilation severity score")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p5<-cobalt::bal.plot(m_plot.out, var.name = "ventilation_duration", which = "both")+ labs(title = "Distribution for ventilation duration")+xlab("Ventilation duration (days)")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p6<-cobalt::bal.plot(m_plot.out, var.name = "rrt", which = "both")+ labs(title = "Distribution for renal replacement therapy")+xlab("Renal replacement therapy")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
p7<-cobalt::bal.plot(m_plot.out, var.name = "inotropes", which = "both")+ labs(title = "Distribution for inotropes/vasopressors")+xlab("Inotropes/vasopressors")+scale_fill_discrete(name = "Steroid treatment")+theme(text=element_text(size=5))
#p8<-gghistogram(match_data, x = "inf_during_icu", bins = 2, facet.by = "cort_steroids", fill= "cort_steroids" )+ labs(title = "Distribution infection")+xlab("Histograms of steroid usage")+theme(text=element_text(size=5))
#p8
figure<- ggarrange(p1,p2,p4,p7,p5,p6,p5,
                   ncol = 2, nrow= 4,
                   common.legend = TRUE, legend = "bottom")+theme(text=element_text(size=6))


figure
tiff("Density graph steroids.tiff", width = 5, height = 5, units = 'in', res = 300)
# Make plot
figure
dev.off()


### analysis of steroid matching
library(radiant.data)

#age
match_data %>% group_by(cort_steroids) %>% summarise(weighted_inf = weighted.mean(age, weights),weighted.sd(age, weights), n= n())
match_data %>% group_by(cort_steroids) %>% summarise(weighted_inf = weighted.mean(age, weights),weighted.sd(age, weights), n= n())
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(weighted_inf = mean(age),sd(age), n= n())

#sex
match_data_sex10<- match_data %>% mutate(sex_10 = if_else(sex =="Male",1,0))
match_data_sex10 %>% group_by(cort_steroids) %>%summarise(weighted.mean(sex_10, weights), weighted.sd(sex_10, weights), n= n())
filtered_cortico_only %>% group_by(cort_steroids, sex) %>% summarise(n= n())

#comorbidity score
library(spatstat)
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.median(comorbidity, weights), quarts = weighted.quantile(comorbidity, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=median(comorbidity), quarts = quantile(comorbidity, probs = seq(0,1,0.25)))


#ventilation severity
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.median(ventilation_severity, weights), quarts = weighted.quantile(ventilation_severity, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=median(ventilation_severity), quarts = quantile(ventilation_severity, probs = seq(0,1,0.25)))

#ventilation duration
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.median(ventilation_duration, weights), quarts = weighted.quantile(ventilation_duration, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=median(ventilation_duration), quarts = quantile(ventilation_duration, probs = seq(0,1,0.25)))
wilcox.test(ventilation_duration~cort_steroids, data=match_data)
wilcox.test(ventilation_duration~cort_steroids, data=filtered_cortico_only)


#rrt
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.mean(rrt, weights), quarts = weighted.quantile(rrt, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=mean(rrt), quarts = quantile(rrt, probs = seq(0,1,0.25)))

#inotropes
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.mean(inotropes, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=mean(inotropes))


#LOS ICU
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.median(OUT_ICU_DURATION_INT, weights), quarts = weighted.quantile(OUT_ICU_DURATION_INT, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=median(OUT_ICU_DURATION_INT, na.rm = TRUE), quarts = quantile(OUT_ICU_DURATION_INT,na.rm= TRUE, probs = seq(0,1,0.25)))

#LOS Hospital
match_data %>% group_by(cort_steroids) %>% summarise(com_med=weighted.median(OUT_HOSP_DURATION_INT, weights), quarts = weighted.quantile(OUT_HOSP_DURATION_INT, weights))
filtered_cortico_only %>% group_by(cort_steroids) %>% summarise(com_med=median(OUT_HOSP_DURATION_INT, na.rm = TRUE), quarts = quantile(OUT_HOSP_DURATION_INT,na.rm= TRUE, probs = seq(0,1,0.25)))



dta_m_ster<- match.data(m_steroid_no_infection_full.out)
dim(dta_m_ster)


inf_severity<- filter_vent %>% dplyr::select(centre, INF_SEVERITY,NEW_SUBJECT_ID)
#inf_severity<- inf_severity %>% rename(age = INC_AGE_INT, ventilation_duration = RESP_DURATION_INV_VENT_INT)
#write.csv(filter_vent, "filter_vent.csv", row.names = FALSE)


dta_m_ster_amended<- dta_m_ster %>% left_join(inf_severity, by =c("centre", "NEW_SUBJECT_ID"))


##analysis of indication for steroids:
saveRDS(dta_m_ster, "dta_m_ster.rds")
saveRDS(dta_m_ster_amended, "dta_m_ster_amended.rds")

dta_m_ster_amended<- save_dta_amended
dta_m_ster_amended<- match_data
# indications for the steroids
indications_ster<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%  group_by(ICU_CORTICO_INDICATION_RAD) %>%  tally()
indications_ster$n<- indications_ster$n/507
indications_ster_no_infinicu <- dta_m_ster_amended %>% filter(cort_steroids==1) %>% filter(inf_during_icu ==0) %>%  group_by(ICU_CORTICO_INDICATION_RAD) %>%  tally()
indications_ster_no_infinicu$n<- indications_ster_no_infinicu$n/147
indications_ster_yes_infinicu <- dta_m_ster_amended %>% filter(cort_steroids==1) %>% filter(inf_during_icu ==1) %>%  group_by(ICU_CORTICO_INDICATION_RAD) %>%  tally()
indications_ster_yes_infinicu$n<- indications_ster_yes_infinicu$n/360




interval_ster_yes_inf<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%filter(inf_during_icu==1) %>%  summarise(mean_int = mean(ICU_CORTICO_INTERV_INT), 
                                                                                                            med_int = median(ICU_CORTICO_INTERV_INT),
                                                                                                            quart = quantile(ICU_CORTICO_INTERV_INT, probs = seq(0,1,0.25)))
interval_ster_no_inf<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%filter(inf_during_icu==0) %>%  summarise(mean_int = mean(ICU_CORTICO_INTERV_INT), 
                                                                                                           med_int = median(ICU_CORTICO_INTERV_INT),
                                                                                                           quant = quantile(ICU_CORTICO_INTERV_INT, probs = seq(0,1,0.25)))
interval_ster_all<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%  summarise(mean_int = mean(ICU_CORTICO_INTERV_INT), 
                                                                                                                med_int = median(ICU_CORTICO_INTERV_INT),
                                                                                                                quant = quantile(ICU_CORTICO_INTERV_INT, probs = seq(0,1,0.25)))
duration_ster_yes_inf<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%filter(inf_during_icu==1) %>%  summarise(mean_int = mean(ICU_CORTICO_DURATION_INT), 
                                                                                                                    med_int = median(ICU_CORTICO_DURATION_INT),
                                                                                                                    quart = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
duration_ster_no_inf<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%filter(inf_during_icu==0) %>%  summarise(mean_int = mean(ICU_CORTICO_DURATION_INT), 
                                                                                                                   med_int = median(ICU_CORTICO_DURATION_INT),
                                                                                                                   quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
duration_ster_all<- dta_m_ster_amended %>% filter(cort_steroids==1) %>%  summarise(mean_int = mean(ICU_CORTICO_DURATION_INT), 
                                                                                   med_int = median(ICU_CORTICO_DURATION_INT),
                                                                                   quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))




#outcome death analysis
outcome_dead<- dta_m_ster_amended %>% group_by(cort_steroids, OUT_DEAD_DURING_ICU_YN) %>% tally()
outcome_sepsis<- dta_m_ster_amended %>% group_by(cort_steroids, INF_SEVERITY) %>% tally()
chisq.test(table(dta_m_ster$OUT_DEAD_DURING_ICU_YN, dta_m_ster$OUT_DEAD_DURING_ICU_YN))
save_dta_amended<- dta_m_ster_amended


###################################################
# Do patients on steroids get more infections?
# change the laebls in the matched df to allow chi-squared analysis - 
dta_m_ster_amended$inf_during_icu<-factor(dta_m_ster_amended$inf_during_icu, c(0,1), labels = c("no infection", "infection"))
dta_m_ster_amended$cort_steroids <-factor(dta_m_ster_amended$cort_steroids, c(0,1), labels = c("no steroids", "steroids"))
cross<- table(dta_m_ster$inf_during_icu, dta_m_ster$cort_steroids)
addmargins(cross)
round(100*prop.table(cross,2),digits = 0)
chisq.test(table(dta_m_ster$inf_during_icu, dta_m_ster$cort_steroids))
#####################################################



# duration of steroids analysis
dta_m_ster_amended<- save_dta_amended
wilcox.test(ICU_CORTICO_DURATION_INT~inf_during_icu, data=dta_m_ster_amended)
m<-dta_m_ster_amended %>% group_by(cort_steroids,inf_during_icu) %>% summarise(duration_med = median(ICU_CORTICO_DURATION_INT), sd_duration = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)),
                                                                                number=length(ICU_CORTICO_DURATION_INT))
j<-dta_m_ster_amended %>% group_by(cort_steroids,inf_during_icu) %>% summarise(duration_med = median(ICU_CORTICO_INTERV_INT), sd_duration = quantile(ICU_CORTICO_INTERV_INT, probs = seq(0,1,0.25)),
                                                                               number=length(ICU_CORTICO_INTERV_INT))


n<-dta_m_ster_amended %>% filter(cort_steroids==1 & ICU_CORTICO_DURATION_INT>0) %>% group_by(inf_during_icu) %>% summarise(median = median(ICU_CORTICO_DURATION_INT),
                                                                                                                            sd_duration = sd(ICU_CORTICO_DURATION_INT), 
                                                                                                                            number=length(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
x<- dta_m_ster_amended %>% group_by(cort_steroids, inf_during_icu, OUT_DEAD_DURING_ICU_YN) %>% tally()
y<-dta_m_ster_amended %>% group_by(cort_steroids) %>% summarise(median_d = median(ventilation_duration), iqrd = quantile(ventilation_duration, probs = seq(0,1,0.25)))


####AMR analysis
steroid_matched<-dta_m_ster_amended %>% rename(amr_inf= INF_MDR_PATHO_YN)
steroid_matched_inf<- steroid_matched %>% filter(inf_during_icu ==1)
steroid_matched$amr_inf<-factor(steroid_matched$amr_inf, c(FALSE,TRUE), labels = c("no amr", "amr"))
steroid_matched$cort_steroids<-factor(steroid_matched$cort_steroids, c(0,1), labels = c("no steroids", "steroids"))


cross<- table(steroid_matched$amr_inf, steroid_matched$cort_steroids)
addmargins(cross)
round(100*prop.table(cross,2),digits = 0)
chisq.test(table(steroid_matched$amr_inf, steroid_matched$cort_steroids))


#### further analysis corticosteroid cohort

df_univariate<- filter_vent_allow_infection # use reduced_c_df for all patients, 
#filter_vent for all that are ventilated and dont have an infection at admission and 
#filter_vent_allow_infectionf or all that are ventilated

df_univariate<- df_univariate %>% mutate(inf_during_icu = if_else(INF_DURING_ICU_YN ==TRUE, 1,0)) 
df_univariate<- df_univariate %>% mutate(cort_steroids = if_else(ICU_CORTICO_YN ==TRUE, 1,0)) 
df_univariate$inf_during_icu<- ifelse(is.na(df_univariate$inf_during_icu), 0, df_univariate$inf_during_icu)
df_univariate %>% dplyr::group_by(inf_during_icu) %>% tally()

# 6 NA in age - will be removed
age_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_age = mean(INC_AGE_INT, na.rm = TRUE), 
                                                                            median_age = median(INC_AGE_INT, na.rm = TRUE),
                                                                            iqr_age = stats::IQR(INC_AGE_INT, na.rm = TRUE),
                                                                            sd_age = sd(INC_AGE_INT, na.rm = TRUE), n = max(row_number()))

# no NA in comorbidity
comorb_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_com = mean(comorb, na.rm = TRUE), 
                                                                               median_com = median(comorb, na.rm = TRUE),
                                                                               iqr_com = stats::IQR(comorb, na.rm = TRUE),
                                                                               sd_com = sd(comorb, na.rm = TRUE), n = max(row_number()))
#no NA in vent severity
ventsev_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_ventsev = mean(ventilation_severity, na.rm = TRUE), 
                                                                                median_ventsev = median(ventilation_severity, na.rm = TRUE),
                                                                                iqr_ventsev = stats::IQR(ventilation_severity, na.rm = TRUE),
                                                                                sd_ventsev = sd(ventilation_severity, na.rm = TRUE), n = max(row_number()))

#several NA here - use imputation to do this
imp_df<- df_univariate %>% dplyr::select(RESP_DURATION_INV_VENT_INT, INC_AGE_INT, comorb, INC_SEX_RAD,cort_steroids, ventilation_severity)
imp_df<- imp_df %>% mutate(cort_steroids = if_else(is.na(cort_steroids), 0, cort_steroids))
imp_df<- imp_df %>% drop_na(INC_AGE_INT)

imp_multi_b <- mice(imp_df, m = 5, method = c("pmm", "", "", "", "", ""))  # Impute missing values multiple times
data_imp_multi_all_b <- complete(imp_multi_b,       # Store multiply imputed data
                               "repeated",
                               include = TRUE)

data_imp_multi_b <- data.frame(                   # Combine imputed Y and X1-X4 (for convenience)
  data_imp_multi_all_b[ , 1:6], imp_df[, 2:6])
data_imp_multi_b<- data_imp_multi_b %>% mutate(ventilation_duration = (RESP_DURATION_INV_VENT_INT.1+RESP_DURATION_INV_VENT_INT.2+RESP_DURATION_INV_VENT_INT.3+RESP_DURATION_INV_VENT_INT.4+RESP_DURATION_INV_VENT_INT.5)/5)

imputed_values<- data_imp_multi_b %>% dplyr::select(age, sex, comorbidity, ventilation_duration, ventilation_severity, cort_steroids)

# add the other things back on

filtered_cortico_only<- filtered_cortico_only %>% dplyr::select(-age,-sex,-comorbidity,-ventilation_severity, -ventilation_duration, -cort_steroids)
filtered_cortico_only<- filtered_cortico_only %>% cbind(imputed_values)




ventdur_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_ventdur = mean(RESP_DURATION_INV_VENT_INT, na.rm = TRUE), 
                                                                                median_ventdur = median(RESP_DURATION_INV_VENT_INT, na.rm = TRUE),
                                                                                iqr_ventdur = stats::IQR(RESP_DURATION_INV_VENT_INT, na.rm = TRUE),
                                                                                sd_ventdur = sd(RESP_DURATION_INV_VENT_INT, na.rm = TRUE), n = max(row_number()))

df_univariate<- df_univariate %>% mutate(sex_m10 = if_else(INC_SEX_RAD =="Male", 1,0)) 
sexm_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_male = mean(sex_m10, na.rm = TRUE))

df_univariate<- df_univariate %>% mutate(icu_cortico = if_else(ICU_CORTICO_YN ==TRUE, 1,0)) 
cortico_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_cortico = mean(icu_cortico, na.rm = TRUE))

df_univariate<- df_univariate %>% mutate(vent_yn= if_else(RESP_INV_VENT_YN ==TRUE, 1,0)) 
vent_assoc<- df_univariate %>% dplyr::group_by(inf_during_icu) %>% summarise(mean_vent = mean(vent_yn, na.rm = TRUE))

wilcox.test(RESP_DURATION_INV_VENT_INT~inf_during_icu, data=df_univariate)
wilcox.test(comorb~inf_during_icu, data=df_univariate)
wilcox.test(ventilation_severity~inf_during_icu, data=df_univariate)
wilcox.test(sex_m10~inf_during_icu, data=df_univariate)
wilcox.test(INC_AGE_INT~inf_during_icu, data=df_univariate)
wilcox.test(vent_yn~inf_during_icu, data = df_univariate)



wilcox.test(icu_cortico~inf_during_icu, data=df_univariate)
dta_m_ster_amended<- save_dta_amended
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0&OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0& OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

dta_m_ster_amended %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0 & OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0 & OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

dta_m_ster_amended %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
dta_m_ster_amended %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
dta_m_ster_amended %>% filter(ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter( OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  filter( OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
dta_m_ster_amended %>%  
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
dta_m_ster_amended %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

# now for the whole cohort:
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0&OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0& OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

filter_vent %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0 & OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0 & OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
filter_vent %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))




filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0&OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0& OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==1 &ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

filtered_cortico_only %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0 & OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0 & OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter(cort_steroids==0) %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))

filtered_cortico_only %>% 
  summarise( median = median(comorbidity), quant = quantile(comorbidity, probs = seq(0,1,0.25)))
filtered_cortico_only %>% 
  summarise( median = median(age), quant = quantile(age, probs = seq(0,1,0.25)))
filtered_cortico_only %>% filter(ICU_CORTICO_DURATION_INT>0) %>% 
  summarise( median = median(ICU_CORTICO_DURATION_INT), quant = quantile(ICU_CORTICO_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter( OUT_HOSP_DURATION_INT>0) %>% 
  summarise( median = median(OUT_HOSP_DURATION_INT), quant = quantile(OUT_HOSP_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  filter( OUT_ICU_DURATION_INT>0) %>% 
  summarise( median = median(OUT_ICU_DURATION_INT), quant = quantile(OUT_ICU_DURATION_INT, probs = seq(0,1,0.25)))
filtered_cortico_only %>%  
  summarise( median = median(ventilation_duration), quant = quantile(ventilation_duration, probs = seq(0,1,0.25)))
filtered_cortico_only %>% 
  summarise( median = median(ventilation_severity), quant = quantile(ventilation_severity, probs = seq(0,1,0.25)))




red_data_centre<- red_data_centre %>%replace_na(list(cort_steroids=0, ICU_CORTICO_DURATION_INT = 0))
red_data_centre<- red_data_centre %>% mutate(centre = paste(NEW_COUNTRY_ID, "_",NEW_SITE_ID))

red_data_centre<- red_data_centre %>% rename(cort_steroids = ICU_CORTICO_YN)
percent_per_centre<-red_data_centre %>% group_by(centre) %>%  summarise(n_steroids = sum(cort_steroids), 
                                                                        n = length(row_number()))
percent_per_centre<- percent_per_centre %>% mutate(percentage_steroids = n_steroids/n*100)
additional_per_centre <- red_data_centre %>% group_by(centre) %>% filter(cort_steroids==1) %>% 
  summarise(n = length(row_number()), med_length = median(ICU_CORTICO_DURATION_INT),iqr_low = quantile(ICU_CORTICO_DURATION_INT, probs = c(0.25)), high = quantile(ICU_CORTICO_DURATION_INT, probs=c(0.75)))


percent_per_centre<-percent_per_centre %>% left_join(additional_per_centre, by="centre")
percent_per_centre<- percent_per_centre %>% rename(lower_limit = iqr_low, upper_limit = high) %>% mutate(iqr_length = upper_limit-lower_limit)
percent_per_centre<- percent_per_centre %>% replace_na(list(med_length=0, iqr_length=0, upper_limit=0, lower_limit=0))



############################################antiviral analysis
#############################################################
#tocilizumab only
# take away all those that received steroids
filter_vent<- read.csv("filter_vent.csv", stringsAsFactors = FALSE)
filtered_antiviral_only<- filter_vent %>% filter(ICU_CORTICO_YN == "FALSE")
#filtered antiviral is either no treatment or treatment with toci in some way
filtered_antiviral_only<- filtered_antiviral_only %>% filter(((ICU_ANTIMALARIAL_YN ==TRUE | ICU_ANTIVIRALS_YN ==TRUE) &
                                                 (ICU_ANTIVIRALS_RAD == "Tocilizumab"| ICU_OTHER_ANTIVIRALS_RAD == "Tocilizumab"))|
                                                (ICU_ANTIMALARIAL_YN ==FALSE & ICU_ANTIVIRALS_YN == FALSE))

filtered_antiviral_only<- filtered_antiviral_only %>% 
  mutate(antivirals = if_else(((ICU_ANTIMALARIAL_YN ==TRUE | ICU_ANTIVIRALS_YN ==TRUE) & ((ICU_ANTIVIRALS_RAD == "Tocilizumab")|(ICU_OTHER_ANTIVIRALS_RAD =="Tocilizumab"))), 1,0), 
         inf_during_icu = if_else(INF_DURING_ICU_YN ==TRUE, 1,0), icu_death = ifelse(OUT_DEATH_DURING_ICU_YN ==TRUE,1,0))


filtered_antiviral_only<- filtered_antiviral_only %>% replace_na(list(icu_death = FALSE))
filtered_antiviral_only<- filtered_antiviral_only %>% dplyr::select(inf_during_icu, comorb, ventilation_severity, antivirals, OUT_ICU_DURATION_INT, OUT_HOSP_DURATION_INT, icu_death,  INC_AGE_INT, INC_SEX_RAD, RESP_DURATION_INV_VENT_INT, ICU_RRT_DIAL_YN, ICU_INOTROPES_YN)

filtered_antiviral_only<- filtered_antiviral_only %>% replace_na(list(antivirals = 0, ICU_RRT_DIAL_YN = 0))
filtered_antiviral_only<- filtered_antiviral_only %>% drop_na(INC_AGE_INT)
filtered_antiviral_only<- filtered_antiviral_only %>% drop_na(inf_during_icu)

filtered_antiviral_only<- filtered_antiviral_only %>% rename(ventilation_duration = RESP_DURATION_INV_VENT_INT)

# impute ventilation duration
impute_vent_duration_viral<- filtered_antiviral_only %>% dplyr::select(ventilation_duration, INC_AGE_INT, comorb, INC_SEX_RAD, antivirals, ventilation_severity)

#install.packages("mice")                        # Install mice package
library("mice") 
##### Impute data via predictive mean matching (single imputation)#####

imp_single <- mice(imputation_vent_duration, m = 1, method = "pmm") # Impute missing values
data_imp_single <- complete(imp_single)         # Store imputed data
# head(data_imp_single)                         # First 6 rows of our imputed data


imp_multi <- mice(impute_vent_duration_viral, m = 5, method = c("pmm", "", "", "", "", ""))  # Impute missing values multiple times
data_imp_multi_all <- complete(imp_multi,       # Store multiply imputed data
                               "repeated",
                               include = TRUE)

data_imp_multi <- data.frame(                   # Combine imputed Y and X1-X4 (for convenience)
  data_imp_multi_all[ , 1:6], impute_vent_duration_viral[, 2:6])
data_imp_multi<- data_imp_multi %>% mutate(ventilation_duration = (+ventilation_duration.1+ventilation_duration.2+ventilation_duration.3+ventilation_duration.4+ventilation_duration.5)/5)

imputed_values<- data_imp_multi %>% dplyr::select(INC_AGE_INT, INC_SEX_RAD, comorb, ventilation_duration, ventilation_severity, antivirals)

# add the other things back on

filtered_antiviral_only<- filtered_antiviral_only %>% dplyr::select(-INC_AGE_INT,-INC_SEX_RAD,-comorb,-ventilation_severity, -ventilation_duration, -antivirals)
filtered_antiviral_only<- filtered_antiviral_only %>% cbind(imputed_values)


#filtered_antiviral_only<- filtered_antiviral_only %>% drop_na()

#####matching
filtered_antiviral_only<- filtered_antiviral_only %>% rename(comorbidity = comorb,age = INC_AGE_INT,sex = INC_SEX_RAD)



m_antiviral.out<-matchit(antivirals ~ comorbidity +ventilation_severity+age+sex+ventilation_duration+ICU_RRT_DIAL_YN+ICU_INOTROPES_YN , data = filtered_antiviral_only, method = "full", ratio = 1)
summary(m_antiviral.out)
m_antiviral.sum<- summary(m_antiviral.out)
plot(m_antiviral.sum, var.order = "data", abs = FALSE, threshold = c(0.05, NA) )


dta_m_vir<- match.data(m_antiviral.out)
dim(dta_m_vir)
los_hospital<- filter_vent %>% dplyr::select(OUT_HOSP_DURATION_INT, NEW_BMI, OUT_ICU_DURATION_INT)
dta_m_vir_amended<- dta_m_vir %>% left_join(los_hospital, by = c("NEW_BMI", "OUT_ICU_DURATION_INT"))
filtered_antiviral_only<- filtered_antiviral_only%>% left_join(los_hospital, by = c("NEW_BMI", "OUT_ICU_DURATION_INT"))

plot(m_antiviral.out,type ="density" , interactive = FALSE, which.xs = c( "comorbidity", "ventilation_duration", "ventilation_severity"))

tiff("Density plot antiviral.tiff", width = 5, height = 4, units = 'in', res = 300)
plot(m_antiviral.out,type ="density" , interactive = FALSE, which.xs = c( "comorbidity", "ventilation_duration", "ventilation_severity"))
# Make plot
dev.off()

saveRDS(dta_m_vir, "dta_m_vir.rds")
dta_m_vir<- readRDS("dta_m_vir.rds")
###table generation for paper
renamed_antivial_only<-filtered_antiviral_only %>% rename(rrt = ICU_RRT_DIAL_YN, inotropes = ICU_INOTROPES_YN)

vir_plot.out<-matchit(antivirals ~ comorbidity +ventilation_severity+sex+age+ventilation_duration+rrt+inotropes , data = renamed_antivial_only, method = "full", ratio = 1)

library(ggpubr)
p1<-cobalt::bal.plot(vir_plot.out, var.name = "age", which = "both")+ labs(title = "Distribution for age")+xlab("Age")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p2<-cobalt::bal.plot(vir_plot.out, var.name = "sex", which = "both")+ labs(title = "Distribution for sex (Male)")+xlab("Sex (Male)")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p3<-cobalt::bal.plot(vir_plot.out, var.name = "comorbidity", which = "both")+ labs(title = "Distribution for comorbidity score")+xlab("Comorbidity score")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p4<-cobalt::bal.plot(vir_plot.out, var.name = "ventilation_severity", which = "both")+ labs(title = "Distribution for ventilation severity score")+xlab("Ventilation severity score")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p5<-cobalt::bal.plot(vir_plot.out, var.name = "ventilation_duration", which = "both")+ labs(title = "Distribution for ventilation duration")+xlab("Ventilation duration (days)")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p6<-cobalt::bal.plot(vir_plot.out, var.name = "rrt", which = "both")+ labs(title = "Distribution for renal replacement therapy")+xlab("Renal replacement therapy")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
p7<-cobalt::bal.plot(vir_plot.out, var.name = "inotropes", which = "both")+ labs(title = "Distribution for inotropes/vasopressors")+xlab("Inotropes/vasopressors")+scale_fill_discrete(name = "Antiviral treatment")+theme(text=element_text(size=5))
#p8<-gghistogram(match_data, x = "inf_during_icu", bins = 2, facet.by = "cort_steroids", fill= "cort_steroids" )+ labs(title = "Distribution infection")+xlab("Histograms of steroid usage")+theme(text=element_text(size=5))
#p8
figure<- ggarrange(p1,p2,p4,p7,p5,p6,p5,
                   ncol = 2, nrow= 4,
                   common.legend = TRUE, legend = "bottom")+theme(text=element_text(size=6))


figure
tiff("Density graph antiviral.tiff", width = 5, height = 5, units = 'in', res = 300)
# Make plot
figure
dev.off()

v<- data.frame(old = c("age", "rrt", "inotropes", "ventilation_severity", "ventilation_duration", "comorbidity", "sex_Male"), new = c("Age", "Renal replacement therapy","Inotropes/vasopressors", "Ventilation severity", "Ventilation duration", "Comorbidity", "Sex (Male)"))


tiff("Love plot antivirals.tiff", width = 6, height = 5, units = 'in', res = 300)
love.plot(bal.tab(vir_plot.out, m.threshold = 0.1,stat="mean.diffs", grid = TRUE, stars="raw", abs=FALSE),var.names = v, var.order = "unadjusted", title = "Covariate Balance Love Plot",sample.names = c("Unmatched", "Matched"),abs=F, stars = "raw" )
dev.off()


dta_m_vir %>% group_by(antivirals,ICU_RRT_DIAL_YN) %>% summarise(max_rrt = max(row_number()))
dta_m_vir %>% group_by(antivirals,ICU_INOTROPES_YN) %>% summarise(max_rrt = max(row_number()))
filtered_antiviral_only %>% group_by(antivirals,ICU_RRT_DIAL_YN) %>% summarise(max_rrt = max(row_number()))
filtered_antiviral_only %>% group_by(antivirals, ICU_INOTROPES_YN) %>% summarise(max_rrt = max(row_number()))


toci_match<- match.data(m_antiviral.out)


#age
toci_match %>% group_by(antivirals) %>% summarise(weighted_inf = weighted.mean(age, weights),weighted.sd(age, weights), n= n())
toci_match %>% group_by(antivirals) %>% summarise(weighted_inf = weighted.mean(age, weights),weighted.sd(age, weights), n= n())
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(weighted_inf = mean(age),sd(age), n= n())

#sex
toci_match_sex10<- toci_match %>% mutate(sex_10 = if_else(sex =="Male",1,0))
toci_match_sex10 %>% group_by(antivirals) %>%summarise(weighted.mean(sex_10, weights), weighted.sd(sex_10, weights), n= n())
filtered_antiviral_only %>% group_by(antivirals, sex) %>% summarise(n= n())

#comorbidity score
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.median(comorbidity, weights), quarts = weighted.quantile(comorbidity, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=median(comorbidity), quarts = quantile(comorbidity, probs = seq(0,1,0.25)))


#ventilation severity
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.median(ventilation_severity, weights), quarts = weighted.quantile(ventilation_severity, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=median(ventilation_severity), quarts = quantile(ventilation_severity, probs = seq(0,1,0.25)))

#ventilation duration
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.median(ventilation_duration, weights), quarts = weighted.quantile(ventilation_duration, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=median(ventilation_duration), quarts = quantile(ventilation_duration, probs = seq(0,1,0.25)))
wilcox.test(ventilation_duration~antivirals, data=toci_match)
wilcox.test(ventilation_duration~antivirals, data=filtered_antiviral_only)


#rrt
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.mean(ICU_RRT_DIAL_YN, weights), quarts = weighted.quantile(ICU_RRT_DIAL_YN, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=mean(ICU_RRT_DIAL_YN), quarts = quantile(ICU_RRT_DIAL_YN, probs = seq(0,1,0.25)))

#inotropes
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.mean(ICU_INOTROPES_YN, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=mean(ICU_INOTROPES_YN))


#LOS ICU
toci_match %>% group_by(antivirals) %>% summarise(com_med=weighted.median(OUT_ICU_DURATION_INT, weights), quarts = weighted.quantile(OUT_ICU_DURATION_INT, weights))
filtered_antiviral_only %>% group_by(antivirals) %>% summarise(com_med=median(OUT_ICU_DURATION_INT, na.rm = TRUE), quarts = quantile(OUT_ICU_DURATION_INT,na.rm= TRUE, probs = seq(0,1,0.25)))

