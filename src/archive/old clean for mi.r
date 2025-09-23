
UNITE_2020_corrected_preimputation <- UNITE_2020_corrected %>%
  select(-c(
    # Steroid-related variables
    "ICU_CORTICO_INDICATION_SHOCK",
    "ICU_CORTICO_INDICATION_HYPERINFLAMMATION",
    "ICU_CORTICO_INDICATION_PNEUMONITIS",
    "ICU_CORTICO_INDICATION_PRE_EXISTING_CONDITION",
    "ICU_CORTICO_INDICATION_OTHER",
    "ICU_CORTICO_ICU_INITIATION_AT_ADMISSION_YN",
    "ICU_CORTICO_ICU_INITIATION_PRIOR_ADMISSION_ICU_YN",
    "ICU_CORTICO_ICU_INITIATION_DURING_ADMISSION_ICU_YN",
    
    # Comorbidity scores
    "cardiac_d",
    "liver_d",
    "neuro_d",
    "diabetes_d",
    "kidney_d",
    "htn_d",
    "asthma_d",
    "resp_d",
    "immunosup_d",
    "hiv_d",
    
    # Subject identifiers
    "NEW_SUBJECT_ID",
    "NEW_CENTRE_ID",
    "NEW_PATIENT_ID",
    "NEW_COUNTRY_ID",

    #  ICU bloods
    "ICU_SUP_WHITE_CELL_INT",
    "ICU_SUP_NEUTRO_INT",
    "ICU_SUP_LYMPH",
    "ICU_SUP_FERRITINE_DEC",
    "COAG_DIMERS_INT",
    "COAG_PLATELETS_INT",

    # resp metrics already used for resp failure score 
    "RESP_PRONE_YN",
    "RESP_ECMO_YN",
    "RESP_NEUROM_BLOC_YN",
    "RESP_INTUBATED_YN",
    "RESP_INTUB_DAYS_AFT_ADM_INT",
    "RESP_DURATION_INV_VENT_INT",
    "RESP_INV_VENT_YN",
    "RESP_NI_VENT_YN",

    # duplicate outcomes for follow up at 60 days
    "OUTCOME_LD_OTHER",
    "OUTCOME_LD_TRANSFERRED",
    "OUTCOME_LD_DISCHARGED",
    "OUT_HOSP_DURATION_OVERALL_INT",

    # Drop varibles which are characters and not suitable for factors 
    ICU_ADM_DIAG_RAD, 
    ICU_CORTICO_INDICATION_RAD,
    ICU_ANTIVIRALS_RAD,
    OUTCOME_LD, 
    ICU_OTHER_ANTIVIRALS_RAD,NEW_CENTRE_ID,
    INC_HIV_YN,
    ICU_SUPP_TYPE_RAD,
    INC_DIABETES_YN,

    #Outcomes not suitable for imputation
    OUTCOME_LD_DEATH,
    OUT_HOSP_DURATION_INT,
    OUT_ICU_DURATION_INT,

    # colinear CRP 
    ICU_CRP_RAD50,
    ICU_CRP_RAD

  ))
