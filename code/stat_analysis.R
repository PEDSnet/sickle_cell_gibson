# survival analysis
# source("code/EDA.R")
# install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
# install.packages("ggsurvfit")
# the wasy survival is coded here is different from table 1
# this is overall survival, death of any causes 

cohort_covars <- results_tbl("analytics_dataset") %>% collect() %>% 
                    mutate(age_at_ce = as.numeric(difftime(transplant_date, birth_date, units = "days"))/365.25, 
                            ethnicity = factor(ethnicity, levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Other", "Unknown"))) %>%
          mutate(survival = if_else(survival, 1, 0), 
                survival_1yr = if_else(death_days_since_ce >= 365.25 | is.na(death_days_since_ce), 0, 1),
                survival_3yr = if_else(death_days_since_ce >= (365.25*3) | is.na(death_days_since_ce), 0, 1),
                survival_5yr = if_else(death_days_since_ce >= (365.25*5) | is.na(death_days_since_ce), 0, 1),
                disease_relapse_days_since_ce = as.numeric(difftime(graft_fail_date, transplant_date, units = "days"))) %>%
          mutate(age_at_ce = round(as.numeric(difftime(transplant_date, birth_date, unit = "days"))/365.25), 
         donor_relation = if_else(is.na(donor_relation), "Unknown", donor_relation),
         match_status = if_else(is.na(match_status), "Unknown", match_status),
         second_transplant = if_else(is.na(second_transplant_date), "No", "Yes"),
         gender = if_else(gender == "FEMALE", "Female", "Male"),
         scd_type = recode(scd_type, "aa" = "aplastic anemia", "dba" = "Diamond-Blackfan anemia", "bta" = "beta thalassemia major", "scd" = "sickle cell disease")) %>%
        mutate(CD3_reconstitution = case_when(CD3_reconstitute_3mon ~ "0-3 months",
                                        CD3_reconstitute_6mon ~ "3-6 months",
                                        CD3_reconstitute_9mon ~ "6-9 months",
                                        CD3_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         CD4_reconstitution = case_when(CD4_reconstitute_3mon ~ "0-3 months",
                                        CD4_reconstitute_6mon ~ "3-6 months",
                                        CD4_reconstitute_9mon ~ "6-9 months",
                                        CD4_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         CD8_reconstitution = case_when(CD8_reconstitute_3mon ~ "0-3 months",
                                        CD8_reconstitute_6mon ~ "3-6 months",
                                        CD8_reconstitute_9mon ~ "6-9 months",
                                        CD8_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         IgM_reconstitution = case_when(IgM_reconstitute_3mon~ "0-3 months",
                                        IgM_reconstitute_6mon ~ "3-6 months",
                                        IgM_reconstitute_3mon ~ "6-9 months",
                                        IgM_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months")) %>%
         mutate(across(ends_with("reconstitution"), ~factor(.x, levels = c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "Never/after 12 months")))) %>%
         mutate(vod_subgroup = if_else(vod, scd_type, "No VOD")) %>% 
         mutate(vod_subgroup = factor(vod_subgroup, levels = c("beta thalassemia major", "sickle cell disease", "No VOD")))%>%
         mutate(busulfan_and_vod = case_when(has_busulfan & vod ~ "busulfan & VOD",
                                            has_busulfan & !vod ~ "busulfan & No VOD", 
                                            TRUE ~ "Other categories")) %>%
         mutate(across(c("has_busulfan"), ~if_else(is.na(.x) | !(.x), "No", "Yes")))

cohort_covars <- cohort_covars %>% 
  mutate(#survival = if_else(survival, "Alive", "Death"), 
         #censoring status 0=censored, 1=dead
         survival_time_1yr_censored = case_when(death_days_since_ce <= 365.25 ~ death_days_since_ce, 
                                                death_days_since_ce > 365.25 ~ 365.25,
                                                is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25, na.rm = TRUE),
                                                TRUE ~ NA_real_),
         survival_time_3yr_censored = case_when(death_days_since_ce <= 365.25*3 ~ death_days_since_ce, 
                                                death_days_since_ce > 365.25*3 ~ 365.25*3,
                                                is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*3, na.rm = TRUE),
                                                TRUE ~ NA_real_),
        survival_time_5yr_censored = case_when(death_days_since_ce <= 365.25*5 ~ death_days_since_ce, 
                                               death_days_since_ce > 365.25*5 ~ 365.25*5, 
                                               is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*5, na.rm = TRUE),
                                               TRUE ~ NA_real_),
         disease_relapse_days_since_ce_5yr_censored = case_when(disease_relapse_days_since_ce <= 365.25*5 ~ disease_relapse_days_since_ce, 
                                                                disease_relapse_days_since_ce > 365.25*5 ~ 365.25*5,
                                                          is.na(disease_relapse_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*5, na.rm = TRUE),
                                                TRUE ~ NA_real_),
         disease_relapse_5yr = if_else(disease_relapse_days_since_ce > 365.25*5 | is.na(disease_relapse_days_since_ce), FALSE, TRUE))

cohort_covars <- cohort_covars %>% mutate(ferritin_pre = factor(ferritin_pre, levels = c("low", "moderate", "high")),
         ferritin_pre_max = factor(ferritin_pre_max, levels = c("low", "moderate", "high"))) %>%
  mutate(bacteremia_time_5yr_censored = case_when(bacteremia_days_since_ce > 365.25*5 ~ 365.25*5, 
                                                  bacteremia_days_since_ce >= 0 & bacteremia_days_since_ce <= 365.25*5 ~ bacteremia_days_since_ce,death_days_since_ce <= 365.25*5 ~ death_days_since_ce, #if die before 5 years, mark as censored 
                                                  is.na(bacteremia_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*5, na.rm = TRUE),
                                               TRUE ~ NA_real_)) %>%
  mutate(bacteremia_5yr = case_when(is.na(bacteremia_days_since_ce) & bacteremia_time_5yr_censored < 365.25*5 ~ 0, 
                                bacteremia_time_5yr_censored >= 0 & bacteremia_time_5yr_censored < 365.25*5 & !is.na(bacteremia_days_since_ce)~ 1, 
                                bacteremia_time_5yr_censored >= (365.25 * 5) ~0)) %>%
  mutate(graft_manip = if_else(is.na(graft_manip), "Unknown", graft_manip)) %>%
  mutate(match_status_reduced = factor(if_else(match_status %in% c("10/10", "9/10", "8/10", "Haploidentical"), "Non Cord Blood", match_status), 
                              levels = c("Haploidentical", "8-10/10", "Cord Blood", "Unknown"))) %>%
  mutate(disease_relapse = if_else(is.na(graft_fail_date), "No", "Yes")) 

cohort_covars %>% output_tbl("survival_analysis_dataset")
# by default, the survival function uses the Kaplan-Meier method
# outcome must be TRUE/FALSE
# Surv( aim_2a_2$death_days_since_ce, aim_2a_2$survival)

# we need a last-day of follow up for censoring
# among the survival ones, which ones were censored? 
# 1-yr survival 
# print(survfit2(Surv(survival_time_1yr_censored, survival_1yr) ~ ferritin_pre, 
#         data = aim_2a_2 %>% filter(!is.na(survival_time_1yr_censored))) %>% 
#   ggsurvfit() +
#   labs(
#     x = "Days",
#     y = "Overall 1yr-survival probability"
#     ) + 
#   add_confidence_interval() +
#   add_risktable()) # the numbers at risk

# 3-yr survival
# print(survfit2(Surv(survival_time_3yr_censored, survival_3yr) ~ ferritin_pre, 
#         data = aim_2a_2 %>% filter(!is.na(survival_time_3yr_censored))) %>% 
#   ggsurvfit() +
#   labs(
#     x = "Days",
#     y = "Overall 3yr-survival probability"
#     ) + 
#   add_confidence_interval() +
#   add_risktable()) # the numbers at risk 

# 5-yr survival
# ggsurvfit::survfit2() returns a survfit object
# survival::survfit()



