# survival analysis
# source("code/EDA.R")
# install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
require(ggplot2)
require(lubridate)
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
                graft_fail_days_since_ce = as.numeric(difftime(graft_fail_date, transplant_date, units = "days"))) %>%
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
         graft_fail_days_since_ce_5yr_censored = case_when(graft_fail_days_since_ce > 365.25*5 ~ 365.25*5, #censored at 5 years
                                                           graft_fail_days_since_ce <= 365.25*5 ~ graft_fail_days_since_ce, #graft failure
                                                           death_days_since_ce <= 365.25*3 ~ death_days_since_ce, #death: competing event
                                                           last_in_person_visit_since_ce <= 365.25*5 & is.na(graft_fail_days_since_ce) ~ last_in_person_visit_since_ce, #censored at last visit
                                                           last_in_person_visit_since_ce > 365.25*5 & is.na(graft_fail_days_since_ce) ~ 365.25*5, #censored at 5 years follow-up
                                                            TRUE ~ NA),
         # graft failure is a competing event: 0: censored, 1: graft failure, 2:death                                       
         graft_fail_5yr = case_when(graft_fail_days_since_ce <= 365.25*5 ~ "1",
                                    death_days_since_ce <= 365.25*3 ~ "2", 
                                    graft_fail_days_since_ce > 365.25*5 | !is.na(last_in_person_visit_since_ce) ~ "0",
                                    TRUE ~ NA))

cohort_covars <- cohort_covars %>% #mutate(ferritin_pre = factor(ferritin_pre, levels = c("low", "moderate", "high")),
         #ferritin_pre_max = factor(ferritin_pre_max, levels = c("low", "moderate", "high"))) %>%
  mutate(bacteremia_time_5yr_censored = case_when(bacteremia_days_since_ce > 365.25*5 ~ 365.25*5, 
                                                  bacteremia_days_since_ce >= 0 & bacteremia_days_since_ce <= 365.25*5 ~ bacteremia_days_since_ce,death_days_since_ce <= 365.25*5 ~ death_days_since_ce, #if die before 5 years, mark as censored 
                                                  is.na(bacteremia_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*5, na.rm = TRUE),
                                               TRUE ~ NA_real_)) %>%
  mutate(bacteremia_5yr = case_when(is.na(bacteremia_days_since_ce) & bacteremia_time_5yr_censored < 365.25*5 ~ 0, 
                                bacteremia_time_5yr_censored >= 0 & bacteremia_time_5yr_censored < 365.25*5 & !is.na(bacteremia_days_since_ce)~ 1, 
                                bacteremia_time_5yr_censored >= (365.25 * 5) ~0)) %>%
  mutate(graft_manip = if_else(is.na(graft_manip), "Unknown", graft_manip)) %>%
  mutate(disease_relapse = if_else(is.na(graft_fail_date), "No", "Yes")) 

cohort_covars %>% output_tbl("survival_analysis_dataset")

# results_tbl("cr_data") %>% collect() %>% output_tbl("cr_data", local = TRUE, file = TRUE)
                              
# stat_dataset %>% filter(chart_completion, eligibility == "Yes") %>%
#       group_by(survival_5yr, transplant_type_combined) %>% summarise(n = n_distinct(record_id)) %>% view()

# now get ferritin time series for aim 3
# there are patients who only had 1 ferritin measurement after transplant
# rslt$potential_aim3 <- results_tbl("analytics_dataset") %>%
#       filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date, IRT) %>%
#       left_join(results_tbl("covar_ferritin_mx"), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>% 
#       filter(transplant_date - ferritin_date <= 365) %>%
#       group_by(person_id) %>% 
#       mutate(n = n()) %>% filter(n >=1) %>% ungroup()

# # we want the base ferritin levels to be closet to the start date of treatment
# rslt$phleb_px <- results_tbl("analytics_dataset") %>%
#       # filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date) %>%
#       find_procedures(procedure_codeset_name = "therapeutic_phleb_px") %>%
#       filter(procedure_date >= transplant_date) %>%
#       filter(procedure_date - transplant_date <= 365) %>%
#       group_by(person_id) %>%
#       collect_new() %>%
#       mutate(start_date = min(procedure_date, na.rm = TRUE), 
#              end_date = max(procedure_date, na.rm = TRUE),
#              duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration) %>% 
#       mutate(route = NA) %>%
#       mutate(IRT_type = "phlebotomy") 

# # patients who got chelation
# dx_codeset<- load_codeset("defibrotide_rx") %>% mutate(type = "defibrotide") %>% #defibrotide
#                 union(load_codeset("deferoxamine_rx") %>% mutate(type = "deferoxamine")) %>% #deferoxamine
#                 union(load_codeset("deferasirox_rx") %>% mutate(type = "deferasirox")) %>% #deferasirox
#                 compute_new(temp = TRUE, name = "drug_id")

# rslt$chelation_rx <- results_tbl("analytics_dataset") %>%
#       # filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date) %>%
#       find_drugs(dx_codeset) %>%
#       filter(drug_exposure_start_date >= transplant_date) %>%
#       filter(drug_exposure_start_date - transplant_date <= 365) %>% 
#       distinct(person_id, drug_exposure_start_date, IRT_type = type, .keep_all = TRUE) 

# # for intravenous drugs, assume a single one per day
# rslt$chelation_iv <- rslt$chelation_rx %>% filter(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous")) %>%
#       group_by(person_id) %>%
#       collect_new() %>%
#       mutate(start_date = min(drug_exposure_start_date, na.rm = TRUE), 
#              end_date = max(drug_exposure_start_date, na.rm = TRUE),
#              duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration, .keep_all = TRUE) %>% 
#       mutate(route = "IV") %>%
#       select(person_id, record_id, start_date, end_date, transplant_date, duration, IRT_type, route)

# # oral drugs
# rslt$chelation_oral <- rslt$chelation_rx %>% filter(!(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous"))) %>%
#       filter(!is.na(quantity)) %>%
#       collect_new() %>%
#       group_by(person_id) %>%
#       mutate(start_date = min(drug_exposure_start_date, na.rm = TRUE), 
#              end_date = drug_exposure_start_date + days(quantity),
#              duration = quantity) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration, .keep_all = TRUE) %>% 
#       select(person_id, record_id, start_date, end_date, transplant_date, duration, IRT_type, route = route_source_value) 

# # Only took the duration of the first prescription, some prescriptions have up to 11 refills 

# # append the chelation and IRT patients
# rslt$irt <- rslt$phleb_px %>%
#           union(rslt$chelation_iv) %>%
#           union(rslt$chelation_oral) 

# # who got more than 1 treatment
# rslt$irt_exclude <- rslt$irt %>% group_by(record_id) %>% 
#       summarise(n = n()) %>% 
#       filter(n> 1) %>% ungroup() %>% view()
      
# rslt$irt_exclude %>% 
#       inner_join(rslt$irt, by = "record_id") %>% 
#       arrange(person_id) %>% view()
      
# # exclude patients with more than 1 treatment for now: 
# rslt$irt <- rslt$irt %>% anti_join(rslt$irt_exclude, by = "record_id")

# rslt$irt %>% distinct_ct()

# # ferritin_0 is the one closest to the start date of treatment
# rslt$ferritin_start <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>%
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       # filter(transplant_date - ferritin_date <= 365) %>%
#       group_by(person_id) %>% 
#       filter(start_date <= ferritin_date) %>%
#       slice_max(start_date, n =1, with_ties = FALSE) %>% ungroup() %>%
#       select(person_id, record_id, transplant_date, start_date, ferritin0_date = ferritin_date, ferritin0 = ferritin) 
      

# rslt$ferritin_end <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>%
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       filter(end_date <= ferritin_date) %>%
#       group_by(person_id) %>% 
#       slice_min(end_date, n =1, with_ties = FALSE) %>% ungroup() %>%
#       select(person_id, record_id, transplant_date, end_date, ferritinf_date = ferritin_date, ferritinf = ferritin)

# rslt$ferritin_all <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>% 
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>%
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       filter(transplant_date - ferritin_date <= 365) %>% collect_new()

# print(rslt$ferritin_all %>% 
#       mutate(ferritin_date = as.numeric(ferritin_date - transplant_date)) %>%
#       mutate(start_date = as.numeric(start_date - transplant_date)) %>%
#       mutate(end_date = as.numeric(end_date - transplant_date)) %>% 
#       filter(duration > 0) %>%
#       ggplot(aes(x = ferritin_date, y = ferritin, group = record_id)) + 
#       geom_point() + 
#       geom_line() +
#       geom_vline(aes(xintercept = start_date), linetype = "dashed", color = "red") +
#       geom_vline(aes(xintercept = end_date), linetype = "dashed", color = "red") +
#       facet_wrap(~record_id, ncol = 3)) 

# # special cases: patients with only 1 treatment: duration == 0
# rslt$ferritin_start %>% inner_join(rslt$ferritin_end, by = c("record_id", "person_id", "transplant_date")) %>% 
#             mutate(duration = end_date - start_date) %>% view()


