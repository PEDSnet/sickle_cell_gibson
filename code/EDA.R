
library(tidyr)
library(lubridate)

library(pdftools)
library(ggplot2)

rslt = list()

rslt$study_cohorts <- results_tbl("study_cohorts") %>% select(person_id, aim_2a_1:aim_3_2) %>%
                    union(results_tbl("multi_transplant_cohort") %>% anti_join(results_tbl("study_cohorts"), by = "person_id") %>% 
                            select(person_id, aim_2a_1:aim_3_2)) %>%
                    collect()
study_cohorts <- results_tbl("study_cohorts") %>% collect() %>% select(person_id, site, gender)
                        left_join(results_tbl("cr_cohorts") %>% 
                                select(person_id, record_id, transplant_date, transplant_type, second_transplant_date) %>% collect(), by = "person_id")

# get iron overload status
# majority with units nanogram per milliliter
# assume missing units are nanogram per milliliter, no unit conversion done  
cohort_covars <- results_tbl("covar_ferritin_mx") %>% 
            filter((aim == "aim_2a_2" & ferritin_type == "pre") |
                    ((aim == "aim_2b_2" | aim == "aim_3_1" | aim == "aim_3_2") & ferritin_type == "post")|
                     aim == "aim_2a_1" | aim == "aim_2b_1",
                    !is.na(ferritin_level)) %>% collect() %>%
            group_by(person_id, transplant_date, aim) %>% 
            slice_min(abs(difftime(ferritin_date, transplant_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup() %>% 
            left_join(results_tbl("no_multi_transplant_px") %>% 
                        group_by(person_id) %>%
                        summarise(no_transplants = n_distinct(transplant_date),
                                disease_relapse_date = max(transplant_date, na.rm = TRUE)) %>%
                        filter(no_transplants > 1) %>% 
                        ungroup() %>% mutate(disease_relapse = TRUE) %>% collect(), by = "person_id") %>%
            mutate(disease_relapse = if_else(is.na(disease_relapse), FALSE, disease_relapse))

cohort_covars %>% view()
# a patient could be included in multiple aims
cohort_covars %>% distinct(person_id, aim) %>% count() # n = 1159
cohort_covars %>% distinct_ct() # n = 723
    
# overall survival
# there are 2 patients with more than 1 death causes
# for patients with more than 1 transplant, death days since ce is the number of days since 
# the second transplants
cohort_covars <- cohort_covars %>% copy_to_new(df = ., name = "sdsd") %>%
                        left_join(cdm_tbl("death") %>% 
                                select(person_id, death_date), by = "person_id") %>% collect() %>%
                        mutate(death = if_else(!is.na(death_date), TRUE, FALSE), 
                                death_days_since_ce = if_else(!is.na(death_date), 
                                                    as.numeric(difftime(death_date, transplant_date, units = "days")), NA)) %>%
                        distinct(person_id, aim, death_date, .keep_all = TRUE)
cohort_covars <- cohort_covars %>% mutate(death_days_since_last_transplant = if_else(disease_relapse, 
                                        as.numeric(difftime(death_date, disease_relapse_date, units = "days")),
                                        as.numeric(difftime(death_date, transplant_date, units = "days"))))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, death_date) %>% count() # n = 1159

# GVHD status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_gvhd_dx") %>% collect() %>%
                                mutate(gvhd_days_since_ce = as.numeric(difftime(condition_start_date, transplant_date, units = "days"))) %>%
                                group_by(person_id, transplant_date) %>%
                                filter(gvhd_days_since_ce >= 0) %>%
                                slice_min(abs(gvhd_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                ungroup() %>% select(person_id, gvhd_start_date = condition_start_date, transplant_date, gvhd_days_since_ce), by = c("person_id", "transplant_date")) %>%
                        mutate(GVHD = if_else(!is.na(gvhd_start_date), TRUE, FALSE)) 
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, gvhd_days_since_ce) %>% count() # n = 1159

# VOD status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_vod_dx") %>%
                                    mutate(vod_date = pmax(bilirubin_date, abd_ultrasound_date, defibrotide_date, na.rm = TRUE)) %>%
                                    collect() %>%
                                mutate(vod_days_since_ce = as.numeric(difftime(vod_date, transplant_date, units = "days"))) %>%
                                group_by(person_id, transplant_date) %>%
                                filter(vod_days_since_ce >= 0) %>%
                                slice_min(abs(vod_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                ungroup() %>% 
                                select(person_id, vod_date, transplant_date, vod_days_since_ce), by = c("person_id", "transplant_date")) %>%
                        mutate(VOD = if_else(!is.na(vod_date), TRUE, FALSE))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, vod_date) %>% count() # n = 1159

# bacteremia status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_bacteremia_dx") %>% collect() %>%
                                    mutate(bacteremia_days_since_ce = as.numeric(difftime(bacteremia_date, transplant_date, units = "days")),
                                           bacteremia = organism_concept_name) %>%
                                    group_by(person_id, transplant_date) %>%
                                    filter(bacteremia_days_since_ce >= 0) %>%
                                    slice_min(abs(bacteremia_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                    ungroup() %>% 
                                    select(person_id, bacteremia_date, transplant_date, bacteremia_days_since_ce, bacteremia), by = c("person_id", "transplant_date"))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, bacteremia_date) %>% count() # n = 1159

# Immune reconstitution
# results_tbl("covar_CD348_mx") %>% union(results_tbl("covar_IgM_mx") %>% 
#                                         mutate(CDtype = "IgM")) %>% 
#                                 select(person_id, transplant_date, 
#                                         immune_date = measurement_date,
#                                         immune_type = CDtype,
#                                             measurement_concept_id, measurement_concept_name, unit_concept_name, range_high, range_low,
#                                         immune_val = value_as_number)

# cases when patients received both phlebotomy and chelation
cohort_covars <- cohort_covars %>% mutate(val = TRUE) %>% pivot_wider(names_from = aim, values_from = val, values_fill = FALSE) %>%
                mutate(IRT = case_when( aim_3_1 & aim_3_2 ~ "phlebotomy & chelation",
                                        aim_3_1 ~ "phlebotomy", 
                                        aim_3_2 ~ "chelation", 
                                       !aim_3_1 & !aim_3_2 ~ "No IRT", 
                                       TRUE ~ NA)) %>% 
                pivot_longer(cols = aim_2a_2:aim_3_1, names_to = "aim", values_to = "val") %>%
                filter(val) %>% select(-val)
cohort_covars %>%view()
# cohort_covars %>% group_by(person_id) %>% summarise(n = n_distinct(IRT)) %>% filter(n >1) %>%
#                 inner_join(cohort_covars, by = "person_id") %>% 
#                 select(person_id, IRT, aim) %>% 
#                 arrange(person_id, IRT) %>% view()
# cohort_covars %>% filter(!(duplicated(person_id) & IRT == "No IRT")) %>% 
#                 select(person_id, transplant_date, IRT, aim) %>%  arrange(person_id, IRT) %>%  view()
                # filter(!(duplicated(person_id, transplant_date) & (IRT == "phlebotomy" | IRT == "chelation"))) %>% view()

message("Compute platelet engraftment dates")
