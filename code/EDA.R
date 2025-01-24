
library(tidyr)
library(lubridate)

library(pdftools)
library(ggplot2)

rslt = list()

 
# get cohort along with dates of 2nd transplants (if any)
# For patients with available chart reviews, use chart reivew transplant dates as index dates
# Filter out false positives from chart reviews, only include cr patients with completed cr
# only include the first LIC measurement for each patient
LIC_cutoffs <- c(5, 10)
LIC_cutoffs_binary <- c(8)
no_LIC <- 3
keep_ties_LIC <- TRUE

# extract most 3 recent LIC values                                 
rslt$cr_data_LIC_recent_3 <- results_tbl("cr_data") %>% collect() %>% 
                                val_extraction(cohort = ., 
                                                no_value = no_LIC, grouping_id = record_id, 
                                                value_name = LIC, value_date = LIC_date, 
                                                cutoffs = LIC_cutoffs,
                                                value_type = "LIC", keep_tie = keep_ties_LIC, 
                                                slice_by = "most_recent")   

# extract most 3 recent LIC values                                 
rslt$cr_data_LIC_recent_1 <- results_tbl("cr_data") %>% collect() %>% 
                                val_extraction(cohort = ., 
                                                no_value = 1, grouping_id = record_id, 
                                                value_name = LIC, value_date = LIC_date, 
                                                cutoffs = LIC_cutoffs,
                                                value_type = "LIC", keep_tie = FALSE, 
                                                slice_by = "most_recent")   

# extract highest LIC values 
rslt$cr_data_LIC_max <- results_tbl("cr_data") %>% collect() %>% 
                                val_extraction(cohort = ., 
                                                no_value = no_LIC, grouping_id = record_id, 
                                                value_name = LIC, value_date = LIC_date, 
                                                cutoffs = LIC_cutoffs,
                                                value_type = "LIC", keep_tie = keep_ties_LIC, 
                                                slice_by = "max") 
# summary of LIC values
LIC_ct <- results_tbl("cr_data") %>% filter(!is.na(LIC_type)) %>%
                        group_by(record_id, LIC_type, transplant_date) %>% 
                        summarise(no_LIC = n(), no_LIC_10 = sum(ifelse(is.na(LIC) | LIC <=10, 0, 1))) %>% ungroup() %>%
                        pivot_wider(names_from = LIC_type, values_from = c("no_LIC", "no_LIC_10")) %>% 
                        rename(no_LIC_pre = `no_LIC_pre-transplant`, no_LIC_post = `no_LIC_post-transplant`,
                               no_LIC_over_10_pre = `no_LIC_10_pre-transplant`,
                               no_LIC_over_10_post = `no_LIC_10_post-transplant`) %>%
                        mutate(across(c("no_LIC_pre", "no_LIC_post", "no_LIC_over_10_pre", "no_LIC_over_10_post"), ~ifelse(is.na(.x), 0, .x))) %>% collect()

rslt$LIC <- rslt$cr_data_LIC_recent_3 %>% 
                        left_join(rslt$cr_data_LIC_recent_1 %>% select(record_id, LIC_1_level, LIC_1, LIC_1_type = LIC_type, 
                                                                LIC_days_since_transplant_1, 
                                                                transplant_date), by = c("record_id", "LIC_type" = "LIC_1_type", "transplant_date")) %>%
                        left_join(rslt$cr_data_LIC_max %>% select(record_id, LIC_max, LIC_3_max, 
                                                                LIC_days_since_transplant_max, LIC_days_since_transplant_3_max,
                                                                LIC_type_max,
                                                                LIC_max_level, LIC_3_max_level, transplant_date), by = c("record_id", "LIC_type" = "LIC_type_max", "transplant_date")) %>%
                        mutate(LIC_level_3_binary = if_else(LIC_3 < LIC_cutoffs_binary[1], "low", "high"),
                                LIC_level_1_binary = if_else(LIC_1 < LIC_cutoffs_binary[1], "low", "high"),
                                LIC_level_3_max_binary = if_else(LIC_3_max < LIC_cutoffs_binary[1], "low", "high"),
                                across(c("LIC_level_3_max_binary", "LIC_level_3_binary", "LIC_level_1_binary"), ~factor(.x, levels = c("low", "high")))) 

rslt$LIC %>% output_tbl("cr_LIC_data")                          

rslt$cr_data <- results_tbl("cr_data") %>% select(record_id, eligibility, transplant_date, 
                                                        transplant_type, donor_relation, match_status, graft_manip, manip_type,
                                                        graft_fail, graft_fail_date, second_transplant_date, chart_completion) %>%  
                                        collect() %>% distinct() %>% 
                        left_join(LIC_ct, by = c("record_id", "transplant_date")) %>%
                        rename(transplant_date_cr = transplant_date, 
                                transplant_type_cr = transplant_type, 
                                second_transplant_date_cr = second_transplant_date,
                                disease_relapse = graft_fail) %>% copy_to_new(df = ., name = "cr_data")             
                                
rslt$study_cohorts <- results_tbl("study_cohorts") %>% 
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        left_join(rslt$cr_data, by = "record_id") %>% 
                        mutate(transplant_date_consistency = ifelse((transplant_date_cr <= transplant_date + days(3) & transplant_date_cr <= transplant_date) | (transplant_date_cr >= transplant_date - days(3) & transplant_date_cr >= transplant_date), 1, 0),
                                transplant_date = ifelse(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                                transplant_type = ifelse(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                                chart_completion = ifelse(is.na(chart_completion), FALSE, chart_completion)) %>%
                        filter((eligibility == "Yes" & chart_completion) | !chart_completion) %>%
                        select(-transplant_date_cr, -transplant_type_cr)

# this check should returns NULL which means all 1st LIC measurements have valid methods and dates
# results_tbl("cr_data") %>% anti_join(results_tbl("cr_data") %>% filter(no_LIC_measurement == 1 | is.na(no_LIC_measurement)) , by = "record_id")

# get second transplant date from EHR
rslt$transplant_px <- rslt$study_cohorts %>% inner_join(results_tbl("no_multi_transplant_px") %>% 
                        select(person_id, second_transplant_date = transplant_date), by = "person_id") %>% 
                        filter(second_transplant_date >= transplant_date) %>% 
                        mutate(second_transplant_date = if_else(second_transplant_date == transplant_date, NA, second_transplant_date)) %>%
                        group_by(person_id) %>%
                        slice_min(second_transplant_date, n = 1, with_ties = FALSE, na_rm = FALSE) %>% 
                        ungroup()

# use second transplant dates from chart review when available
rslt$transplant_px <- rslt$transplant_px %>% 
                        mutate(second_transplant_date = case_when(!is.na(second_transplant_date_cr) ~ second_transplant_date_cr, 
                                                                (is.na(second_transplant_date_cr) & chart_completion) ~ NA, #overwrite second transplant dates from EHR
                                                                TRUE ~ second_transplant_date)) %>% 
                        distinct(record_id, second_transplant_date)

rslt$study_cohorts <- rslt$study_cohorts %>% left_join(rslt$transplant_px, by = "record_id") %>% 
                        group_by(person_id) %>%
                        slice_min(transplant_date, n = 1, with_ties = FALSE, na_rm = FALSE) %>% ungroup()

# get FERRITIN iron overload status, only for the 1st transplant
# majority with units nanogram per milliliter
# assume missing units are nanogram per milliliter, no unit conversion done  

ferritin_cutoff1 <- 1000
ferritin_cutoff2 <- 2500

# missing pre-ferritin values: 23839237 and 4988181

# let's take the average of the 3 most-recent ferritins 
no_ferritin <- 3
keep_tie <- FALSE

rslt$ferritin = results_tbl("covar_ferritin_mx")
# most 3 recent ferritin values 
ferritin_pre <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "pre", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)
# count number of ferritin measurements before transplant
ferritin_pre_ct <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                        left_join(rslt$ferritin, by = "person_id") %>% 
                        filter(transplant_date > ferritin_date) %>% 
                        filter(transplant_date - ferritin_date <= 365) %>%
                        group_by(person_id, transplant_date) %>% 
                        summarise(no_ferritin_pre = n()) %>% ungroup() %>% collect()

ferritin_post <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "post", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)

# count number of ferritin measurements after transplant
ferritin_post_ct <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                        left_join(rslt$ferritin, by = "person_id") %>% 
                        filter(transplant_date <= ferritin_date) %>% 
                        filter(ferritin_date - transplant_date <= 365) %>%
                        group_by(person_id, transplant_date) %>% 
                        summarise(no_ferritin_post = n()) %>% ungroup() %>% collect()

# slice by highest values 
ferritin_pre_max <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "pre_max", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)

ferritin_post_max <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "post_max", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)


ferritin <- ferritin_pre %>% select(person_id, transplant_date, 
                                ferritin_pre = ferritin_level, 
                                ferritin_days_pre = ferritin_days_since_transplant) %>%
                full_join(ferritin_post %>% select(person_id, transplant_date, 
                                                ferritin_post = ferritin_level, 
                                                ferritin_days_post = ferritin_days_since_transplant), by =  c("person_id", "transplant_date")) %>%
                full_join(ferritin_pre_max %>% select(person_id, transplant_date, 
                                                ferritin_pre_max = ferritin_level, 
                                                ferritin_days_pre_max = ferritin_days_since_transplant), by = c("person_id", "transplant_date")) %>%
                full_join(ferritin_post_max %>% select(person_id, transplant_date, 
                                                ferritin_post_max = ferritin_level, 
                                                ferritin_days_post_max = ferritin_days_since_transplant), by = c("person_id", "transplant_date")) %>%
                mutate(across(c("ferritin_post", "ferritin_pre", "ferritin_post_max", "ferritin_pre_max"), ~factor(.x, levels = c("low", "moderate", "high")))) %>%
                left_join(ferritin_pre_ct, by = c("person_id", "transplant_date")) %>%
                left_join(ferritin_post_ct, by = c("person_id", "transplant_date"))

ferritin %>% output_tbl("cr_ferritin_data")

rslt$study_cohorts <- rslt$study_cohorts %>% #select(person_id, transplant_date, site, aim_2a_2, aim_2b_2,
                                                # record_id, site_id, birth_date, gender, race, ethnicity,
                                                # scd_concept_name, scd_type, transplant_concept_name, transplant_type,
                                                # multi_transplants, eligibility, donor_relation, match_status,
                                                # graft_manip, manip_type, disease_relapse, graft_fail_date, second_transplant_date_cr,
                                                # chart_completion, transplant_date_consistency) %>% 
                        left_join(ferritin %>% select(person_id, transplant_date, no_ferritin_pre, no_ferritin_post) %>% 
                                        copy_to_new(df = ., name = "sdsd"), by = c("person_id", "transplant_date"))

# cohort_covars %>% view()
# patients that belong to aim_2a_2 with missing pre-ferritin should be excluded from that aim
# no patients not belong to any aims 
# rm("ferritin_pre", "ferritin_post", "ferritin_pre_max", "ferritin_post_max")
cohort_covars <- rslt$study_cohorts %>% collect() %>% 
                                mutate(aim_2a_2 = if_else(is.na(no_ferritin_pre) | no_ferritin_pre == 0, FALSE, TRUE)) %>%
                                mutate(aim_2b_2 = if_else(is.na(no_ferritin_post) | no_ferritin_post == 0, FALSE, FALSE)) %>%
                                mutate(aim_2a_1 = if_else(is.na(no_LIC_pre) | no_LIC_pre == 0, FALSE, TRUE)) %>%
                                mutate(aim_2b_1 = if_else(is.na(no_LIC_post) | no_LIC_post == 0, FALSE, TRUE)) %>%
                                mutate(total_aims = rowSums(select(., starts_with("aim_")))) %>%
                                filter(total_aims > 0) %>% select(-total_aims) %>% copy_to_new(df = ., name = "cohort_covars")
                
# cohort_covars %>% view()

# a patient could be included in multiple aims
# cohort_covars %>% distinct(person_id, aim) %>% count() # n = 1159
# cohort_covars %>% distinct_ct() # n = 820

# last in-person visit date
last_visits <- cdm_tbl("visit_occurrence") %>% inner_join(cohort_covars %>% 
                                                        select(person_id, transplant_date) %>%
                                                        copy_to_new(df =., name = "sdfsd"), by = "person_id") %>%
            select(person_id, visit_occurrence_id, visit_start_date, visit_end_date, visit_concept_id, transplant_date) %>%
            filter(visit_start_date >= transplant_date) %>% collect() %>%
            mutate(last_visit_since_ce = as.numeric(difftime(visit_end_date, transplant_date, units = "days"))) 
            
last_visits_in_person <- last_visits %>% filter(!(visit_concept_id %in% c(44814711, 44814653, 44814649, 44814650))) %>% 
                        group_by(person_id) %>%
                        slice_max(last_visit_since_ce, n = 1, with_ties = FALSE) %>% ungroup()          
            
last_visits <- last_visits %>% group_by(person_id) %>%
               slice_max(last_visit_since_ce, n = 1, with_ties = FALSE) %>% ungroup() %>% 
               left_join(last_visits_in_person %>% select(person_id, 
                                                        last_in_person_visit_since_ce = last_visit_since_ce), by = "person_id") %>%
                select(person_id, last_visit_since_ce, last_in_person_visit_since_ce)
# last_visits %>% view()
cohort_covars <- cohort_covars %>% left_join(last_visits %>% copy_to_new(df =., name = "sdfsd"), by = "person_id")

# overall survival
# there are 2 patients with more than 1 death causes
# for patients with more than 1 transplant, death days since ce is the number of days since the 1st transplants
cohort_covars <- cohort_covars %>% left_join(cdm_tbl("death") %>% 
                                select(person_id, death_date), by = "person_id") %>% collect() %>%
                                mutate(survival = if_else(!is.na(death_date), FALSE, TRUE), 
                                death_days_since_ce = if_else(!is.na(death_date), as.numeric(difftime(death_date, transplant_date, units = "days")), NA)) %>%
                                group_by(person_id) %>%
                                slice_min(abs(death_days_since_ce), n = 1, with_ties = FALSE) %>% ungroup()

# cohort_covars %>% distinct_ct() # n = 820
# cohort_covars %>% distinct(person_id, death_date) %>% count() # n = 820

# GVHD status
gvhd <- results_tbl("covar_gvhd_dx") %>% select(person_id, gvhd_start_date = condition_start_date) %>% 
        inner_join(rslt$study_cohorts %>% select(person_id, transplant_date), by = "person_id") %>% collect() %>%
        mutate(gvhd_days_since_ce = as.numeric(difftime(gvhd_start_date, transplant_date, units = "days"))) %>%
        filter(gvhd_days_since_ce >= 0, gvhd_days_since_ce <= 365) %>%
        group_by(person_id) %>%
        slice_min(gvhd_days_since_ce, n = 1, with_ties = FALSE) %>% ungroup() %>%
        mutate(GVHD = TRUE)

cohort_covars <- cohort_covars %>% left_join(gvhd, by = c("person_id", "transplant_date")) %>%
                        mutate(GVHD = if_else(!is.na(gvhd_start_date), TRUE, FALSE)) 
# cohort_covars %>% distinct_ct() # n = 723
# cohort_covars %>% distinct(person_id, aim, gvhd_days_since_ce) %>% count() # n = 1159

# VOD status, n = 203 already checked for distinct person_id
# using the 90 days window for VOD
vod <- results_tbl("covar_vod_dx") %>% 
        left_join(rslt$study_cohorts %>% select(person_id, transplant_date), by = "person_id") %>%
        mutate(vod_date := pmin(min_defi_abd_date, min_bili_abd_date, min_bili_abd_weight_date, na.rm = TRUE)) %>%
        collect() 
        
vod_90 <- vod %>% filter((as.numeric(difftime(vod_date, transplant_date, units = "days")) <= 90 &
                                vod_date >= transplant_date) | 
                                abs(as.numeric(difftime(abd_date, transplant_date, units = "days"))) <= 7) %>% 
        distinct(person_id, vod_cat_1, vod_cat_2, vod_cat_3, vod_cat_4, vod_cat_5, vod) 
vod_cat2_ids <- vod_90 %>% filter(vod_cat_2 == 2, is.na(vod_cat_1), is.na(vod_cat_3), is.na(vod_cat_4), is.na(vod_cat_5))
vod_90 <- vod_90 %>% anti_join(vod_cat2_ids, by = "person_id") %>%
            rename("vod_cat_1_90" = "vod_cat_1", "vod_cat_2_90" = "vod_cat_2", "vod_cat_3_90" = "vod_cat_3",
                  "vod_cat_4_90" = "vod_cat_4",
                  "vod_cat_5_90" = "vod_cat_5", "vod_90" = "vod")

# using the 30 days window for VOD
vod_30 <- vod %>% filter((as.numeric(difftime(vod_date, transplant_date, units = "days")) <= 30 &
                                vod_date >= transplant_date) | 
                                abs(as.numeric(difftime(abd_date, transplant_date, units = "days"))) <= 7) %>% 
        distinct(person_id, vod_cat_1, vod_cat_2, vod_cat_3, vod_cat_4, vod_cat_5, vod) 

vod_cat2_ids_30 <- vod_30 %>% filter(vod_cat_2 == 2, is.na(vod_cat_1), is.na(vod_cat_3), is.na(vod_cat_4), is.na(vod_cat_5))
vod_30 <- vod_30 %>% anti_join(vod_cat2_ids_30, by = "person_id")


vod %>% distinct_ct()
# vod %>% arrange(person_id) %>% view()        
cohort_covars <- cohort_covars %>% left_join(vod_90, by = c("person_id")) %>% 
                        mutate(vod_90 = if_else(!is.na(vod_90), TRUE, FALSE)) %>%
                        left_join(vod_30 #%>% 
                              #   rename(c("vod_cat_1_30" = "vod_cat_1",
                              #   "vod_cat_2_30" = "vod_cat_2",
                              #   "vod_cat_3_30" = "vod_cat_3",
                              #   "vod_cat_4_30" = "vod_cat_4",
                              #   "vod_cat_5_30" = "vod_cat_5","vod_30" = "vod",))
                              , by = c("person_id")) %>% 
                        mutate(vod = if_else(!is.na(vod), TRUE, FALSE)) 

# cohort_covars %>% distinct_ct() # n = 723
# cohort_covars %>% distinct(person_id, aim, vod_date) %>% count() # n = 1159

# bacteremia status
bacteremia <- rslt$study_cohorts %>% select(person_id, transplant_date) %>%
                inner_join(results_tbl("covar_bacteremia_dx"), by = "person_id") %>% collect() %>%
                mutate(bacteremia_days_since_ce = as.numeric(difftime(bacteremia_date, transplant_date, units = "days")),
                        bacteremia = TRUE,
                        bacteria = organism_concept_name) %>%
                filter(bacteremia_days_since_ce >= 0, bacteremia_days_since_ce <= 360) %>%
                group_by(person_id, bacteremia) %>%
                summarise(bacteremia_days_since_ce = min(abs(bacteremia_days_since_ce), na.rm = TRUE)) %>% ungroup() 

cohort_covars <- cohort_covars %>% left_join(bacteremia, by = c("person_id"))
# cohort_covars %>% distinct_ct() # n = 723
# cohort_covars %>% distinct(person_id, aim, bacteremia_date) %>% count() # n = 1159

# Immune reconstitution
# results_tbl("covar_CD348_mx") %>% union(results_tbl("covar_IgM_mx") %>% 
#                                         mutate(CDtype = "IgM")) %>% 
#                                 select(person_id, transplant_date, 
#                                         immune_date = measurement_date,
#                                         immune_type = CDtype,
#                                             measurement_concept_id, measurement_concept_name, unit_concept_name, range_high, range_low,
#                                         immune_val = value_as_number)

# cases when patients received both phlebotomy and chelation
cohort_covars <- cohort_covars %>% #mutate(val = TRUE) %>% pivot_wider(names_from = aim, values_from = val, values_fill = FALSE) %>%
                mutate(IRT = case_when( aim_3_1 & aim_3_2 ~ "phlebotomy & chelation",
                                        aim_3_1 ~ "phlebotomy", 
                                        aim_3_2 ~ "chelation", 
                                       !aim_3_1 & !aim_3_2 ~ "No IRT", 
                                       TRUE ~ NA))

# platelet engraftment dates
# ANC engraftment dates

# immune reconstitution
CD3_threshold <- 500
CD3_reconstitute <- get_immune_reconstitution(threshold = CD3_threshold, 
                                      cohort_tbl = cohort_covars, 
                                      immune_tbl = results_tbl("covar_CD3_mx"),
                                      immune_var = CD3_ct) %>%
                    select(person_id, CD3_reconstitute_3mon = reconstitute_3mon, 
                            CD3_reconstitute_6mon = reconstitute_6mon, 
                            CD3_reconstitute_9mon = reconstitute_9mon, 
                            CD3_reconstitute_12mon = reconstitute_12mon, 
                            CD3_ct, CD3_measurement_date = measurement_date)
# CD3_reconstitute <- cohort_covars %>% left_join(results_tbl("covar_CD3_mx") %>% 
#                                 select(person_id, CD3_ct, CD3_date = measurement_date) %>% collect(), by = "person_id") %>%
#                 mutate(CD3_days_since_ce = as.numeric(difftime(CD3_date, transplant_date, units = "days"))) %>%
#                 filter(CD3_days_since_ce >= 0, CD3_days_since_ce <= 365, 
#                         is.na(second_transplant_date) | CD3_days_since_ce <= second_transplant_date,
#                         CD3_ct >= CD3_threshold) %>% 
#                 group_by(person_id) %>%
#                 slice_min(CD3_days_since_ce, n = 1, with_ties = FALSE) %>% ungroup() %>%
#                 mutate(CD3_reconstitute_3mon = if_else(CD3_days_since_ce <= 90, TRUE, FALSE),
#                         CD3_reconstitute_6mon = if_else(CD3_days_since_ce <= 180, TRUE, FALSE),
#                         CD3_reconstitute_9mon = if_else(CD3_days_since_ce <= 270, TRUE, FALSE), 
#                         CD3_reconstitute_12mon = if_else(CD3_days_since_ce <= 365, TRUE, FALSE)) %>%
#                         select(person_id, starts_with("CD3"))

CD4_threshold <- 200
CD4_reconstitute <- get_immune_reconstitution(threshold = CD4_threshold, 
                                      cohort_tbl = cohort_covars, 
                                      immune_tbl = results_tbl("covar_CD4_mx"),
                                      immune_var = CD4_ct) %>%
                    select(person_id, CD4_reconstitute_3mon = reconstitute_3mon, 
                            CD4_reconstitute_6mon = reconstitute_6mon, 
                            CD4_reconstitute_9mon = reconstitute_9mon, 
                            CD4_reconstitute_12mon = reconstitute_12mon, 
                            CD4_ct, CD4_measurement_date = measurement_date)

CD8_threshold <- 200
CD8_reconstitute <- get_immune_reconstitution(threshold = CD8_threshold, 
                                      cohort_tbl = cohort_covars, 
                                      immune_tbl = results_tbl("covar_CD8_mx"),
                                      immune_var = CD8_ct) %>%
                    select(person_id, CD8_reconstitute_3mon = reconstitute_3mon, 
                            CD8_reconstitute_6mon = reconstitute_6mon, 
                            CD8_reconstitute_9mon = reconstitute_9mon, 
                            CD8_reconstitute_12mon = reconstitute_12mon, 
                            CD8_ct, CD8_measurement_date = measurement_date)

IgM_threshold <- 25
IgM_reconstitute <- get_immune_reconstitution(threshold = IgM_threshold, 
                                      cohort_tbl = cohort_covars, 
                                      immune_tbl = results_tbl("covar_IgM_mx"),
                                      immune_var = IgM_ct) %>%
                    select(person_id, IgM_reconstitute_3mon = reconstitute_3mon, 
                            IgM_reconstitute_6mon = reconstitute_6mon, 
                            IgM_reconstitute_9mon = reconstitute_9mon, 
                            IgM_reconstitute_12mon = reconstitute_12mon, 
                            IgM_ct, IgM_measurement_date = measurement_date)

cohort_covars <- cohort_covars %>% left_join(CD3_reconstitute, by = "person_id") %>%
                        left_join(CD4_reconstitute, by = "person_id") %>%
                        left_join(CD8_reconstitute, by = "person_id") %>%
                        left_join(IgM_reconstitute, by = "person_id")

# resolve cases for patients from multiple subgroups
cohort_covars <- cohort_covars %>% left_join(results_tbl("multi_subgroup_resolved") %>% collect(), by = "person_id") %>%
                mutate(scd_type = if_else(!is.na(scd_type.y), scd_type.y, scd_type.x)) %>%
                select(-scd_type.x, -scd_type.y) 
                
# leukemia_dx_px <- find_conditions(load_codeset("leukemia_dx_px")) %>% collect()
# leukemia_dx_px %>% output_tbl(name = "leukemia_dx_px")
leukemia_dx_px <- results_tbl(name = "leukemia_dx_px") %>% collect() %>% 
                        group_by(person_id) %>%
                        summarise(leukemia_start_date = min(condition_start_date, na.rm = TRUE), .groups = "drop") %>%
                        mutate(leukemia_dx = TRUE) %>%
                        ungroup()

cohort_covars <- cohort_covars %>% left_join(leukemia_dx_px, by = "person_id") %>%
                mutate(leukemia_dx = if_else(is.na(leukemia_dx), FALSE, leukemia_dx)) %>%
                mutate(leukemia_prior = if_else(leukemia_start_date <= transplant_date, TRUE, FALSE))

                
# get conditioning agent
conditioning_agents <- cohort_covars %>% distinct(person_id, transplant_date) %>%
                        left_join(results_tbl("conditioning_agents_rx_px") %>% collect(), by = "person_id") %>% 
                        filter(exposure_date <= transplant_date, exposure_date >= transplant_date - days(30)) %>%
                        group_by(person_id, transplant_date, drug_type, conditioning_type) %>%
                        slice_max(exposure_date, n = 1, with_ties = FALSE) %>% ungroup() 

conditioning_agents_ct <-  conditioning_agents %>% group_by(person_id, transplant_date) %>% 
                                summarise(no_conditioning_agents = n_distinct(conditioning_type)) %>% 
                                ungroup()

conditioning_drugs_ct <-  conditioning_agents %>% filter(conditioning_type == "drugs") %>% 
                                mutate(has_busulfan = if_else(grepl("busulfan", drug_type, ignore.case = TRUE), 1, 0)) %>%
                                group_by(person_id, transplant_date) %>% 
                                summarise(no_conditioning_drugs = n_distinct(drug_type), 
                                        has_busulfan = sum(has_busulfan, na.rm = TRUE)) %>% 
                                ungroup() %>%
                                mutate(has_busulfan = if_else(has_busulfan >= 1, TRUE, FALSE))

conditioning_agents_tbl <- cohort_covars %>% distinct(person_id, transplant_date) %>% 
                        left_join(conditioning_agents_ct, by = c("person_id", "transplant_date")) %>%
                        left_join(conditioning_drugs_ct, by = c("person_id", "transplant_date")) %>% 
                        # distinct(no_conditioning_agents, no_conditioning_drugs) %>%
                        mutate(conditioning_type = case_when(is.na(no_conditioning_agents) ~ "None",
                                                            no_conditioning_agents == 1 & is.na(no_conditioning_drugs) ~ "tbi only",
                                                            no_conditioning_agents == 1 & no_conditioning_drugs >= 1 ~ "drugs only",
                                                            no_conditioning_agents == 2 ~ "drugs & tbi",
                                                            TRUE ~ NA))

cohort_covars <- cohort_covars %>%  
                        left_join(conditioning_agents_tbl, by = c("person_id", "transplant_date"))

# combine match_status, transplant_type, and donor_relation into a single category
cohort_covars <- cohort_covars %>% mutate(transplant_type_combined_old = case_when(transplant_type == "Autologous/Gene therapy" ~ "autologous/genetherapy",
                                                match_status %in% c("8/8", "10/10") & donor_relation == "Related" ~ "matched related",
                                                match_status %in% c("8/8", "10/10") & donor_relation == "Unrelated" ~ "matched unrelated",
                                                match_status %in% c("8/10", "9/10", "7/8") & donor_relation == "Unrelated" ~ "mismatched unrelated",
                                                match_status %in% c("8/10", "9/10", "7/8") & donor_relation == "Related" ~ "mismatched related",
                                                match_status %in% c("Haploidentical", "5/8") ~ "Haploidentical",
                                                match_status == "Cord Blood" ~ "Allogeneic: Cord Blood",
                                                TRUE ~ "Unknown")) 
# new criteria for transplant type is
cohort_covars <- cohort_covars %>% mutate(transplant_type_combined = case_when(transplant_type == "Autologous/Gene therapy" ~ "autologous/genetherapy",
                                                match_status %in% c("8/8", "10/10") & donor_relation == "Related" ~ "matched related",
                                                match_status == "Allogeneic: Cord Blood" ~ "matched unrelated/mismatched unrelated/mismatched related/cord blood", 
                                                !is.na(match_status) & !is.na(donor_relation) ~ "matched unrelated/mismatched unrelated/mismatched related/cord blood", 
                                                TRUE ~ "Unknown")) 

# cohort_covars %>% distinct(transplant_type_combined)

# get the list of patients with haploidentical unrelated status
# cohort_covars %>% filter(match_status == "Haploidentical" & donor_relation == "Unrelated") %>% select(record_id) 
# cohort_covars %>% filter(match_status %in% c("8/10", "9/10", "7/8"), donor_relation == "Related" ) %>% select(record_id)
# get the record_id for the haploidentical autologous patient
# cohort_covars %>% filter(transplant_type == "Autologous/Gene therapy" & match_status == "Haploidentical") %>% select(record_id)                                                
# get the record for unknown match_status and allogeneic
# cohort_covars %>% filter(transplant_type == "Allogeneic", is.na(match_status), eligibility == "Yes", chart_completion) %>% select(record_id)
# who had unknown status
cohort_covars %>% filter(transplant_type_combined == "Unknown") %>% select(record_id, match_status, donor_relation) %>% view()

# patients with prior leukemia diagnosis
# cohort_covars %>% filter(leukemia_prior, eligibility == "Yes") %>% select(record_id) %>% view()
cohort_covars %>% output_tbl("analytics_dataset")   

# number of transplant patients with no qualifying LIC or ferritin measuremnets in the first year
LIC <- results_tbl("cr_data") %>% select(record_id, site, eligibility, LIC_date, LIC) %>%
            left_join(results_tbl("study_cohorts") %>% select(person_id, record_id), by = "record_id") %>%
            filter(!is.na(LIC))
ferritin <- results_tbl("covar_ferritin_mx") %>% select(person_id, ferritin, ferritin_date)

LIC_patients <- results_tbl("no_multi_transplant_px") %>% left_join(LIC, by = "person_id") %>% 
            filter(transplant_date >= LIC_date - days(365),
                transplant_date <= LIC_date + days(365)) %>%
            distinct(person_id)

ferritin_patients <- results_tbl("no_multi_transplant_px") %>% left_join(ferritin, by = "person_id") %>% 
            filter(transplant_date >= ferritin_date - days(365),
                    transplant_date <= ferritin_date + days(365)) %>%
            distinct(person_id)

results_tbl("no_multi_transplant_px") %>% distinct(person_id) %>% 
                anti_join(LIC_patients %>% full_join(ferritin_patients)) %>% distinct_ct()

# get vod patients for Nora
results_tbl("analytics_dataset") %>% 
        filter(no_LIC_post>=1 | no_LIC_pre >=1 | no_ferritin_pre >= 1| no_ferritin_post >= 1) %>% 
        filter(vod) %>% filter(site == "chop") %>% 
        distinct(person_id) %>%
        left_join(results_tbl("study_cohorts") %>% select(person_id, record_id), by = "person_id") %>%
        output_tbl("vod_patients_chop", file = TRUE, local = TRUE)


# write out attrition data for aim 3
init_sum(cohort = 'Start', persons = 0)

# ferritin reductions for aim_3
# potential aim 3, n = 158
aim_3 <- results_tbl("study_cohorts") %>% filter(aim_3_2 | aim_3_1)
append_sum(cohort = 'No. patients identified by initial feasibility counts for aim 3',
             persons = distinct_ct(aim_3))

# exclude false positive patients from chart reviews, n = 139
aim_3 <- aim_3 %>% anti_join(results_tbl("cr_data") %>%
                    filter(eligibility == "No") %>% select(record_id), by = "record_id")
append_sum(cohort = 'Above counts excluding false positives from chart reviews',
             persons = distinct_ct(aim_3))
no_ferritin_aim3_ct <- aim_3 %>% distinct_ct()

# if the patient had chart reviews, use the transplant_date from chart reviews
aim_3 <- aim_3 %>% left_join(results_tbl("analytics_dataset") %>% filter(chart_completion, eligibility == "Yes") %>% 
                                  distinct(record_id, transplant_date) %>% mutate(is_cr = TRUE), by = "record_id") %>%
                    mutate(transplant_date.x = ifelse(is_cr, transplant_date.y, transplant_date.x)) %>%
                    rename(transplant_date = transplant_date.x) %>%
                    select(-transplant_date.y)

# get ferritin data, n = 139 patients with valid ferritin 
aim3 <- aim_3 %>%
      select(record_id, person_id, transplant_date) %>%
      left_join(results_tbl("covar_ferritin_mx"), by = "person_id") %>% 
      filter(transplant_date <= ferritin_date) %>% 
      filter(ferritin_date - transplant_date <= 600) %>%
      filter(!is.na(ferritin)) %>%
      group_by(person_id) %>% 
      mutate(n = n()) %>% filter(n >=1) %>% ungroup()
      
aim3_ct <- aim_3 %>% distinct_ct()
append_sum(cohort = 'Above counts with at least 1 ferritin measurement within 1 year of transplant',
             persons = distinct_ct(aim_3))

# find therapeutic phlebotomy, n = 43
phleb_px <- aim_3 %>%
      # filter(chart_completion, eligibility == "Yes") %>%
      select(record_id, person_id, transplant_date) %>%
      find_procedures(procedure_codeset_name = "therapeutic_phleb_px") %>% 
      filter(procedure_date >= transplant_date)

append_sum(cohort = 'No. patients with therapeutic phlebotomy after transplant',
             persons = distinct_ct(phleb_px))

phleb_px_1yr <- phleb_px %>%      
      filter(procedure_date - transplant_date <= 365) %>% 
      group_by(person_id) %>%
      collect_new() %>%
      mutate(start_date = min(procedure_date, na.rm = TRUE), 
             end_date = max(procedure_date, na.rm = TRUE),
             duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
      distinct(person_id, record_id, start_date, end_date, transplant_date, duration) %>% 
      mutate(route = NA) %>%
      mutate(IRT_type = "phlebotomy") 

append_sum(cohort = 'No. patients with therapeutic phlebotomy within 1 year of transplant',
             persons = distinct_ct(phleb_px_1yr))

# patients who got chelation, n = 20
dx_codeset<- load_codeset("deferoxamine_rx") %>% mutate(type = "deferoxamine") %>% #deferoxamine
                union(load_codeset("deferasirox_rx") %>% mutate(type = "deferasirox")) %>% #deferasirox
                compute_new(temp = TRUE, name = "drug_id")

# patients who got chelation, n = 20
chelation_rx <- aim_3 %>%
      # filter(chart_completion, eligibility == "Yes") %>%
      select(record_id, person_id, transplant_date) %>%
      find_drugs(dx_codeset) %>% 
      filter(drug_exposure_start_date >= transplant_date) %>% collect() 

append_sum(cohort = 'No. patients received deferoxamine/deferasirox after transplant',
             persons = distinct_ct(chelation_rx))

chelation_rx <- chelation_rx %>%
      filter(drug_exposure_start_date - transplant_date <= 365) %>% 
      distinct(person_id, drug_exposure_start_date, IRT_type = type, .keep_all = TRUE) 

append_sum(cohort = 'No. patients received deferoxamine/deferasirox within 1 year after transplant',
             persons = distinct_ct(chelation_rx))

# for intravenous drugs, assume a single one per day, n = 3
chelation_iv <- chelation_rx %>% filter(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous")) %>%
      group_by(person_id, IRT_type) %>%
      arrange(person_id, drug_exposure_start_date) %>% 
      # collect_new() %>% 
      mutate(start_date = min(drug_exposure_start_date, na.rm = TRUE), 
             end_date = max(drug_exposure_start_date, na.rm = TRUE),
             duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
      distinct(person_id, record_id, start_date, end_date, transplant_date, duration, .keep_all = TRUE) %>% 
      mutate(route = "IV") %>%
      select(person_id, record_id, start_date, end_date, transplant_date, duration, IRT_type, route)

# oral drugs, n = 14
chelation_oral <- chelation_rx %>% filter(!(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous")) | 
                                          is.na(route_source_value)) %>%
      filter(!is.na(quantity) | !is.na(days_supply)) %>% 
      # collect_new() %>%
      mutate(refills = if_else(is.na(refills), 0, refills)) %>%
      arrange(person_id, drug_exposure_start_date) %>%
      # select(person_id, drug_exposure_start_date, quantity, refills, IRT_type, route_source_value) %>%
      group_by(person_id, IRT_type) %>% 
      mutate(start_date = drug_exposure_start_date, # min(drug_exposure_start_date, na.rm = TRUE), 
             end_date = drug_exposure_start_date + days(quantity*(refills+1)),
             duration = quantity) %>% ungroup() %>%
      distinct(person_id, record_id, start_date, end_date, transplant_date, duration, quantity, refills, .keep_all = TRUE) %>% 
      select(person_id, record_id, start_date, end_date, transplant_date, duration, quantity, refills, IRT_type, route = route_source_value) %>% 
      arrange(person_id) 

# Only took the duration of the first prescription, some prescriptions have up to 11 refills 
# patients with missing both days supply and quantity, n = 4
missing_duration <- chelation_rx %>% anti_join(chelation_oral, by = "person_id") %>% anti_join(chelation_iv, by = "person_id") %>% 
      select(person_id, record_id, drug_exposure_start_date, transplant_date, days_supply, quantity, refills, IRT_type, route = route_source_value) 

# append the chelation and IRT patients
irt <- phleb_px_1yr %>% collect() %>%
          union(chelation_iv) %>% mutate(quantity = NA, refills = NA) %>% 
          union(chelation_oral) 

append_sum(cohort = 'No. patients received deferoxamine/deferasirox within 1 year of transplant',
             persons = distinct_ct(irt %>% filter(IRT_type != "phlebotomy")))
append_sum(cohort = 'No. patients received deferoxamine/deferasirox within 1 year of transplant with missing duration',
             persons = distinct_ct(missing_duration))

irt_exclude <- irt %>% group_by(record_id) %>% 
      summarise(n = n()) %>% 
      filter(n> 1) %>% ungroup()

append_sum(cohort = 'No. patients received both phlebotomy and deferoxamine/deferasirox with non-overlapping period (excluded)',
             persons = irt_exclude %>% distinct_ct("record_id"))

irt_exclude %>% 
      inner_join(irt, by = "record_id") %>% 
      arrange(record_id) %>%
      select(-person_id) %>% view()

# there are patients who got more than 1 IRT treaments, we do have to double check
irt <- irt %>% group_by(record_id, IRT_type) %>%
          mutate(start_date = min(start_date, na.rm = TRUE)) %>%
          mutate(end_date = max(end_date, na.rm = TRUE))

irt <- irt %>% anti_join(irt_exclude, by = "record_id")

ferritin_start <- irt %>% copy_to_new(df = ., name = "irt") %>%
      inner_join(results_tbl("covar_ferritin_mx") %>%
                  select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
      filter(transplant_date <= ferritin_date) %>%
      filter(!is.na(ferritin)) %>%
      # filter(transplant_date - ferritin_date <= 365) %>%
      group_by(person_id) %>% 
      filter(start_date <= ferritin_date) %>%
      slice_max(start_date, n =1, with_ties = FALSE) %>% ungroup() %>%
      select(person_id, record_id, transplant_date, start_date, ferritin0_date = ferritin_date, ferritin0 = ferritin) 
      

ferritin_end <- irt %>% copy_to_new(df = ., name = "irt") %>%
      inner_join(results_tbl("covar_ferritin_mx") %>%
                  select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
      filter(transplant_date <= ferritin_date) %>%
      filter(!is.na(ferritin)) %>%
      filter(end_date <= ferritin_date) %>%
      group_by(person_id) %>% 
      slice_min(end_date, n =1, with_ties = FALSE) %>% ungroup() %>%
      select(person_id, record_id, transplant_date, end_date, ferritinf_date = ferritin_date, ferritinf = ferritin)

ferritin_all <- irt %>% copy_to_new(df = ., name = "irt") %>% 
      inner_join(results_tbl("covar_ferritin_mx") %>%
                  select(person_id, ferritin_date, ferritin), by = "person_id") %>%
      filter(transplant_date <= ferritin_date) %>%
      filter(!is.na(ferritin)) %>%
      filter(transplant_date - ferritin_date <= 365) %>% collect_new()

ferritin_0 <- ferritin_all %>% group_by(person_id, IRT_type) %>%
          slice_min(abs(start_date - ferritin_date), n = 1, with_ties = FALSE) %>% ungroup() %>%
          select(record_id, person_id, IRT_type, start_date, ferritin_date_0 = ferritin_date, ferritin_0 = ferritin)

ferritin_f <- ferritin_all %>% group_by(person_id, IRT_type) %>%
          slice_min(abs(end_date - ferritin_date), n = 1, with_ties = FALSE) %>% ungroup() %>%
          select(record_id, person_id, IRT_type, start_date, ferritin_date_f = ferritin_date, ferritin_f = ferritin)
ferritin_change <- ferritin_0 %>% left_join(ferritin_f, by = c("record_id", "person_id", "IRT_type", "start_date")) %>%
          mutate(ferritin_change_rate = (ferritin_f - ferritin_0)/as.numeric((ferritin_date_f - ferritin_date_0))) 

table1(~ferritin_change_rate| IRT_type, 
      data = ferritin_change,
      caption = "Post-transplant iron reduction therapy",
      overall = FALSE)

# 47 patients, 
# only 37 patients with duration > 0 (the other 10 patients only had a single phlebotomy procedure) 
# remove 3 more for missing ferritin values

p1 <- ferritin_all %>% 
      mutate(ferritin_date = as.numeric(ferritin_date - transplant_date)) %>%
      mutate(start_date = as.numeric(start_date - transplant_date)) %>%
      mutate(end_date = as.numeric(end_date - transplant_date)) %>% 
      filter(duration > 0) %>% 
      ggplot(aes(x = ferritin_date, y = ferritin, group = record_id)) + 
      geom_point() + 
      geom_line() +
      geom_vline(aes(xintercept = start_date), linetype = "dashed", color = "red") +
      geom_vline(aes(xintercept = end_date), linetype = "dashed", color = "red") +
      xlim(c(0, 500)) + 
      facet_wrap(~person_id, ncol = 3)

ggsave("reporting/irt_ferritin_trjectories.png", p1, width = 20, height = 30, units = "cm", dpi = 300)

# now slice the ferritin values
# ferritin_0 is the closet value before the start date within 30 days
# 23 patients 
# test has 34 patients after removing 3 patients with missing ferritin values 
# limit values to just 1 year, 29 patients with valid ferritin values 
append_sum(cohort = 'No. patients with only 1 ferritin measuremnt (excluded)',
             persons = 3)

test <- ferritin_all %>% 
      mutate(ferritin_date = as.numeric(ferritin_date - transplant_date)) %>%
      mutate(start_date = as.numeric(start_date - transplant_date)) %>%
      mutate(end_date = as.numeric(end_date - transplant_date)) %>% 
      filter(!(person_id %in% c(4502732, 599793, 9720611))) %>%  #remove patients with missing ferritin values
      filter(duration > 0) %>%
      filter(ferritin_date <= 365) %>% 
      mutate(end_date = case_when(person_id == 298603 ~ start_date + quantity, 
                                TRUE ~ end_date)) # fix the durations for this patients
      # filter(ferritin_date <= start_date) %>%
      # filter(ferritin_date >= start_date - 30) %>%
      # group_by(person_id) %>%
      # slice_max(ferritin_date, n = 1, with_ties = FALSE) %>% ungroup() 

append_sum(cohort = 'No. patients with at least 2 ferritin measurements after transplant and treatment durations > 0',
             persons = distinct_ct(test))

ferritin_start <- linear_interpolate_per_patient(df = test, df_interpolated = test %>% distinct(person_id, start_date) %>% rename(c("interpolated_t" = "start_date")))
ferritin_start %>% view()

ferritin_end <- linear_interpolate_per_patient(df = test, df_interpolated = test %>% distinct(person_id, end_date) %>% rename(c("interpolated_t" = "end_date")))
ferritin_end %>% view()

# fit a linear line to the ferritin change
final_ferritin_trjectories <- ferritin_start %>% select(person_id, t0 = interpolated_t, f0 = interpolated_ferritin) %>%
      full_join(ferritin_end %>% select(person_id, t1 = interpolated_t, f1 = interpolated_ferritin), by = "person_id") %>%
      filter(!is.na(f0), !is.na(f1)) %>%
      mutate(f_change = (f0-f1)/(t1-t0)/30) 

append_sum(cohort = 'No. patients with at least 2 ferritin measurements, each is within 60 days of treatment start and end dates',
             persons = distinct_ct(final_ferritin_trjectories))

# output attrition
output_sum(name = "attrition_aim3", local = TRUE, file = TRUE)

final_ferritin_trjectories %>% ggplot(aes(x = t1-t0, y = f_change)) + geom_point() + geom_smooth() 
      
ferritin_all %>% output_tbl("ferritin_trajectories_IRT")


# moves table from db 5.1 to db 5.5
select(person_id, record_id, transplant_date)
cohort <- results_tbl("study_cohorts") %>% collect()
cohort_covars <- results_tbl("analytics_dataset") %>% collect() 
cr_data <- results_tbl("cr_data") %>% collect() 

cohort %>% output_tbl("study_cohorts")
cohort_covars %>% output_tbl("analytics_dataset")
cr_data %>% output_tbl("cr_data")
