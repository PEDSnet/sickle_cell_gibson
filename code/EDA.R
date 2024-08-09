
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

# extract most recent LIC values                                 
rslt$cr_data_LIC_recent <- results_tbl("cr_data") %>% collect() %>% 
                                val_extraction(cohort = ., 
                                                no_value = no_LIC, grouping_id = record_id, 
                                                value_name = LIC, value_date = LIC_date, 
                                                cutoffs = LIC_cutoffs,
                                                value_type = "LIC", keep_tie = keep_ties_LIC, 
                                                slice_by = "most_recent")   
# extract highest LIC values 
rslt$cr_data_LIC_max <- results_tbl("cr_data") %>% collect() %>% 
                                val_extraction(cohort = ., 
                                                no_value = no_LIC, grouping_id = record_id, 
                                                value_name = LIC, value_date = LIC_date, 
                                                cutoffs = LIC_cutoffs,
                                                value_type = "LIC", keep_tie = keep_ties_LIC, 
                                                slice_by = "max") 

# merge the most recent and highest LIC values 
rslt$cr_data <- rslt$cr_data_LIC_recent %>% 
                        left_join(rslt$cr_data_LIC_max %>% select(record_id, LIC_max, LIC_3_max, 
                                                                LIC_days_since_transplant_max, LIC_days_since_transplant_3_max,
                                                                LIC_type_max,
                                                                LIC_max_level, LIC_3_max_level), by = c("record_id", "LIC_type" = "LIC_type_max"))                              

rslt$cr_data <- results_tbl("cr_data") %>% select(-c("LIC_other_unit", "LIC_date", 
                                                "LIC_raw", "LIC_est_method", "is_ms", "LIC", "LIC_type", "site", "no_LIC_measurement")) %>% 
                                        collect() %>% distinct() %>%
                        left_join(rslt$cr_data %>% select("record_id", "LIC_type", "LIC_days_since_transplant",
                                                        "LIC_3", "LIC_days_since_transplant_3", "LIC_level", "LIC_3_level", 
                                                        "LIC_max", "LIC_3_max", "LIC_days_since_transplant_max", "LIC_days_since_transplant_3_max",
                                                        "LIC_max_level", "LIC_3_max_level"), by = "record_id") %>%                               
                                rename(transplant_date_cr = transplant_date, 
                                transplant_type_cr = transplant_type, 
                                second_transplant_date_cr = second_transplant_date,
                                disease_relapse = graft_fail) %>% # Disease relapse: 
                                mutate(LIC_level_3_binary = if_else(LIC_3 < LIC_cutoffs_binary[1], "low", "high"),
                                LIC_level_3_max_binary = if_else(LIC_3_max < LIC_cutoffs_binary[1], "low", "high"),
                                across(c("LIC_level_3_max_binary", "LIC_level_3_binary"), ~factor(.x, levels = c("low", "high")))) %>%
                                copy_to_new(df = ., name = "cr_data") 
# rslt$cr_data %>% view()                                
                                
rslt$study_cohorts <- results_tbl("study_cohorts") %>% 
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        left_join(rslt$cr_data, by = "record_id") %>% 
                        mutate(transplant_date_consistency = if_else(transplant_date_cr == transplant_date, 1, 0),
                                transplant_date = if_else(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                                transplant_type = if_else(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                                chart_completion = if_else(is.na(chart_completion), FALSE, chart_completion)) %>%
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

ferritin_post <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "post", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)

# slice by highest values 
ferritin_pre_max <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "pre_max", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)

ferritin_post_max <- rslt$study_cohorts %>% select(person_id, transplant_date) %>% 
                left_join(rslt$ferritin, by = "person_id") %>%
                ferritin_classification(ferritin_type = "post_max", no_ferritin = no_ferritin, keep_tie = keep_tie,
                                        ferritin_cutoff1 = ferritin_cutoff1, ferritin_cutoff2 = ferritin_cutoff2)


ferritin <- ferritin_pre %>% select(person_id, ferritin_pre = ferritin_level, ferritin_days_pre = ferritin_days_since_transplant) %>%
                full_join(ferritin_post %>% select(person_id, ferritin_post = ferritin_level, ferritin_days_post = ferritin_days_since_transplant), by = "person_id") %>%
                full_join(ferritin_pre_max %>% select(person_id, ferritin_pre_max = ferritin_level, ferritin_days_pre_max = ferritin_days_since_transplant), by = "person_id") %>%
                full_join(ferritin_post_max %>% select(person_id, ferritin_post_max = ferritin_level, ferritin_days_post_max = ferritin_days_since_transplant), by = "person_id") %>%
                mutate(across(c("ferritin_post", "ferritin_pre", "ferritin_post_max", "ferritin_pre_max"), ~factor(.x, levels = c("low", "moderate", "high"))))

cohort_covars <- rslt$study_cohorts %>% left_join(ferritin %>% copy_to_new(df = ., name = "sdsd"), by = c("person_id"))

# cohort_covars %>% view()
# patients that belong to aim_2a_2 with missing pre-ferritin should be excluded from that aim
# no patients not belong to any aims 
# rm("ferritin_pre", "ferritin_post", "ferritin_pre_max", "ferritin_post_max")
cohort_covars <- cohort_covars %>% collect() %>% 
                                mutate(aim_2a_2 = if_else(aim_2a_2 & is.na(ferritin_pre), FALSE, aim_2a_2)) %>%
                                mutate(aim_2b_2 = if_else(aim_2b_2 & is.na(ferritin_post), FALSE, aim_2b_2)) %>%
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
vod <- results_tbl("covar_vod_dx") %>% 
        left_join(rslt$study_cohorts %>% select(person_id, transplant_date), by = "person_id") %>%
        mutate(vod_date := pmin(min_defi_abd_date, min_bili_abd_date, min_bili_abd_weight_date, na.rm = TRUE)) %>%
        collect() %>% filter(as.numeric(difftime(vod_date, transplant_date, units = "days")) <= 90 |
                                (vod_date >= transplant_date) | 
                                abs(as.numeric(difftime(abd_date, transplant_date, units = "days"))) <= 7) %>% 
        distinct(person_id, vod_cat_1, vod_cat_2, vod_cat_3, vod_cat_4, vod_cat_5, vod) 
vod_cat2_ids <- vod %>% filter(vod_cat_2 == 2, is.na(vod_cat_1), is.na(vod_cat_3), is.na(vod_cat_4), is.na(vod_cat_5))
vod <- vod %>% anti_join(vod_cat2_ids, by = "person_id")
vod %>% distinct_ct()
# vod %>% arrange(person_id) %>% view()        
cohort_covars <- cohort_covars %>% left_join(vod, by = "person_id") %>%
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

# patients from multiple subgroups

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
                        slice_max(exposure_date, n = 1, with_ties = FALSE) %>% ungroup() %>% view()

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
cohort_covars <- cohort_covars %>% mutate(transplant_type_combined = case_when(transplant_type == "Autologous/Gene therapy" ~ "Autologous/Gene therapy",
                                                match_status == "10/10" & donor_relation == "Related" ~ "Allogeneic: 10/10 Related",
                                                match_status == "10/10" & donor_relation == "Unrelated" ~ "Allogeneic: 10/10 Unrelated",
                                                match_status %in% c("8/10", "9/10") & donor_relation == "Related" ~ "Allogeneic: 8/10 or 9/10",
                                                match_status == "Cord Blood" ~ "Allogeneic: Cord Blood",
                                                match_status == "Haploidentical" ~ "Allogeneic: Haploidentical",
                                                TRUE ~ "other")) 

# get the record_id for the haploidentical autologous patient
cohort_covars %>% filter(transplant_type == "Autologous/Gene therapy" & match_status == "Haploidentical") %>% select(record_id)                                                
# get the record for unknown match_status and allogeneic
cohort_covars %>% filter(transplant_type == "Allogeneic", is.na(match_status), eligibility == "Yes", chart_completion) %>% select(record_id)

# patients with prior leukemia diagnosis
cohort_covars %>% filter(leukemia_prior, eligibility == "Yes") %>% select(record_id) %>% view()
# cohort_covars %>% output_tbl("analytics_dataset")                        