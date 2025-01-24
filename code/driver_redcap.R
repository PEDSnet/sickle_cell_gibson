.run  <- function() {

    library(tidyr)
    library(ggplot2)
    
    setup_pkgs()

    # Stanford always have some weird data
    # 08_0017 has an LIC date but no LIC_raw value
    # 08_0019_cr measurement on 2019-06-06 R2* require no conversion
    test <- read.csv("redcap/IronOverloadChartRev_DATA_2024-11-01_1059.csv")
    test %>% filter(redcap_data_access_group == "cchmc") %>% select(record_id, starts_with("iron_val_")) %>% view()

    # get the most recent redcap data, filter out old patient ids, extract completed records
    cr_data <- get_redcap_data(redcap_filename = "redcap/IronOverloadChartRev_DATA_2024-11-26_2219.csv")
    cr_data %>% view()
    
    # number of completed chart reviews n = 600 
    num_cr_completed <- cr_data %>% distinct_ct("record_id")
    num_false_positives <- cr_data %>% filter(eligibility == "No") %>% distinct_ct("record_id") # 102
    num_tru_positives <- cr_data %>% filter(eligibility == "Yes") %>% distinct_ct("record_id") # 498
                            
    # is there any record without graft manipulation status, the best fix is to manually fix their redcap form 
    cr_data %>% filter(is.na(graft_manip), eligibility == "Yes") %>% view()

    # check for consistency between graft failure and graft failure dates:
    # all graft failures should have a date
    cr_data %>% filter(graft_fail == "Yes", is.na(graft_fail_date)) %>% distinct(record_id, graft_fail, graft_fail_date) %>% view()

    # is there any record without graft failure status, the best fix is to manually fix their redcap form
    cr_data %>% filter(is.na(graft_fail), eligibility == "Yes") %>% distinct(record_id, graft_fail, graft_fail_date, eligibility) %>% view()

    # filter out false postives and patients without LIC values
    # LIC unit conversion: Note all LIC values are expected in mg Fe / g dry weight
    # R2* from chop and colorado: no conversion needed, other institutions: ask Nora
    cr_data %>% filter(eligibility == "Yes", !is.na(LIC_date)) %>% 
                group_by(site, LIC_other_unit, LIC_est_method) %>% 
                summarise(n = n_distinct(record_id)) %>%  
                ungroup() %>% filter(LIC_est_method == "R2*") %>% #is.na(LIC_other_unit) | is.na(LIC_est_method)
                view()

    # T2* from stanford: 
    # T2* from colorado: no conversion needed
    # other institutions: ask Nora
    cr_data %>% filter(eligibility == "Yes", !is.na(LIC_date)) %>% 
                group_by(site, LIC_other_unit, LIC_est_method) %>% 
                summarise(n = n_distinct(record_id)) %>%  
                ungroup() %>% filter(LIC_est_method == "T2*") %>% #is.na(LIC_other_unit) | is.na(LIC_est_method)
                view()

    # to check for a specific value
    cr_data %>% filter(eligibility == "Yes", !is.na(LIC_date), LIC_est_method == "T2*", site == "stanford", is.na(LIC_other_unit)) %>% view()
    
    cr_data <- cr_data %>% mutate(is_ms = case_when(grepl("milliseconds", LIC_other_unit, ignore.case = TRUE) ~ TRUE, 
                                        grepl("ms", LIC_other_unit, ignore.case = TRUE) ~ TRUE, 
                                        is.na(LIC_other_unit) ~ NA,
                                        TRUE ~ FALSE),
                        is_Hz = case_when(record_id == "04_0001_cr" & LIC_date == as.Date("2007-07-19") ~ TRUE,
                                          record_id == "04_0001_cr" & LIC_date == as.Date("2017-03-23") & LIC_est_method == "R2" ~ TRUE, 
                                          record_id == "04_0002_cr" & LIC_date == as.Date("2007-09-04") ~ TRUE, 
                                          record_id == "04_0002_cr" & LIC_date == as.Date("2011-12-13") & LIC_est_method == "R2" ~ TRUE, 
                                          record_id == "04_0003_cr" & LIC_date == as.Date("2009-08-07") ~ TRUE, 
                                          record_id == "04_0004_cr" & LIC_date == as.Date("2007-08-21") ~ TRUE, 
                                          record_id == "04_0005_cr" & LIC_date == as.Date("2007-08-07") ~ TRUE, 
                                          record_id == "04_0008_cr" & LIC_date == as.Date("2009-09-02") ~ TRUE, 
                                          record_id == "04_0009_cr" & LIC_date == as.Date("2007-07-24") ~ TRUE, 
                                          record_id == "04_0013_cr" & LIC_date == as.Date("2009-07-22") & LIC_est_method == "R2*"~ TRUE, 
                                          record_id == "04_0014_cr" & LIC_date == as.Date("2009-08-07") & LIC_est_method == "R2*"~ TRUE, 
                                          record_id == "04_0021_cr" & LIC_date == as.Date("2019-03-31") & LIC_est_method == "R2"~ TRUE, 
                                          record_id == "04_0024_cr" & LIC_date == as.Date("2009-09-04") ~ TRUE, 
                                          record_id == "04_0044_cr" & LIC_date == as.Date("2016-11-29.0254*33") & LIC_est_method == "R2"~ TRUE, 
                                          TRUE ~ FALSE),                                
                        LIC = case_when(record_id == "08_0019_cr" & LIC_est_method == "R2*" ~ LIC_raw,
                                        site %in% c("chop", "colorado") ~ LIC_raw, # no conversion needed for these sites 
                                        !is_ms | is.na(is_ms) ~ LIC_raw, 
                                        is_ms ~ 1000/LIC_raw*0.0254 + 0.202, 
                                        TRUE ~ NA),
                        LIC = ifelse(is_Hz, LIC_raw*0.0254 + 0.202, LIC), 
                        LIC_type = case_when(LIC_date >= transplant_date ~ "post-transplant", 
                                             LIC_date < transplant_date ~ "pre-transplant", 
                                             TRUE ~ NA)) 
    # double check
    cr_data %>% filter(eligibility == "Yes", !is.na(LIC_date)) %>% 
                select(site, record_id, LIC_date, LIC_est_method, LIC_other_unit, is_ms, LIC_raw, LIC) %>% 
                arrange(site, record_id, LIC_date, LIC_est_method, is_ms) %>% 
                view()
                # group_by(LIC_other_unit, LIC_est_method, is_ms) %>% summarise(n = n_distinct(record_id)) %>% view()
    cr_data %>% filter(eligibility == "Yes", site == "lurie", LIC_est_method != "R2", !(LIC_other_unit %in% c("NA", "N/A", "na", "n/a"))) %>% distinct(record_id) %>% view()

    # check for lurie value
    lurie <- cr_data %>% filter(eligibility == "Yes") %>%
                filter(site == "lurie") %>%
                filter(!is.na(LIC_date)) 
                # select(record_id, LIC_date, LIC_raw, LIC_est_method, is_ms, is_Hz, LIC, LIC_type, LIC_other_unit) 
    
    # output Lurie data for Nora
    lurie %>% select(record_id, LIC_date, LIC_raw, LIC_est_method, is_ms, is_Hz, LIC, LIC_type, LIC_other_unit)  %>%
                arrange(record_id, LIC_date) %>% output_tbl(name = "lurie_LIC_data", file = TRUE, local = TRUE)

    # for dates with multiple LIC values, we will use the R2 values, aka exclude R2* values
    # for dates with multiple R2 values, we will use average values            
    lurie_multi_per_day <- lurie %>% group_by(record_id, LIC_date) %>% 
                summarise(n = n()) %>% ungroup() %>% filter(n >= 2) %>% 
                left_join(lurie, by = c("record_id", "LIC_date")) %>% 
                filter(LIC_est_method == "R2*") #%>% 
                # group_by(record_id, LIC_date, LIC_est_method ) %>% 
                # summarise(LIC = mean(LIC, na.rm = TRUE)) %>% ungroup() 
    lurie <- lurie %>% anti_join(lurie_multi_per_day, by = c("record_id", "LIC_date", "LIC_est_method"))        

     cr_data <- cr_data %>% anti_join(lurie, by = c("record_id", "LIC_date")) %>%
                union(lurie)

    # write cr data to a table
    # patients with match status that were not included in the study

    cr_data <- cr_data %>% mutate(match_status = case_when(record_id %in% c("01_0094_cr", "01_0046_cr", "01_0003_cr", "01_0138_cr", "01_0110_cr", "01_0124_cr", "01_0008_cr", "01_0034_cr") ~ "8/8", 
                                                           record_id %in% c("01_0147_cr", "01_0050_cr", "01_0023_cr") ~ "7/8", 
                                                           record_id == "01_0188_cr" ~ "5/8", 
                                                           record_id == "01_0106_cr" ~ "8/10", 
                                                           TRUE ~ match_status))                                                             
    cr_data %>% output_tbl(name = "cr_data")

    # completion by site
    cr_data %>% filter(iron_overload_chart_review_complete == 2) %>% 
        mutate(eligibility = case_when(eligibility == 1 ~ "Eligible", 
                                        eligibility == 0 ~ "Not eligible", 
                                        TRUE ~ "Missing")) %>%
        group_by(site, eligibility) %>% 
        summarise(num_completed = n()) %>% 
        ungroup() %>%
        pivot_wider(names_from = eligibility, values_from = num_completed)  

    # let's check who did not have a transplant
    # aka false positives
    cr_data %>% filter(is.na(transplant_date)) %>% 
        select(record_id) %>%
        inner_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>% collect(), by = c("record_id")) %>% 
        inner_join(results_tbl("study_cohorts") %>% collect(), by = c("person_id", "record_id", "site")) %>% 
        select(record_id, person_id, transplant_date, transplant_type, aim_2a_2:aim_3_2) %>%
        output_tbl(name = "false_positive_ids")
    
    # let's see how well the transplant date and status match
    cr_concord <- cr_data %>% filter(!is.na(transplant_date)) %>% 
                        select(record_id, transplant_date_cr = transplant_date, 
                        transplant_type_cr = transplant_type, 
                        second_transplant_date_cr = second_transplant_date) %>% distinct() %>%
                        inner_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>%
                                    select(person_id, record_id) %>% collect(), by = c("record_id")) %>% 
                        left_join(results_tbl("study_cohorts") %>% collect() %>%
                                mutate(transplant_type = case_when(grepl("Allogeneic", transplant_concept_name, ignore.case = TRUE) ~ "Allogeneic",
                                                                    grepl("Autologous", transplant_concept_name, ignore.case = TRUE) ~ "Autologous/Gene therapy",
                                                                    TRUE ~ "other")) %>%
                                select(person_id, transplant_date, transplant_type), by = c("person_id"))
    cr_concord %>% view()

    # data quality check
    # patients needing second transplant from EHR
    transplant_px <- results_tbl("no_multi_transplant_px") %>% collect_new() 
    transplant_px <- transplant_px %>% group_by(person_id) %>%
                    mutate(transplant_days = as.numeric(difftime(transplant_date, min(transplant_date), units = "days"))) %>%
                    summarise(transplant_1 = min(transplant_days),
                            transplant_2 = max(transplant_days),
                            second_transplant_date = max(transplant_date),
                            transplant_date = min(transplant_date)) %>%
                    mutate(transplant_2 = if_else(transplant_2 == transplant_1, NA, transplant_2)) %>%
                    select(person_id, transplant_date, transplant_days = transplant_1, second_transplant_days = transplant_2, second_transplant_date) %>% 
                    ungroup() %>% filter(!is.na(second_transplant_days))
    
    # check for mismatches, only 3 out of 371
    # the case without transplant date or type are the extra cases that we could not verify
    cr_concord %>% left_join(transplant_px, by = c("person_id", "transplant_date")) %>%
        filter(transplant_date_cr != transplant_date | 
                transplant_type_cr != transplant_type | 
                second_transplant_date_cr != second_transplant_date) %>% view()    

    # check for exact matches  
    cr_concord %>% inner_join(transplant_px, by = c("person_id", "transplant_date")) %>%
        filter(transplant_date_cr == transplant_date | is.na, 
                transplant_type_cr == transplant_type, 
                second_transplant_date_cr == second_transplant_date | is.na(second_transplant_date_cr) | is.na(transplant_date_cr)) %>% view()  

    # get the new list of patient ids
    study_cohorts <- results_tbl("study_cohorts") %>% select(person_id, aim_2a_1:aim_3_2) %>%
                    union(results_tbl("multi_transplant_cohort") %>% anti_join(results_tbl("study_cohorts"), by = "person_id") %>% 
                            select(person_id, aim_2a_1:aim_3_2) ) %>%
                    collect()

    results_tbl("study_cohorts") %>% collect() %>% # no multiple transplants
                filter(!(aim_3_1| aim_3_2)) %>% 
                distinct(person_id, site) %>%
                union(results_tbl("multi_transplant_cohort") %>% collect() %>% #patients with unverifiable transplant status
                        # filter(site == sites[s]) %>% #, (aim_2a_2| aim_2b_2| aim_3_1| aim_3_2)) %>% 
                        distinct(person_id, site)) %>%
                distinct(person_id, site)

                
    study_cohorts <- cr_data %>% filter(!is.na(transplant_date)) %>% 
                        select(record_id, transplant_date, transplant_type, second_transplant_date) %>% distinct() %>%
                        inner_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>%
                                    select(person_id, record_id) %>% collect(), by = c("record_id")) %>% 
                        left_join(study_cohorts %>% select(person_id, aim_2a_1:aim_3_2), by = c("person_id")) 
                        
    study_cohorts %>% output_tbl(name = "cr_cohorts")
}
