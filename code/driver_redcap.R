.run  <- function() {

    setup_pkgs()

    cr_data <- read.csv("redcap/IronOverloadChartRev_DATA_2024-05-02_1145.csv")
    # filter out old patient ids
    cr_data <- cr_data %>% filter(!grepl("^[A-Za-z]", record_id))
    # extract completed records 
    cr_data <- cr_data %>% rename(site = redcap_data_access_group, 
                        chart_completion = iron_overload_chart_review_complete,
                        transplant_date = date_tpx,
                        transplant_type = tpx_type) %>%
                mutate(transplant_type = case_when(transplant_type == 1 ~ "Allogeneic", 
                                            transplant_type == 2 ~ "Autologous/Gene therapy", 
                                            TRUE ~ NA),
                        chart_completion = if_else(chart_completion == 2, TRUE, FALSE), 
                        across(iron_date_1:est_method_15, ~if_else(.x == "", NA, .x)),
                        across(c("transplant_date", "second_transplant_date"), ~if_else(.x == "", NA, .x)),
                        across(iron_val_1:iron_val_15, ~as.numeric(.x)),
                        # transplant_date = if_else(transplant_date == "", NA, as.Date(transplant_date)),
                        eligibility = case_when(is.na(transplant_date ) & is.na(eligibility) ~ 0,
                                                !is.na(transplant_date ) & is.na(eligibility) ~ 1,  
                                                TRUE ~ eligibility)) %>%
                filter(chart_completion) 

    # any eligible cases without any LIC measurements would be filtered out
    cr_data_eligible <- cr_data %>% pivot_longer(cols = c(starts_with("iron_date_")), names_to = "no_LIC_measurement", 
                            values_to = "LIC_date", values_drop_na = TRUE) %>% 
                mutate(no_LIC_measurement = as.numeric(sub(".*_([0-9]+)$", "\\1", no_LIC_measurement)),
                    LIC_date = as.Date(LIC_date)) %>% 
                pivot_longer(cols = starts_with("iron_val_"), names_to = "no_LIC_measurement1", 
                            values_to = "LIC", values_drop_na = TRUE) %>%
                mutate(no_LIC_measurement1 = as.numeric(sub(".*_([0-9]+)$", "\\1", no_LIC_measurement1))) %>% 
                filter(no_LIC_measurement == no_LIC_measurement1) %>%
                pivot_longer(cols = starts_with("est_method_"), names_to = "no_LIC_measurement2", 
                            values_to = "LIC_est_method", values_drop_na = TRUE) %>%
                mutate(no_LIC_measurement2 = as.numeric(sub(".*_([0-9]+)$", "\\1", no_LIC_measurement2))) %>% 
                filter(no_LIC_measurement == no_LIC_measurement2) %>% 
                    # ((!is.na(LIC_date) | !is.na(LIC) | !is.na(LIC_est_method))& eligibility != 0) | eligibility == 0) %>%
                select(-no_LIC_measurement1, -no_LIC_measurement2)

    cr_data_ineligible <- cr_data %>% select(-c(iron_date_1:est_method_15)) %>% 
                            anti_join(cr_data_eligible %>% distinct(record_id), by = "record_id") %>%
                            # filter(eligibility == 0) %>%
                            mutate(LIC = NA, LIC_date = NA, LIC_est_method = NA, no_LIC_measurement = NA)
    cr_data <- bind_rows(cr_data_eligible, cr_data_ineligible) %>% arrange(record_id, no_LIC_measurement)
    cr_data %>% view()

    # number of completed chart reviews n = 271 
    num_cr_completed <- cr_data %>% distinct_ct("record_id")

    # number of false positives n = 40
    num_false_positives <- cr_data %>% filter(is.na(transplant_date)) %>% distinct_ct("record_id") # 40
    num_tru_positives <- cr_data %>% filter(!is.na(transplant_date)) %>% distinct_ct("record_id") # 231
                            
    # completion by site
    cr_data %>% filter(iron_overload_chart_review_complete == 2) %>% 
        mutate(eligibility = case_when(eligibility == 1 ~ "Eligible", 
                                        eligibility == 0 ~ "Not eligible", 
                                        TRUE ~ "Missing")) %>%
        group_by(site, eligibility) %>% 
        summarise(num_completed = n()) %>% 
        ungroup() %>%
        pivot_wider(names_from = eligibility, values_from = num_completed)

    # chop data
    cr_chop <- cr_data %>% 
                filter(site == "chop") %>% distinct(record_id) %>% view()        

    # let's check who did not have a transplant
    # aka false positives
    cr_data %>% filter(is.na(transplant_date)) %>% 
        select(record_id) %>%
        inner_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>% collect(), by = c("record_id")) %>% 
        inner_join(results_tbl("study_cohorts") %>% collect(), by = c("person_id")) %>% 
        select(record_id, person_id, transplant_date, transplant_type, aim_2a_2:aim_3_2) %>%
        output_tbl(name = "false_positive_ids", local = TRUE, file = TRUE)
    
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
