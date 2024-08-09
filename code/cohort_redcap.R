#' This function update the most recent redcap data
#' filter out old patient ids, extract completed records and push it to the database
get_redcap_data <- function(redcap_filename){

    # pull the most recent version of the data
    cr_data <- read.csv(redcap_filename)
    # filter out old patient ids
    cr_data <- cr_data %>% filter(!grepl("^[A-Za-z]", record_id))
    # only stanford records
    # stanford_cr <- cr_data %>% filter(redcap_data_access_group == "stanford", eligibility == 1) %>% 
    #       select(record_id, iron_val_1:iron_val_10) %>% 
    #       mutate(across(iron_val_1:iron_val_10, ~as.character(.x))) %>%
    #       pivot_longer(cols = starts_with("iron_val_"), names_to = "no_LIC_measurement1", 
    #                         values_to = "LIC_raw", values_drop_na = TRUE) %>% 
    #       separate(LIC_raw, into = c("start", "end"), sep = "-", remove = FALSE) %>%
    #       mutate(start = as.numeric(start), 
    #              end = as.numeric(end)) %>% 
    #       mutate(LIC_raw = rowMeans(select(., c("start", "end")), na.rm = TRUE)) %>% 
    #       select(-c(start, end)) %>%
    #       pivot_wider(names_from = no_LIC_measurement1, values_from = LIC_raw, values_fill = NA) %>%
    #       view()

    # ok have to take care of this "" as well

    # extract completed records 
    cr_data <- cr_data %>% rename(site = redcap_data_access_group, 
                        chart_completion = iron_overload_chart_review_complete,
                        transplant_date = date_tpx,
                        transplant_type = tpx_type,
                        LIC_other_unit = lic_other_unit) %>%
                mutate(transplant_type = case_when(transplant_type == 1 ~ "Allogeneic", 
                                            transplant_type == 2 ~ "Autologous/Gene therapy", 
                                            TRUE ~ NA),
                        chart_completion = if_else(chart_completion == 2, TRUE, FALSE), 
                        across(iron_date_1:est_method_15, ~if_else(.x == "", NA, .x)),
                        across(c("transplant_date", "second_transplant_date", "graft_fail_date"), ~if_else(.x == "", NA, as.Date(.x, "%Y-%m-%d"))),
                        across(iron_val_1:iron_val_15, ~as.character(.x)),
                        # transplant_date = if_else(transplant_date == "", NA, as.Date(transplant_date)),
                        eligibility = case_when(is.na(transplant_date ) & is.na(eligibility) ~ 0,
                                                !is.na(transplant_date ) & is.na(eligibility) ~ 1,  
                                                TRUE ~ eligibility),
                        eligibility = if_else(eligibility ==0, "No", "Yes"),
                        donor_relation = case_when( donor_relation == 1 ~ "Related",
                                                    donor_relation == 2 ~ "Unrelated",
                                                    TRUE ~ NA),
                        match_status = case_when(match_status == 1 ~ "9/10",
                                                match_status == 2 ~ "10/10",
                                                match_status == 4 ~ "8/10",
                                                match_status == 3 ~ "Haploidentical",
                                                match_status == 5 ~ "Cord Blood",
                                                TRUE ~ NA),
                        manip_type = case_when(manip_type == 1 ~ "Alpha/Beta T cell depletion",
                                               manip_type == 2 ~ "Post-transplant cyclophosphmide",
                                               TRUE ~ NA),
                        graft_manip = case_when(graft_manip == 0 ~ "No", 
                                                graft_manip == 1 ~ "Yes", 
                                                TRUE ~ NA),
                        graft_fail = case_when(graft_fail == 0 ~ "No", 
                                                graft_fail == 1 ~ "Yes", 
                                                TRUE ~ NA),
                        LIC_other_unit = if_else(LIC_other_unit %in% c(" ", "n/a", "N/A", "na", ""), NA, LIC_other_unit)) 

    # DQ check
    # after this step, everything should be fixxed                    
    cr_dq_check(cr_data)

    # filter out incomplete charts
    cr_data <- cr_data %>% filter(chart_completion)
    # any eligible cases without any LIC measurements would be filtered out
    cr_data_eligible <- cr_data %>% pivot_longer(cols = c(starts_with("iron_date_")), names_to = "no_LIC_measurement", 
                            values_to = "LIC_date", values_drop_na = TRUE) %>% 
                mutate(no_LIC_measurement = as.numeric(sub(".*_([0-9]+)$", "\\1", no_LIC_measurement)),
                    LIC_date = as.Date(LIC_date)) %>% 
                pivot_longer(cols = starts_with("iron_val_"), names_to = "no_LIC_measurement1", 
                            values_to = "LIC_raw", values_drop_na = TRUE) %>%
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
                            mutate(LIC_raw = NA, LIC_date = NA, LIC_est_method = NA, no_LIC_measurement = NA)

    cr_data <- bind_rows(cr_data_eligible, cr_data_ineligible) %>% arrange(record_id, no_LIC_measurement) %>%
                mutate(LIC_est_method = case_when(LIC_est_method == 1 ~ "R2",
                                                LIC_est_method == 2 ~ "R2*",
                                                LIC_est_method == 3 ~ "T2*",
                                                TRUE ~ NA)) 

    # only stanford records come in ranges                                            
    cr_data <- cr_data %>% separate(LIC_raw, into = c("start", "end"), sep = "-", remove = FALSE) %>%
                mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
                mutate(LIC_raw = rowMeans(select(., c("start", "end")), na.rm = TRUE)) %>% 
                select(-c(start, end))
                

    return(cr_data)
} 

cr_dq_check <- function(cr_data){
    # count incomplete charts
    incomplete_cr <- cr_data %>% mutate(t = rowSums(!is.na(select(., -c("record_id", "site"))))) %>%  filter(t>0, !chart_completion, eligibility == "Yes") %>% 
        select(record_id) %>% pull()
        if(length(incomplete_cr) == 0){
            print("All completed charts are marked as complete")
        } else {
            print(paste0(length(incomplete_cr), " incomplete chart"))
            for (i in incomplete_cr){
            print(paste0(i, " incomplete chart"))}
        }

    # count chart with missing LIC dates or methods    
    missing_LIC_info <- cr_data %>% mutate(t = rowSums(!is.na(select(., iron_date_1:est_method_15)))) %>%  
                        filter((t %% 3) != 0) %>% select(record_id) %>% pull()
    if(length(missing_LIC_info) == 0){
            print("All completed charts have all LIC info")
        } else {
            print(paste0(length(missing_LIC_info), " charts with missing LIC info"))
            for (i in missing_LIC_info){
            print(paste0(i, " missing LIC info"))}
        }                    
}