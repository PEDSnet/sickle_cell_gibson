.run() <- function(){
    setup_pkgs()

    message("get immune variables from lab data")
    rslt$immune <- results_tbl("study_cohorts") %>% distinct(person_id, transplant_date) %>% 
                inner_join(cdm_tbl("measurement_labs"), by = "person_id") %>% 
                filter(grepl("CD3", measurement_concept_name, ignore.case = TRUE) | 
                grepl("CD4", measurement_concept_name, ignore.case = TRUE) |
                grepl("CD8", measurement_concept_name, ignore.case = TRUE) |
                grepl("IgM", measurement_concept_name, ignore.case = TRUE) |
                grepl("lymphocyte", measurement_concept_name, ignore.case = TRUE)) %>%
                select(person_id, measurement_concept_id, measurement_concept_name, value_as_number, 
                        measurement_date, unit_concept_name, range_high, range_low) %>%
                collect()
   rslt$immune %>% output_tbl(name = "covar_immune_mx")

    rslt$WBC <- results_tbl("study_cohorts") %>% distinct(person_id, transplant_date) %>% 
                inner_join(cdm_tbl("measurement_labs"), by = "person_id") %>% 
                filter(grepl("white", measurement_concept_name, ignore.case = TRUE), 
                grepl("blood cell", measurement_concept_name, ignore.case = TRUE),
                !is.na(value_as_number)) %>%
                select(person_id, measurement_concept_id, measurement_concept_name, value_as_number, 
                        measurement_date, unit_concept_name, range_high, range_low) %>%  
                collect() 
    rslt$WBC <- rslt$WBC %>% mutate(WBC_ct = if_else(value_as_number < 100, value_as_number * 1000, value_as_number)) 

  # acordding to google, CD3 is broken down to CD4 and CD8 
  # check units and add type of measurement
  # the units are mostly cells per microliter, some are cells per cubic millimeter
  
  message("get transplant dates from chart review and exclude false positive patients")
  rslt$study_cohorts <- results_tbl("study_cohorts") %>% 
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        left_join(results_tbl("cr_data") %>% select(record_id, transplant_date_cr = transplant_date, 
                                                                transplant_type_cr = transplant_type, 
                                                                second_transplant_date_cr = second_transplant_date,
                                                                disease_relapse = graft_fail, chart_completion, eligibility) %>% distinct(), by = "record_id") %>% 
                        mutate(transplant_date = if_else(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                                transplant_type = if_else(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                                chart_completion = if_else(is.na(chart_completion), FALSE, chart_completion)) %>%
                        filter((eligibility == "Yes" & chart_completion) | !chart_completion) %>%
                        select(-transplant_date_cr, -transplant_type_cr)

  # get dates of 2nd transplants (if any)
  # For patients with available chart reviews, use chart reivew transplant dates as index dates
  # Filter out false positives from chart reviews
  rslt$transplant_px <- results_tbl("study_cohorts") %>%
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        left_join(results_tbl("cr_data") %>% 
                        select(record_id, record_id, transplant_date_cr = transplant_date, 
                                transplant_type_cr = transplant_type, 
                                second_transplant_date_cr = second_transplant_date,
                                chart_completion, eligibility) %>% distinct(), by = "record_id") %>% 
                        mutate(transplant_date = if_else(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                                transplant_type = if_else(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                                chart_completion = if_else(is.na(chart_completion), FALSE, chart_completion)) %>%
                        filter((eligibility == "Yes" & chart_completion) | !chart_completion) %>%
                        select(-transplant_date_cr, -transplant_type_cr)
  rslt$transplant_px %>% view()

  rslt$transplant_px <- rslt$transplant_px %>% inner_join(results_tbl("no_multi_transplant_px") %>% 
                        select(person_id, second_transplant_date = transplant_date), by = "person_id") %>% collect_new() %>%
                        filter(second_transplant_date >= transplant_date) %>% 
                        mutate(second_transplant_date = if_else(second_transplant_date == transplant_date, NA, second_transplant_date)) %>%
                        group_by(person_id) %>%
                        slice_min(second_transplant_date, n = 1, with_ties = FALSE, na_rm = FALSE) %>% 
                        ungroup()

  # use second transplant dates from chart review 
  rslt$transplant_px <- rslt$transplant_px %>% 
                        mutate(second_transplant_date = case_when(!is.na(second_transplant_date_cr) ~ second_transplant_date_cr, 
                                                                (is.na(second_transplant_date_cr) & chart_completion) ~ NA, #overwrite second transplant dates from EHR
                                                                TRUE ~ second_transplant_date)) %>% 
                        select(-second_transplant_date_cr)

  rslt$CD3 <- results_tbl("covar_immune_mx") %>% 
                filter(grepl("CD3", measurement_concept_name, ignore.case = TRUE),
                        measurement_concept_id %in% c("3010993", "3011412", "3014686", "3018064", "3022533", "3025632", "4209239"),
                        !grepl("CD4", measurement_concept_name, ignore.case = TRUE), 
                        !grepl("CD8", measurement_concept_name, ignore.case = TRUE)) %>% collect()

  rslt$CD3_valid <- rslt$CD3 %>% filter(range_low > 100 | (is.na(range_low) & unit_concept_name == "per microliter"),
                                        (value_as_number >= 10), !is.na(value_as_number)) 
  
  rslt$lymphocytes <- results_tbl("covar_immune_mx") %>% 
                                    filter(grepl("lymphocyte", measurement_concept_name, ignore.case = TRUE),
                                        !grepl("kappa", measurement_concept_name, ignore.case = TRUE),
                                        !grepl("lambda", measurement_concept_name, ignore.case = TRUE),
                                        measurement_concept_id %in% c("3003215", "3004327", "3019198", "4254663")) %>% 
                                mutate(conversion_factor = case_when(unit_concept_name == "billion per liter" ~ 1000,
                                                unit_concept_name == "per microliter" ~ 1,
                                                unit_concept_name == "thousand per microliter" ~ 1000,
                                                unit_concept_name == "Kelvin per microliter" ~ 1000,
                                                (unit_concept_name == "No information" | unit_concept_name == "No matching concept") & measurement_concept_id== "3004327" & range_low <100 ~ 1000,
                                                unit_concept_name == "cells per microliter" ~ 1,
                                                (unit_concept_name == "percent" | is.na(unit_concept_name)) & measurement_concept_id== "3004327" ~ -99,
                                                unit_concept_name == "No information" & is.na(range_low) & measurement_concept_id== "3019198" ~ -99,
                                                unit_concept_name == "No matching concept" & range_low < 100 & measurement_concept_id== "3019198" ~ 1000,
                                                unit_concept_name == "microliter" ~ 1,
                                                unit_concept_name == "per cubic millimeter" ~ 1,
                                                unit_concept_name == "thousand per cubic millimeter" ~ 1000,
                                                unit_concept_name == "percent" & measurement_concept_id== "3019198" ~ 100,
                                                unit_concept_name == "percent" & measurement_concept_id== "4254663" ~ 100,
                                                TRUE ~ NA),
                                       lymphocyte_ct = value_as_number * conversion_factor) %>%
                                filter(!is.na(lymphocyte_ct), lymphocyte_ct > 0)
                                # group_by(measurement_concept_id, measurement_concept_name,  unit_concept_name, range_high, range_low) %>% summarise(n = n()) %>% view()
 
  # need a lymphocyte counts
  rslt$CD3_invalid <- rslt$CD3 %>% filter(unit_concept_name %in% c("No information", "No matching concept", "percent") | is.na(unit_concept_name),
                                        range_high < 100 | is.na(range_high),
                                        range_low < 100, measurement_concept_id != "4210706") %>% 
                                mutate(conversion_factor = if_else(measurement_concept_id == "3011412", 100, NA)) 
                                # group_by(measurement_concept_id, measurement_concept_name,  unit_concept_name, range_high, range_low) %>% summarise(n = n()) %>% view()
  
  #multiply the CD3 perc by lymphocyte counts to get CD3 counts
  rslt$CD3_invalid <- rslt$CD3_invalid %>% filter(measurement_concept_id == "3022533" | measurement_concept_id == "3014686") %>%
                        rename(CD3_perc = value_as_number) %>% 
                        inner_join(rslt$lymphocytes %>% select(person_id, 
                                                                lymphocyte_ct, 
                                                                lymphocyte_date = measurement_date) %>% collect(), by = "person_id", relationship = "many-to-many") %>%
                        filter(abs(as.numeric(difftime(measurement_date, lymphocyte_date, units = "days"))) <= 2) %>% 
                        mutate(CD3_ct = round(CD3_perc * lymphocyte_ct / 100)) %>%
                        distinct(person_id, measurement_date, .keep_all = TRUE) %>%
                        select(-c("CD3_perc", "conversion_factor", "lymphocyte_ct", "lymphocyte_date")) %>%
                        union(rslt$CD3_invalid %>% filter(measurement_concept_id == "3011412") %>%
                                                    mutate(CD3_ct = conversion_factor*value_as_number) %>%
                                                    select(-value_as_number, -conversion_factor))

  rslt$CD3 <- rbind(rslt$CD3_valid %>% rename(CD3_ct = value_as_number), rslt$CD3_invalid) 
  rslt$CD3 %>% output_tbl(name = "covar_CD3_mx")

   rslt$CD4 <- rslt$immune %>% filter(grepl("CD4", measurement_concept_name, ignore.case = TRUE),
                                measurement_concept_id %in% c("4209260", "21494815")) 
   # need a lymphocyte counts
   rslt$CD4 <- rslt$CD4 %>% filter(measurement_concept_id %in% c("21494815")) %>%
                rename(CD4_perc = value_as_number) %>%
                inner_join(rslt$lymphocytes %>% select(person_id, 
                                                                lymphocyte_ct, 
                                                                lymphocyte_date = measurement_date), by = "person_id", relationship = "many-to-many") %>%
                        filter(abs(as.numeric(difftime(measurement_date, lymphocyte_date, units = "days"))) <= 2,
                                !is.na(CD4_perc)) %>% 
                        mutate(CD4_ct = round(CD4_perc * lymphocyte_ct / 100)) %>%
                        distinct(person_id, measurement_date, .keep_all = TRUE) %>% 
                        select(-c("CD4_perc", "lymphocyte_ct", "lymphocyte_date")) %>%
                        union(rslt$CD4 %>% filter(measurement_concept_id == "4209260", !is.na(value_as_number)) %>%
                                            mutate(CD4_ct = value_as_number) %>%
                                            select(-value_as_number))
   rslt$CD4 %>% output_tbl(name = "covar_CD4_mx")                                
                                   
   rslt$CD8 <- results_tbl("covar_immune_mx") %>% filter(grepl("CD8", measurement_concept_name, ignore.case = TRUE), 
                                measurement_concept_id %in% c("21494814", "44807257"))                     
                                # group_by(measurement_concept_id, measurement_concept_name,  unit_concept_name, range_high, range_low) %>% summarise(n = n()) %>% view()
   # need a lymphocyte counts
    rslt$CD8 <- rslt$CD8 %>% filter(measurement_concept_id %in% c("21494814")) %>%
                rename(CD8_perc = value_as_number) %>%
                inner_join(rslt$lymphocytes %>% select(person_id, 
                                                                lymphocyte_ct, 
                                                                lymphocyte_date = measurement_date), by = "person_id") %>% collect() %>%
                        filter(abs(as.numeric(difftime(measurement_date, lymphocyte_date, units = "days"))) <= 2,
                                !is.na(CD8_perc)) %>% 
                        mutate(CD8_ct = round(CD8_perc * lymphocyte_ct / 100)) %>%
                        distinct(person_id, measurement_date, .keep_all = TRUE) %>% 
                        select(-c("CD8_perc", "lymphocyte_ct", "lymphocyte_date")) %>%
                        union(rslt$CD8 %>% filter(measurement_concept_id == "44807257", !is.na(value_as_number)) %>%
                                            mutate(CD8_ct = value_as_number) %>%
                                            select(-value_as_number) %>% collect())
    rslt$CD8 %>% output_tbl(name = "covar_CD8_mx")                                       
    
    rslt$IgM <- rslt$immune %>% filter(grepl("IgM", measurement_concept_name, ignore.case = TRUE), 
                                measurement_concept_id %in% c("3003039", "3028026")) 
    # need a lymphocyte counts
    rslt$IgM <- rslt$IgM %>% filter(measurement_concept_id %in% c("3003039")) %>%
                rename(IgM_perc = value_as_number) %>%
                inner_join(rslt$lymphocytes %>% select(person_id, 
                                                                lymphocyte_ct, 
                                                                lymphocyte_date = measurement_date), by = "person_id", relationship = "many-to-many") %>%
                        filter(abs(as.numeric(difftime(measurement_date, lymphocyte_date, units = "days"))) <= 2,
                                !is.na(IgM_perc)) %>% 
                        mutate(IgM_ct = round(IgM_perc * lymphocyte_ct / 100)) %>%
                        distinct(person_id, measurement_date, .keep_all = TRUE) %>% 
                        select(-c("IgM_perc", "lymphocyte_ct", "lymphocyte_date")) %>%
                        union(rslt$IgM %>% filter(measurement_concept_id == "3028026", !is.na(value_as_number)) %>%
                                            mutate(IgM_ct = value_as_number) %>%
                                            select(-value_as_number))                   
                                # group_by(measurement_concept_id, measurement_concept_name,  unit_concept_name, range_high, range_low) %>% summarise(n = n()) %>% view()
    rslt$IgM %>% output_tbl(name = "covar_IgM_mx")

}