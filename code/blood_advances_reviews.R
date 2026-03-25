# Additional analyses required by reviewers
results_tbl(db = config("vocabulary"), in_schema("vocabulary", "concept")) %>%
 filter(grepl("transferrin|TFN", concept_name, ignore.case = TRUE)) %>%
 filter(grepl("saturation", concept_name, ignore.case = TRUE)) %>% 
 filter(domain_id == "Measurement") %>% 
 collect() %>% view()
 write.csv("specs/transferrin_saturation_labs.csv", row.names = FALSE)

cohort_covars <- results_tbl("analytics_dataset") 
LIC_data <- results_tbl("cr_LIC_data") %>% 
                  select(record_id, transplant_date, LIC, LIC_type, LIC_date, LIC_days_since_transplant, LIC_level) %>%  
                  # LIC_1, LIC_days_since_transplant_1, LIC_1_level, LIC_level_1_binary,
                  # LIC_3, LIC_days_since_transplant_3, LIC_3_level,
                  # LIC_max, LIC_3_max, LIC_days_since_transplant_max, LIC_days_since_transplant_3_max,
                  # LIC_max_level, LIC_3_max_level, LIC_level_3_binary, LIC_level_3_max_binary) %>%
                  collect()

cohort_covars %>%
  select(person_id, transplant_date) %>%
  inner_join(results_tbl("cdm_measurement") %>%
            inner_join(load_codeset("transferrin_saturation_labs"), by = c("measurement_concept_id" = "concept_id")), by = "person_id") %>% collect() %>% view()
  # filter(measurement_date <= transplant_date) %>% collect() %>%
  filter(abs(as.numeric(difftime(transplant_date, measurement_date, units = "days"))) <= 365) %>% view()
  
results_tbl("cdm_measurement") %>%
  filter(grepl("transferrin|TFN", measurement_source_value, ignore.case = TRUE) | grepl("transferrin", measurement_concept_name, ignore.case = TRUE)) %>%
  inner_join(cohort_covars %>% select(person_id, transplant_date), by = "person_id") %>%
  filter(measurement_date <= transplant_date) %>% collect() %>%
  filter(abs(as.numeric(difftime(transplant_date, measurement_date, units = "days"))) <= 365) %>% 
  group_by(measurement_concept_name) %>% 
  summarise(n = n_distinct(measurement_id))


# identify by bacteria type
# for all measurements prior to transplant
bacteremia_tbl <- results_tbl("covar_bacteremia_dx") %>% 
    mutate(bacteria_name = case_when(grepl("Enterobacter", organism_concept_name, ignore.case = TRUE) ~ "Enterobacter",
                                    grepl("Escherichia coli", organism_concept_name, ignore.case = TRUE) ~ "Escherichia coli",
                                    grepl("Staphylococcus", organism_concept_name, ignore.case = TRUE) ~ "Staphylococcus (including MRSA/MSSA)",
                                    grepl("Streptococcus", organism_concept_name, ignore.case = TRUE) ~ "Streptococcus & Related Genera",
                                    grepl("Salmonella", organism_concept_name, ignore.case = TRUE) ~ "Salmonella",
                                    grepl("Enterococcus", organism_concept_name, ignore.case = TRUE) ~ "Enterococcus",
                                    grepl("Corynebacterium|Coryneform|Bacillus|Clostridium|Actinomyces|Cutibacterium acnes", organism_concept_name, ignore.case = TRUE) ~ "Gram-Positive Bacilli",
                                    grepl("Acinetobacter|Achromobacter xylosoxidans|Rhizobium|Agrobacterium group|Klebsiella pneumoniae|
                                    Stenotrophomonas maltophilia|Citrobacter freundii|Enterobacter cloacae complex|
                                    Non-lactose fermenting gram-negative bacillus|Lactose fermenter|Pseudomonas oryzihabitans|
                                    putida|Pseudomonas fluorescens|Pseudomonas aeruginosa|Pseudomonas species not aeruginosa|
                                    Haemophilus influenzae|Serratia marcescens|Pantoea|Burkholderia gladioli|Citrobacter amalonaticus|
                                    Gram-negative bacillus", organism_concept_name, ignore.case = TRUE) ~ "Gram-Negative Bacilli",
                                    grepl("Gram-negative coccus|Moraxella catarrhalis|Neisseria", organism_concept_name, ignore.case = TRUE) ~ "Gram-negative cocci",
                                    TRUE ~ "Unspecified (no matching concept name)")) 
                                    
bacteremia_tbl %>%
    group_by(bacteria_name) %>%       
    summarise(n = n_distinct(person_id)) %>%
    view()

# bacteria counts by aims
stat_dataset <- results_tbl("survival_analysis_dataset") 

ferritin_data <- results_tbl("cr_ferritin_data") 

stat_dataset_bac <-  stat_dataset %>%
            left_join(ferritin_data, by = c("person_id", "transplant_date")) %>%       
                filter(chart_completion, eligibility == "Yes", !is.na(ferritin_pre)) %>%
                mutate(ferritin_pre = if_else(ferritin_pre %in% c("low", "moderate"), "low-moderate", "high")) %>%
                filter(!is.na(survival_time_5yr_censored)) %>%
                mutate(bacteremia = ifelse(is.na(bacteremia), FALSE, bacteremia))

stat_dataset_bac %>% 
    select(person_id, ferritin_pre, transplant_date) %>%
    left_join(bacteremia_tbl %>% 
            select(person_id, bacteria_name, bacteremia_date), by = "person_id") %>%
    filter(bacteremia_date >= transplant_date) %>%
    collect() %>% 
    mutate(ferritin_pre = factor(ferritin_pre, levels = c("low-moderate", "high"))) %>%
    group_by(person_id) %>%
    slice_min(bacteremia_date, n = 1, with_ties = FALSE) %>%
    ungroup() %>% group_by(ferritin_pre, bacteria_name) %>%
    summarise(n = n_distinct(person_id)) %>% ungroup() %>%
    pivot_wider(names_from = ferritin_pre, values_from = n, values_fill = 0) %>% view()
