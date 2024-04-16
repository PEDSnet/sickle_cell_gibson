.run  <- function() {

    setup_pkgs() # Load runtime packages as specified above

    message('Starting execution with framework version ',
    config('framework_version'))

    rslt <- list()  
    # transplant patients
    # type = "all" means both allogeneic and autologous transplants
    # this functions only return patients with fewer than 2 transplants 
    rslt$transplant_cohort <- find_transplants(procedure_codeset_name = "transplant_px") 
    rslt$transplant_cohort[[1]] %>% output_tbl(name = "no_multi_transplant_px", indexes = list('person_id'), local = FALSE, results_tag = FALSE)

    rslt$transplant_cohort <- results_tbl("transplant_px") %>% collect_new() %>% group_by(transplant_type) %>% summarise(n = n_distinct(person_id))
    append_sum(cohort = 'transplant_px',
                persons = distinct_ct(rslt$transplant_cohort)) 

    # filter out patients with multiple transplants that are not identifiable 
    # rslt$transplant_cohort <- rslt$transplant_cohort %>%
    #                                 anti_join(read.csv("local/multi_transplants_patid.csv"), by = c("person_id")) %>%
    #                                 inner_join(read.csv("local/transplant_patid.csv") %>% 
    #                                             mutate(transplant_date = as.Date(transplant_date)), 
    #                                             by = c("person_id", "procedure_occurrence_id"= "transplant_occurrence_id", "transplant_date"))
    
    # patients with transplant dates outside first and last
    tbl_1 <- tbl %>% inner_join(cdm_tbl("visit_occurrence") %>% select(person_id,
                                                                visit_occurrence_id, 
                                                                visit_start_date, 
                                                                visit_end_date), by = c("visit_occurrence_id", 
                                                                                    "person_id")) %>% 
                                                            filter(procedure_date < visit_start_date |
                                                                procedure_date > visit_end_date) %>% collect_new()

    # 47 cases with procedure dates outside visit durations are already in the transplant cohort
    tbl_1 %>% distinct(person_id) %>% inner_join(rslt$transplant_cohort %>% distinct(person_id), by = c("person_id")) %>% distinct(person_id) %>% view()

    # n = 1003 with procedure dates outside visit durations are NOT in the transplant cohort
    tbl_1 %>% anti_join(rslt$transplant_cohort %>% distinct(person_id), by = c("person_id")) %>% arrange(person_id) %>% 
        output_tbl(name = "transplant_outside_visit_duration_ids", local = TRUE, file = TRUE)
        
    # group by site
    # tbl_1 %>% anti_join(rslt$transplant_cohort %>% distinct(person_id), by = c("person_id")) %>% 
    #     group_by(site) %>% summarise(n = n_distinct(person_id)) %>% view()    

    # # group by number of procedures 
    # tbl_1 %>% anti_join(rslt$transplant_cohort %>% distinct(person_id), by = c("person_id")) %>% 
    #         group_by(person_id) %>% summarise(n = n_distinct(procedure_date)) %>% view() 
            
    # we need a special cohort for cases when:
    # patients had more than 2 transplants/transplant dates outside first and last visit
    # AND meeting the criteria for the other aims 
    # n_rows = 2044
    rslt$transplant_cohort[[2]] %>% anti_join(rslt$transplant_cohort[[1]], by = "person_id") %>% distinct(person_id, transplant_date) %>%
         union(read_csv("local/transplant_outside_visit_duration_ids.csv") %>% distinct(person_id, procedure_date) %>% rename(transplant_date = procedure_date)) %>%
         distinct(person_id, transplant_date) %>% output_tbl(name = "multi_transplant_px")

    rslt$aim_2a_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("mri_px"), 
                                                                procedure_date = mri_date,
                                                                end_date = as.Date("2023-06-30"),
                                                                is_pre = TRUE, 
                                                                cut_off = cut_off) %>% collect()

    rslt$aim_2a_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                            transplant_cohort = results_tbl("multi_transplant_px"), 
                                                            procedure_cohort = results_tbl("ferritin_lab"), 
                                                            procedure_date = ferritin_date,
                                                            end_date = as.Date("2023-06-30"),
                                                            is_pre = TRUE, 
                                                            cut_off = cut_off) %>% collect()
    

    rslt$aim_2b_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("mri_px"), 
                                                                procedure_date = mri_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off) %>% collect()    

    rslt$aim_2b_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("ferritin_lab"), 
                                                                procedure_date = ferritin_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off, min_days = 180) %>% collect()

    rslt$aim_3_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("therapeutic_phleb_px"), 
                                                                procedure_date = phleb_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off) %>% collect()

    # 2. Number of patients who were prescribed deferasirox or deferoxamine after date of transplant.
    rslt$aim_3_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                            transplant_cohort = results_tbl("multi_transplant_px"), 
                                                            procedure_cohort = results_tbl("drugs_rx"), 
                                                            procedure_date = drug_exposure_start_date,
                                                            end_date = as.Date("2022-12-31"),
                                                            is_pre = FALSE, 
                                                            cut_off = cut_off) %>% collect()

    col_names <- c("person_id", "transplant_date", "transplant_concept_id","transplant_concept_name",
                "transplant_type", "scd_occurrence_id", "scd_concept_id", "scd_concept_name",
                "scd_start_date", "scd_type", "scd_start_age_in_months", "site",
                "birth_date","gender","age_at_transplant")

    # n = 168
    rslt$multi_transplant_cohort <- rslt$aim_2a_1_cohort %>% select(person_id, transplant_date, site) %>% mutate(study_aim = "aim_2a_1") %>% 
                    union(rslt$aim_2a_2_cohort %>% select(person_id, transplant_date, site)%>% mutate(study_aim = "aim_2a_2")) %>%
                    union(rslt$aim_2b_1_cohort %>% select(person_id, transplant_date, site)%>% mutate(study_aim = "aim_2b_1")) %>%
                    union(rslt$aim_2b_2_cohort %>% select(person_id, transplant_date, site)%>% mutate(study_aim = "aim_2b_2")) %>%
                    union(rslt$aim_3_1_cohort %>% select(person_id, transplant_date, site)%>% mutate(study_aim = "aim_3_1")) %>%
                    union(rslt$aim_3_2_cohort %>% select(person_id, transplant_date, site)%>% mutate(study_aim = "aim_3_2")) %>%
                    mutate(val = TRUE)

    rslt$multi_transplant_cohort %>% pivot_wider(id_cols = c(person_id, transplant_date, site), names_from = study_aim, values_from = val, values_fill = FALSE) %>% 
                                    output_tbl(name = "multi_transplant_cohort")

    # group by site
    rslt$multi_transplant_cohort %>% 
        group_by(site) %>% summarise(n = n_distinct(person_id)) %>% view()

}