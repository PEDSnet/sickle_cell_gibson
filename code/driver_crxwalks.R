.run  <- function() {

    setup_pkgs() # Load runtime packages as specified above

    message('Starting execution with framework version ',
    config('framework_version'))

    library(tidyr)

    message('Generating redcap ids for chart reviews') 
    sites <- results_tbl("study_cohorts") %>% distinct(site) %>% arrange(site) %>% pull(site)

    master_xwalk_ids <- results_tbl("master_xwalk_ids") %>% collect()

    # get PEDSnet IDs
    tbl <- results_tbl("study_cohorts") %>% collect() %>% # no multiple transplants
                filter(!(aim_3_1| aim_3_2)) %>% 
                distinct(person_id, site) %>%
                union(results_tbl("multi_transplant_cohort") %>% collect() %>% #patients with unverifiable transplant status
                        # filter(site == sites[s]) %>% #, (aim_2a_2| aim_2b_2| aim_3_1| aim_3_2)) %>% 
                        distinct(person_id, site)) %>%
                distinct(person_id, site)
    # new patients
     new_ids <- tbl %>% anti_join(master_xwalk_ids, by = "person_id") %>% select(person_id, site)
     site_ct <- master_xwalk_ids %>% group_by(site) %>% summarise(site_ct = n()) %>% ungroup()

    # map PEDSnet IDs to site IDs
    tbl <- new_ids %>% copy_to_new(df = ., name = "jbi") %>%
            inner_join(cdm_tbl("person") %>% select(person_id, person_site_id = site_id), by = "person_id") %>% 
            collect()
    tbl <- tbl %>% group_by(site) %>% 
                mutate(site_id = cur_group_id(), 
                        site_id = if_else(site_id < 10, paste0("0", site_id), as.character(site_id))) %>%
                inner_join(site_ct, by = "site") %>%
                group_by(site) %>%
                mutate(record_id = case_when((row_number()  + site_ct)< 10 ~ paste0(site_id, "_000", row_number() + site_ct, "_cr"),
                                    (row_number()  + site_ct) < 100 ~ paste0(site_id, "_00", row_number()+ site_ct, "_cr"),
                                    (row_number()  + site_ct) < 1000 ~ paste0(site_id, "_0", row_number()+ site_ct, "_cr"))) %>%
            select(person_id, person_site_id, record_id, site_id, site)
    master_xwalk_ids <- tbl %>% union(master_xwalk_ids) %>% arrange(site_id, record_id)
    master_xwalk_ids %>% view()
#     master_xwalk_ids %>% output_tbl(name = "master_xwalk_ids_ALL_AIMS")

    # generate crosswalks ids for chart review for sites
    for (s in 1:length(sites)){
      # get PEDSnet IDs
        master_xwalk_ids %>% filter(site == sites[s]) %>% mutate(MRN = "") %>%
                select(site, site_id = person_site_id, record_id, MRN) %>% 
                output_tbl(name = paste0("cr_ids/", sites[s], "_FORCHART_REVIEW_FULL_LIST"), local = TRUE, file = TRUE)
    }

#     results_tbl("master_xwalk_ids_ALL_AIMS") %>% anti_join(results_tbl("master_xwalk_ids"), by = "record_id") %>%
#             output_tbl(name = "master_xwalk_ids_NEW", local = TRUE, file = TRUE)

    # just for convenience: make a new study_cohorts table with record_id and aims including multitransplant patients
    rslt$study_cohorts <- results_tbl("study_cohorts_old") %>% 
        select(person_id, transplant_date, site, aim_2a_2:aim_3_2) %>%
        union(results_tbl("multi_transplant_cohort")) %>% collect() %>%
        distinct() %>% filter(!duplicated(as.character(person_id))) %>%
        copy_to_new(df = ., name = "fsdsd")

    rslt$study_cohorts <- rslt$study_cohorts %>% left_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>% 
                                                        select(person_id, record_id, site_id), by = "person_id")

    # add demographics info: 
    rslt$study_cohorts <- rslt$study_cohorts %>% left_join(results_tbl("cdm_person") %>% 
                                                            select(person_id, birth_date, 
                                                                gender = gender_concept_name, 
                                                                race = race_concept_name, 
                                                                ethnicity = ethnicity_concept_name), by = "person_id") 
    rslt$study_cohorts <- rslt$study_cohorts %>% left_join(results_tbl("scd_dx") %>% 
                                                            select(person_id, scd_concept_name, scd_type) %>%
                                                            distinct(person_id, .keep_all = TRUE), by = "person_id")
    rslt$study_cohorts <- rslt$study_cohorts %>% left_join(results_tbl("no_multi_transplant_px") %>% 
                                                            select(person_id, transplant_concept_name, transplant_date) %>%
                                                            distinct(person_id, .keep_all = TRUE), 
                                                            by =c("person_id", "transplant_date")) %>%
                                                mutate(transplant_type = case_when(grepl("Allogeneic", transplant_concept_name, ignore.case = TRUE) ~ "Allogeneic",
                                                                    grepl("Autologous", transplant_concept_name, ignore.case = TRUE) ~ "Autologous/Gene therapy",
                                                                    TRUE ~ "other"))
     rslt$study_cohorts <- rslt$study_cohorts %>% left_join(results_tbl("study_cohorts_old") %>% select(person_id) %>% 
                                                            mutate(multi_transplants = FALSE), by = "person_id") %>%
                                                mutate(multi_transplants = if_else(is.na(multi_transplants), TRUE, multi_transplants))
     rslt$study_cohorts %>% view()
     rslt$study_cohorts %>% distinct_ct()  
     rslt$study_cohorts %>% collect() %>% output_tbl(name = "study_cohorts") 

}
