.run  <- function() {

    setup_pkgs() # Load runtime packages as specified above

    message('Starting execution with framework version ',
    config('framework_version'))

    library(tidyr)

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
    master_xwalk_ids %>% output_tbl(name = "master_xwalk_ids_ALL_AIMS")

    for (s in 1:length(sites)){
      # get PEDSnet IDs
        master_xwalk_ids %>% filter(site == sites[s]) %>% mutate(MRN = "") %>%
                select(site, site_id = person_site_id, record_id, MRN) %>% 
                output_tbl(name = paste0("cr_ids/", sites[s], "_FORCHART_REVIEW_FULL_LIST"), local = TRUE, file = TRUE)
    }

    results_tbl("master_xwalk_ids_ALL_AIMS") %>% anti_join(results_tbl("study_cohorts"), by = "person_id") %>% 
            view()

}

