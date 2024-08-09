.run <- function(){
    # this function is used to resolve patients that belong to multiple subgroups
    source("code.EAD.R")

    # load condition codesets
    dx <- load_codeset("aplastic_anemia_dx") %>% mutate(scd_type = "aa") %>%
        union(load_codeset("diamond_blackfan_anemia_dx") %>% mutate(scd_type = "dba")) %>% 
        union(load_codeset("beta_thal_major_dx") %>% mutate(scd_type = "bta")) %>%
        union(load_codeset("scd_dx") %>% mutate(scd_type = "scd"))

    # 747 patients including true positive patients from chart reviews + patients in aim 3 without chart reviews 
    subgroup_ct <- cohort_covars %>% #filter(chart_completion) %>% 
                distinct(person_id, record_id) %>% 
                copy_to_new(df=., name = "njh") %>%
                inner_join(results_tbl("cdm_condition_occurrence") %>% 
                                select(person_id, condition_concept_name, condition_concept_id, 
                                        condition_type_concept_id, condition_status_concept_id, condition_start_date), by = "person_id") %>%
                inner_join(dx, by = c("condition_concept_id" = "concept_id")) %>% collect()

    # 105 patients with more than 2 subgroups 
    subgroup_ct %>% group_by(person_id) %>% 
                summarise(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>% filter(subgroup_ct > 1) %>% distinct_ct()

    # after filtering by final diagnosis or in problem list, 90 patients
    # require dianogsis code to be final diagnosis or in problem list, down to 90 patients
    # 615 patients with 1 subgroup, so clearly some patients did not have a condition in the problem list
    subgroup_ct %>% filter(condition_status_concept_id == 4230359 |
                        condition_type_concept_id %in% c(2000000089, 2000000090, 2000000091)) %>% 
                group_by(person_id) %>% 
                summarise(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>% filter(subgroup_ct > 1) %>% distinct_ct()

    # For the 90 patients, let's just list out the unique diseases in the order of diagnosis
    multi_subgroup_tbl_confirmed <- subgroup_ct %>% filter(condition_status_concept_id == 4230359 |
                condition_type_concept_id %in% c(2000000089, 2000000090, 2000000091)) %>% 
                group_by(person_id) %>% 
                mutate(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>%
                filter(subgroup_ct > 1) %>% 
                arrange(desc(condition_start_date)) %>%
                group_by(person_id, scd_type) %>% mutate(subgroup_ct = n()) %>% ungroup() %>%
                distinct(person_id, scd_type, .keep_all = TRUE) %>%
                arrange(person_id, condition_start_date, condition_concept_name) %>% 
                select(person_id, condition_start_date, scd_type, subgroup_ct) 

    # tie cases that need Nora's review
    tie_ids_confirmed <- c(547751, 4272899, 6748949, 8247950, 8555392) 
    tie_ids_unconfirmed <- c(848566, 4141896, 4141896)

    # the rule is if the most recent ones have more codes, then we will use the most recent subgroup 
    # and most of the time, they have more diagnosis as well, n = 85 patients 
    multi_subgroup_confirmed_resolved <- multi_subgroup_tbl_confirmed %>% 
                                filter(!(person_id %in% tie_ids_confirmed)) %>%
                                group_by(person_id) %>%
                                slice_max(subgroup_ct, n = 1, with_ties = FALSE) %>% ungroup() %>%
                                select(person_id, scd_type) %>% distinct() %>% view()

    # now onto multiple subgroups that had no final diagnosis or in problem list
    multi_subgroup_tbl_unconfirmed <- subgroup_ct %>% 
                group_by(person_id) %>% 
                mutate(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>%
                filter(subgroup_ct > 1) %>% 
                arrange(desc(condition_start_date)) %>%
                group_by(person_id, scd_type) %>% mutate(subgroup_ct = n()) %>% ungroup() %>%
                distinct(person_id, scd_type, .keep_all = TRUE) %>%
                arrange(person_id, condition_start_date, condition_concept_name) %>% 
                select(person_id, condition_start_date, scd_type, subgroup_ct) %>%
                anti_join(multi_subgroup_resolved, by = "person_id") %>%
                filter(!(person_id %in% tie_ids_confirmed)) %>% view()

    # applying the same rule that the most recent ones have more codes, then we will use the most recent subgroup 
    # and most of the time, they have more diagnosis as well, n = 13 patients 
    multi_subgroup_unconfirmed_resolved <- multi_subgroup_tbl_unconfirmed %>% 
                                filter(!(person_id %in% tie_ids_unconfirmed)) %>%
                                group_by(person_id) %>%
                                slice_max(subgroup_ct, n = 1, with_ties = FALSE) %>% ungroup() %>%
                                select(person_id, scd_type) %>% distinct()

    # no take care of the unconfirmed cases
    subgroup_ct %>% #filter(condition_status_concept_id == 4230359 | condition_type_concept_id %in% c(2000000089, 2000000090, 2000000091)) %>% 
                filter(person_id %in% c(tie_ids_confirmed, tie_ids_unconfirmed)) %>%
                group_by(person_id, condition_concept_name) %>% 
                mutate(ct = n()) %>%
                slice_max(condition_start_date, n = 1, with_ties = FALSE) %>% ungroup() %>%
                select(person_id, condition_concept_name, scd_type, last_diagonsis_date = condition_start_date, ct) %>% 
                arrange(person_id, last_diagonsis_date, condition_concept_name, scd_type) %>% 
                output_tbl(name = "multi_subgroup_unresolved", file = TRUE, local = TRUE)
    
    multi_subgroup_confirmed_resolved %>% union(multi_subgroup_unconfirmed_resolved) %>%
        output_tbl(name = "multi_subgroup_resolved", file = TRUE, local = TRUE)

}
