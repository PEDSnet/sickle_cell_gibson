.run <- function(){
    # this function is used to resolve patients that belong to multiple subgroups
    source("code.EAD.R")

    # load condition codesets
    dx <- load_codeset("aplastic_anemia_dx") %>% mutate(scd_type = "aa") %>%
        union(load_codeset("diamond_blackfan_anemia_dx") %>% mutate(scd_type = "dba")) %>% 
        union(load_codeset("beta_thal_major_dx") %>% mutate(scd_type = "bta")) %>%
        union(load_codeset("scd_dx") %>% mutate(scd_type = "scd"))
    dx %>% output_tbl(name = "disgnosis_dx", file = TRUE, local = TRUE)

    # 725 true positive patients 
    subgroup_ct <- results_tbl("study_cohorts") %>% 
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        select(person_id, record_id) %>%
                        left_join(results_tbl("cr_data"), by = "record_id") %>%
                        filter(chart_completion, eligibility == "Yes") %>% 
                        distinct(person_id, record_id) %>% 
                        copy_to_new(df=., name = "njh") %>%
                        inner_join(results_tbl("cdm_condition_occurrence") %>% 
                                select(person_id, condition_concept_name, condition_concept_id, 
                                        condition_type_concept_id, condition_status_concept_id, condition_start_date), by = "person_id") %>%
                        inner_join(dx, by = c("condition_concept_id" = "concept_id")) %>% collect()

    # 620 patients with just 1 subgroup
    # 105 patients with more than 2 subgroups 
    one_generic_subgroup <- subgroup_ct %>% group_by(person_id) %>% 
                summarise(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>% filter(subgroup_ct == 1) 
    one_generic_subgroup %>% distinct_ct()

    more_than_one_generic_subgroup <- subgroup_ct %>% group_by(person_id) %>% 
                summarise(subgroup_ct = n_distinct(scd_type)) %>% ungroup() %>% filter(subgroup_ct > 1) 
    more_than_one_generic_subgroup %>% distinct_ct()
    
    # require dianogsis code to be final diagnosis or in problem list
    # 90 patients still have 2 subgroups and 599 patients have 1 subgroups 
    # that means 36 patients were missed out when we filtered out by problem list or final diagnosis
    one_final_subgroup <- subgroup_ct %>% filter(condition_status_concept_id == 4230359 |
                        condition_type_concept_id %in% c(2000000089, 2000000090, 2000000091)) %>% 
                group_by(person_id) %>% 
                summarise(no_subgroup = n_distinct(scd_type)) %>% ungroup() %>% filter(no_subgroup == 1) 
    one_final_subgroup %>% distinct_ct()
    one_final_subgroup %>% inner_join(one_generic_subgroup, by = "person_id") %>% distinct_ct()

    # For the 90 patients with more than 2 subgroups, 
    # list out the unique diseases in the order of diagnosis and number of times they were diagnosed
    # for that scd_type
    multi_subgroup_tbl_confirmed <- subgroup_ct %>% filter(condition_status_concept_id == 4230359 |
                condition_type_concept_id %in% c(2000000089, 2000000090, 2000000091)) %>% 
                group_by(person_id) %>% 
                mutate(no_subgroup = n_distinct(scd_type)) %>% ungroup() %>%
                filter(no_subgroup > 1) %>% 
                arrange(desc(condition_start_date)) %>%
                group_by(person_id, scd_type) %>% 
                mutate(no_dx_per_subgroup = n()) %>% ungroup() %>%
                group_by(person_id, scd_type) %>%
                mutate(first_diagonsis_date = min(condition_start_date, na.rm = TRUE)) %>%
                mutate(last_diagonsis_date = max(condition_start_date, na.rm = TRUE)) %>% ungroup() %>%
                distinct(person_id, scd_type, .keep_all = TRUE) %>%
                arrange(person_id, first_diagonsis_date, condition_concept_name) %>% 
                select(person_id, first_diagonsis_date, last_diagonsis_date, scd_type, no_subgroup, no_dx_per_subgroup) 
    
    multi_subgroup_tbl_confirmed %>% view()
    # tie cases that need Nora's review

    # the rule is if the most recent ones have more codes, then we will use the most recent subgroup 
    # and we resolve 85 patients (who also had more diagnosis)
    # 5 patients could not be resolved because the most recent diagnosis had fewer counts (need Nora to confirm)
    tie_ids_confirmed <- c(547751, 4272899, 6748949, 8247950, 8555392) 
    multi_subgroup_confirmed_resolved <- multi_subgroup_tbl_confirmed %>% 
                                filter(!(person_id %in% tie_ids_confirmed)) %>%
                                group_by(person_id) %>%
                                slice_max(no_dx_per_subgroup, n = 1, with_ties = FALSE) %>% ungroup() %>%
                                distinct(person_id, scd_type, .keep_all = TRUE) %>% view()

    # now onto multiple subgroups that had no final diagnosis or in problem list
    # n = 15
    multi_subgroup_tbl_unconfirmed <- subgroup_ct %>% 
                group_by(person_id) %>% 
                mutate(no_subgroup = n_distinct(scd_type)) %>% ungroup() %>%
                filter(no_subgroup > 1) %>% 
                arrange(desc(condition_start_date)) %>%
                group_by(person_id, scd_type) %>% 
                mutate(no_dx_per_subgroup = n()) %>% ungroup() %>%
                group_by(person_id, scd_type) %>%
                mutate(first_diagonsis_date = min(condition_start_date, na.rm = TRUE)) %>%
                mutate(last_diagonsis_date = max(condition_start_date, na.rm = TRUE)) %>% ungroup() %>%
                distinct(person_id, scd_type, .keep_all = TRUE) %>%
                arrange(person_id, first_diagonsis_date, condition_concept_name) %>% 
                select(person_id, first_diagonsis_date, last_diagonsis_date, scd_type, no_subgroup, no_dx_per_subgroup) %>%
                anti_join(multi_subgroup_confirmed_resolved, by = "person_id") %>%
                filter(!(person_id %in% tie_ids_confirmed)) %>% view()

                # group_by(person_id) %>% 
                # mutate(no_subgroup = n_distinct(scd_type)) %>% ungroup() %>%
                # filter(no_subgroup > 1) %>% 
                # arrange(desc(condition_start_date)) %>%
                # group_by(person_id, scd_type) %>% 
                # mutate(no_dx_per_subgroup = n()) %>% ungroup() %>%
                # distinct(person_id, scd_type, .keep_all = TRUE) %>%
                # group_by(person_id, scd_type) %>%
                # mutate(first_diagonsis_date = min(condition_start_date, na.rm = TRUE)) %>%
                # mutate(last_diagonsis_date = max(condition_start_date, na.rm = TRUE)) %>% ungroup() %>%

                # arrange(person_id, condition_start_date, condition_concept_name) %>% 
                # select(person_id, condition_start_date, scd_type, subgroup_ct) %>%
                # anti_join(multi_subgroup_resolved, by = "person_id") %>%
                # filter(!(person_id %in% tie_ids_confirmed)) %>% view()

    # applying the same rule that the most recent ones have more codes, then we will use the most recent subgroup 
    # and most of the time, they have more diagnosis as well, 
    # we resolved n = 13 patients and 2 cases need Nora's review
    tie_ids_unconfirmed <- c(848566, 4141896)
    multi_subgroup_unconfirmed_resolved <- multi_subgroup_tbl_unconfirmed %>% 
                                filter(!(person_id %in% tie_ids_unconfirmed)) %>%
                                group_by(person_id) %>%
                                slice_max(no_dx_per_subgroup, n = 1, with_ties = FALSE) %>% ungroup() %>%
                                select(person_id, scd_type) %>% distinct()

    # now take care of the unconfirmed cases
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
    # all good, now need to add the 7 cases that nora confirms 

}
