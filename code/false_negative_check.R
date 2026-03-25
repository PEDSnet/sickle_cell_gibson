# use dcc_pedsnet_v53 for this analysis

# find all potential procedures 
read.csv("local/false_negative_ids.csv") %>% copy_to_new(df =., name = "false_neg") %>%
    inner_join(cdm_tbl("person"), by = c("person_id" = "site_id", "site")) %>% 
    select(site_id = person_id, person_id = person_id.y, site) %>%
        inner_join(cdm_tbl("procedure_occurrence"), by = c("person_id", "site")) %>%
        filter(grepl("transplant", procedure_concept_name, ignore.case = TRUE)|
                grepl("transfusion", procedure_concept_name, ignore.case = TRUE),
                grepl("stem cell", procedure_concept_name, ignore.case = TRUE) | 
                grepl("stem-cell", procedure_concept_name, ignore.case = TRUE)) %>% 
        distinct(procedure_concept_id, procedure_concept_name, person_id) %>%
        view()

# only these 2 had transplants: c(4672304, 4906306)
read.csv("local/false_negative_ids.csv") %>% copy_to_new(df =., name = "false_neg") %>%
    inner_join(cdm_tbl("person"), by = c("person_id" = "site_id", "site")) %>% 
    select(site_id = person_id, person_id = person_id.y, site) %>%
        inner_join(cdm_tbl("procedure_occurrence"), by = c("person_id", "site")) %>%
        filter(grepl("transfusion", procedure_concept_name, ignore.case = TRUE),
                # procedure_concept_id %in% c("2000087", "2002367", "40486968"),
                !(person_id %in% c("4672304", "4906306"))) %>% 
        distinct(procedure_concept_id, procedure_concept_name, person_id) %>%
        view()

# filter out relevant px_codes
px <- c("2788068", "2788074", "2841327", "2885573", "2000087", 
        "2002367", "40486968", "1781154", "2813368", "2000087", "2002367", "40486968")

# 13 patients 
read.csv("local/false_negative_ids.csv") %>% copy_to_new(df =., name = "false_neg") %>%
    inner_join(cdm_tbl("person"), by = c("person_id" = "site_id", "site")) %>% 
    select(site_id = person_id, person_id = person_id.y, site) %>%
    inner_join(cdm_tbl("procedure_occurrence") %>% 
            filter(procedure_concept_id %in% px), by = c("person_id", "site")) %>%
            distinct(procedure_concept_id, procedure_concept_name) %>% view()

dx <- load_codeset("aplastic_anemia_dx") %>% 
        union(load_codeset("diamond_blackfan_anemia_dx")) %>% 
        union(load_codeset("beta_thal_major_dx")) 

# 1 patient is in our transplant cohort but with a missing disease code Beta thalassemia, 4278669 (we missed this)
cdm_tbl("condition_occurrence") %>% 
    filter(person_id == "7524422") %>% distinct(condition_concept_name, condition_concept_id) %>% view()


