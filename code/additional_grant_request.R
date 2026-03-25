
# additional grant application request
# Total number of patients with beta thalassemia (ok if no iron markers)
# n = 641

btm_tbl <- results_tbl("scd_dx") %>% 
    inner_join(load_codeset("beta_thal_major_dx"), by = c("scd_concept_id"= "concept_id")) 
    
btm_tbl %>% distinct_ct()

# Total number of patients with beta thalassemia who had stemcell transplant (again fine if no iron markers)

btm_transplant_tbl <- btm_tbl %>% 
    find_procedures("transplant_px") 

# Total number of patients with beta thalassemia who had an autologous/gene therapy transplant transplant (again fine if no iron markers)
btm_transplant_tbl %>% 
    filter(type == "autologous") %>%
    distinct_ct()

# Total number of patients with beta thalassemia who had an allogeneic transplant transplant (again fine if no iron markers)
btm_transplant_tbl %>% 
    filter(type == "allogeneic") %>%
    distinct_ct()

btm_transplant_tbl %>% 
    distinct_ct()
