setwd("~/Projects/sickle_cell_Gibson/")
source("site/run.R")

.run  <- function() {

    setup_pkgs() # Load runtime packages as specified above

    message('Starting execution with framework version ',
    config('framework_version'))

    library(tidyr)
    library(lubridate)
    # Set up the step log with as many attrition columns as you need.
    # For example, this call sets up the log with a `persons` count that will be
    # required at each step.
    init_sum(cohort = 'Start', persons = 0)

    # By convention, accumulate execution results in a list rather than as
    # independent variables, in order to make returning the entire set easier
    rslt <- list()

    start_date <- as.Date("2010-01-01")
    end_date <- as.Date("2023-12-31")
    cut_off = 365 #days 

    # patients with relevant non-malignant blood disorders include: sickle cell disease, beta thalassemia major, aplastic anemia, Diamond-Blackfan anemia (congenital pure red cell aplasia)
    rslt$scd_cohort <- find_scd_conditions() %>% collect() 
    append_sum(cohort = 'scd_dx',
                persons = distinct_ct(rslt$scd_cohort))
    rslt$scd_cohort %>% output_tbl(name = "scd_dx", indexes = list('person_id'), local = FALSE, results_tag = FALSE)
    
    rslt$mri_cohort <- find_procedures("mri_px") %>% 
                        rename(mri_concept_id = procedure_concept_id, 
                               mri_concept_name = procedure_concept_name, 
                               mri_date = procedure_date, 
                               mri_site = site) %>% collect_new()
    append_sum(cohort = 'mri_px',
                persons = distinct_ct(rslt$mri_cohort))
    rslt$mri_cohort %>% output_tbl(name = "mri_px", indexes = list('person_id'), local = FALSE, results_tag = FALSE)

    rslt$phleb_cohort <- find_procedures("therapeutic_phleb_px") %>% 
                        rename(phleb_concept_id = procedure_concept_id, 
                               phleb_concept_name = procedure_concept_name, 
                               phleb_date = procedure_date, 
                               phleb_site = site) %>% collect_new()
    append_sum(cohort = 'therapeutic_phleb_px',
                persons = distinct_ct(rslt$phleb_cohort))
    rslt$phleb_cohort %>% output_tbl(name = "therapeutic_phleb_px", indexes = list('person_id'), local = FALSE, results_tag = FALSE)

    # find patients with deferasirox or deferoxamine
    dx_codeset<- load_codeset("defibrotide_rx") %>% mutate(type = "defibrotide") %>% #defibrotide
                union(load_codeset("deferoxamine_rx") %>% mutate(type = "deferoxamine")) %>% #deferoxamine
                union(load_codeset("deferasirox_rx") %>% mutate(type = "deferasirox")) %>% #deferasirox
                compute_new(temp = TRUE, name = "drug_id")
    rslt$drugs_cohort <- find_drugs(dx_codeset) %>% collect()
    append_sum(cohort = 'drugs_rx',
                persons = distinct_ct(rslt$drugs_cohort))
    rslt$drugs_cohort %>% output_tbl(name = "drugs_rx", indexes = list('person_id'), local = FALSE, results_tag = FALSE)

    # find patients with pre-transplant ferritin (within 1 year before the date of transplant)
    mx_ferritin <- load_codeset("ferritin_lab") %>% compute_new()
    rslt$ferritin_cohort <- find_measurements(mx_ferritin) %>% collect()
    append_sum(cohort = 'ferritin_lab',
                persons = distinct_ct(rslt$ferritin_cohort))
    rslt$ferritin_cohort %>% output_tbl(name = "ferritin_lab", indexes = list('person_id'), local = FALSE, results_tag = FALSE)

    # find patients with pre-transplant abdominal or liver MRI (within 1 year)
    rslt$aim_2a_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("mri_px"), 
                                                                procedure_date = mri_date,
                                                                end_date = as.Date("2023-06-30"),
                                                                is_pre = TRUE, 
                                                                cut_off = cut_off) %>% collect()
    append_sum(cohort = 'aim_2a_1', persons = distinct_ct(rslt$aim_2a_1_cohort))
    rslt$aim_2a_1_cohort %>% output_tbl(name = "aim_2a_1", indexes = list('person_id'))

    # Aim 2a:
    # 2. Number of patients who have received a bone marrow transplant AND had a non-malignant blood disorder AND had pre-transplant ferritin (within 1 year before the date of transplant).
    # Same disease criteria. 
    rslt$aim_2a_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                            transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                            procedure_cohort = results_tbl("ferritin_lab"), 
                                                            procedure_date = ferritin_date,
                                                            end_date = as.Date("2023-06-30"),
                                                            is_pre = TRUE, 
                                                            cut_off = cut_off) %>% collect()
    append_sum(cohort = 'aim_2a_2', persons = distinct_ct(rslt$aim_2a_2_cohort))
    rslt$aim_2a_2_cohort %>% output_tbl(name = "aim_2a_2", indexes = list('person_id'))

    # Aim 2b:   
    # 1. Number of patients who have received a bone marrow transplant AND have a non-malignant blood disorder AND have a post-transplant abdominal or liver MRI (any date after the date of transplant). 
    #       a. Relevant non-malignant blood disorders include: sickle cell disease, beta thalassemia major, aplastic anemia, Diamond-Blackfan anemia
    #       b. Both allogeneic and autologous transplants should be included.
    # 2. Number of patients who have received a bone marrow transplant AND a non-malignant blood disorder AND with a post-transplant ferritin (within 1 year after the date of transplant). Same disease criteria.
    # for patients with more than 1 transplant, keep the earliest one
    rslt$aim_2b_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("mri_px"), 
                                                                procedure_date = mri_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off) %>% collect()
    append_sum(cohort = 'aim_2b_1', persons = distinct_ct(rslt$aim_2b_1_cohort))
    rslt$aim_2b_1_cohort %>% output_tbl(name = "aim_2b_1", indexes = list('person_id'))

    rslt$aim_2b_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("ferritin_lab"), 
                                                                procedure_date = ferritin_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off, min_days = 180)  %>% collect()
    append_sum(cohort = 'aim_2b_2', persons = distinct_ct(rslt$aim_2b_2_cohort))
    rslt$aim_2b_2_cohort %>% output_tbl(name = "aim_2b_2", indexes = list('person_id'))

    # Aim 3: 
    # 1. Number of patients who had a phlebotomy procedure after date of transplant.
    rslt$aim_3_1_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                                transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                                procedure_cohort = results_tbl("therapeutic_phleb_px"), 
                                                                procedure_date = phleb_date,
                                                                end_date = as.Date("2022-12-31"),
                                                                is_pre = FALSE, 
                                                                cut_off = cut_off) %>% collect()
    append_sum(cohort = 'aim_3_1', persons = distinct_ct(rslt$aim_3_1_cohort))
    rslt$aim_3_1_cohort %>% output_tbl(name = "aim_3_1", indexes = list('person_id'))

    # 2. Number of patients who were prescribed deferasirox or deferoxamine after date of transplant.
    rslt$aim_3_2_cohort <- find_scd_transplant_with_procedure(scd_cohort = results_tbl("scd_dx"), 
                                                            transplant_cohort = results_tbl("no_multi_transplant_px"), 
                                                            procedure_cohort = results_tbl("drugs_rx"), 
                                                            procedure_date = drug_exposure_start_date,
                                                            end_date = as.Date("2022-12-31"),
                                                            is_pre = FALSE, 
                                                            cut_off = cut_off) %>% collect()
    append_sum(cohort = 'aim_3_2', persons = distinct_ct(rslt$aim_3_2_cohort))
    rslt$aim_3_2_cohort %>% output_tbl(name = "aim_3_2", indexes = list('person_id'))

    output_sum(name = "attrition", local = TRUE, file = TRUE)

}
