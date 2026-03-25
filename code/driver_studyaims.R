.run  <- function() {


    library(tidyr)
    setup_pkgs() # Load runtime packages as specified above

    message('Starting execution with framework version ',
          config('framework_version'))

  # Set up the step log with as many attrition columns as you need.
  # For example, this call sets up the log with a `persons` count that will be
  # required at each step.
  init_sum(cohort = 'Start', persons = 0)

  # By convention, accumulate execution results in a list rather than as
  # independent variables, in order to make returning the entire set easier
  rslt <- list()
  # Aim_2_1: 
  min_days = 0
  max_days = 365
  is_pre = TRUE
  # get ferritin measurements closet to the CE date
  # ce-date is the date of 1st transplant and age_at_ce is age at 1st transplant
  tbl <- results_tbl("study_cohorts") %>% filter(aim_2a_1 == TRUE | aim_2a_2 == TRUE) %>%
                pivot_longer(cols = c(aim_2a_1, aim_2a_2), names_to = "aim", values_to = "aim_value") %>%
                select(person_id, ce_date, site, age_at_ce, transplant_type, scd_type, gender, aim) %>%
                inner_join(results_tbl("covar_ferritin_mx") %>% 
                    select(person_id, ce_date, 
                            ferritin_date = measurement_date, 
                            ferritin_value = value_as_number, 
                            ferritin_level), 
                by = c("person_id", "ce_date")) %>% 
                filter(!is.na(ferritin_value)) %>%
                collect_new()
  if(is_pre) {
    tbl <- tbl %>% filter(ferritin_date <= ce_date- days(min_days), ferritin_date >= (ce_date - days(max_days)))
  } else {    
    tbl <- tbl %>% filter(ferritin_date >= (ce_date + days(min-days)), ferritin_date <= (ce_date + days(max_days)))
  }
    # get ferritin measurements closet to the CE date
    tbl <- tbl %>% mutate(ferritin_date = as.numeric(difftime(ferritin_date, ce_date, units = "days"))) %>%
                group_by(person_id, ce_date) %>%
                slice_min(abs(ferritin_date), n =1, with_ties = FALSE) 
    # get relapse status
    tbl_relapse <- tbl %>% left_join(results_tbl("no_multi_transplant_px") %>% collect_new() %>%
                        select(person_id, transplant_relapse_date = transplant_date), by = "person_id") %>%
                    filter(transplant_relapse_date > ce_date) %>% 
                    mutate(transplant_relapse = TRUE,
                           transplant_relapse_day = as.numeric(difftime(transplant_relapse_date, ce_date, units = "days"))) %>%
                    select(-transplant_relapse_date)
    # get death status
    tbl <- tbl %>% left_join(tbl_relapse %>% 
                            select(person_id, transplant_relapse, transplant_relapse_day), 
                            by = c("person_id", "ce_date"))
    # n = 663
    tbl <- tbl %>% left_join(results_tbl("covar_death") %>% collect_new() %>%
                            mutate(death = TRUE) %>%
                            distinct(person_id, death_date, death), by = "person_id") %>%
                    mutate(death_date = as.numeric(difftime(death_date, ce_date, units = "days")))                      
    
    # conditioning agents before transplant???
    tbl <- tbl %>% left_join(results_tbl("conditioning_agents_rx_px") %>% collect_new() %>%
                            mutate(conditioning_agent = TRUE) %>%
                            distinct(person_id, conditioning_agent, transplant_date), by = c("person_id", "ce_date" = "transplant_date"))
    # GVHD
    tbl <- tbl %>% left_join(results_tbl("covar_gvhd_dx") %>% collect_new() %>%
                            mutate(gvhd = TRUE) %>%
                            distinct(person_id, gvhd, transplant_date), by = c("person_id", "ce_date" = "transplant_date"))
    # get VOD status
    tbl <- tbl %>% left_join(results_tbl("covar_vod_dx") %>% collect_new() %>%
                            mutate(vod = TRUE) %>%
                            distinct(person_id, vod, transplant_date), by = c("person_id", "ce_date" = "transplant_date"))  
    # neutrophil engraftment
    ANC <- results_tbl("covar_anc_mx") %>% collect_new() %>%
            mutate(anc = TRUE) %>%
            distinct(person_id, anc, ce_date, anc_date) %>%
            select(person_id, anc, anc_date)

    # get blood culture
    tbl %>% view()

  output_sum(name = "outcome_covariates_attrition", local = TRUE, file = TRUE)  

  message('Done.')

  invisible(rslt)

}