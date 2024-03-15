# Vector of additional packages to load before executing the request
config_append('extra_packages', c())

#' Execute the request
#'
#' This function presumes the environment has been set up, and executes the
#' steps of the request.
#'
#' In addition to performing queries and analyses, the execution path in this
#' function should include periodic progress messages to the user, and logging
#' of intermediate totals and timing data through [append_sum()].
#'
#' @return The return value is dependent on the content of the request, but is
#'   typically a structure pointing to some or all of the retrieved data or
#'   analysis results.  The value is not used by the framework itself.
#' @md
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
  col_names <- c("person_id", "transplant_date", "transplant_concept_id","transplant_concept_name",
                "transplant_type", "scd_occurrence_id", "scd_concept_id", "scd_concept_name",
                "scd_start_date", "scd_type", "scd_start_age_in_months", "site",
                "birth_date","gender","age_at_transplant")

  # already double checked, for patients included in more than one aim, the transplant date is the same
  # and the condition_start_date is the same
  # there is one patinet person_id = 4820766 with 2 death causes that might cause problem 
  # n = 531
  rslt$cohort <- results_tbl("aim_2a_1", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_2a_1") %>%
                    union(results_tbl("aim_2a_2", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_2a_2")) %>%
                    union(results_tbl("aim_2b_1", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_2b_1")) %>%
                    union(results_tbl("aim_2b_2", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_2b_2")) %>%
                    union(results_tbl("aim_3_1", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_3_1")) %>%
                    union(results_tbl("aim_3_2", results_tag = FALSE) %>% select(all_of(col_names)) %>% mutate(study_aim = "aim_3_2")) %>%
                    collect_new() %>%
                    mutate(ce_date = as.Date(transplant_date, format = "%Y-%m-%d"), 
                            val = TRUE,
                            study_aim = factor(study_aim)) %>% 
                    rename(age_at_ce = age_at_transplant) 

  # n = 849
  rslt$cohort <- rslt$cohort %>% pivot_wider(id_cols = c(person_id, ce_date), names_from = study_aim, values_from = val, values_fill = FALSE) %>% 
                inner_join(rslt$cohort %>% select(-study_aim, -val, -ce_date) %>% distinct(person_id, .keep_all = TRUE), by = c("person_id"))
  rslt$cohort %>% view()
  append_sum(cohort = 'Cohort',
             persons = rslt$cohort %>% distinct_ct())
  rslt$cohort %>% output_tbl(name = "study_cohorts")
  
  message('Get relevant data using the predefined cohorts')

  message("CD, CD4, CD8 measurements")
  rslt$CD348 <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("cd3_cd4_cd8_mx"),
                                        is_pre = FALSE, is_closet = FALSE) 
#   rslt$CD348 %>% view()
  # check units and add type of measurement
  # the units are mostly cells per microliter, some are cells per cubic millimeter
  # exclude the ones with percent as unit
#   rslt$CD348 %>% mutate(CDtype = case_when(#measurement_concept_id %in% c("3011412", "3022533") ~ "CD3",
#                                 grepl("T4 helper", measurement_concept_name, ignore.case = TRUE) ~ "CD4",
#                                 grepl("T8 suppressor", measurement_concept_name, ignore.case = TRUE) ~ "CD8",
#                                 grepl("CD3 cells", measurement_concept_name, ignore.case = TRUE) ~ "CD3")) %>%
#                 filter(!is.na(value_as_number), unit_concept_name != "percent") %>%
#                 group_by(measurement_concept_id, measurement_concept_name, unit_concept_name, CDtype, range_high, range_low) %>% summarise(n = n()) %>% view()         

  rslt$CD348 %>% mutate(CDtype = case_when(grepl("T4 helper", measurement_concept_name, ignore.case = TRUE) ~ "CD4",
                                grepl("T8 suppressor", measurement_concept_name, ignore.case = TRUE) ~ "CD8",
                                grepl("CD3 cells", measurement_concept_name, ignore.case = TRUE) ~ "CD3")) %>%
                                filter(!is.na(value_as_number), unit_concept_name != "percent") %>%
                                output_tbl(name = "covar_CD348_mx")
  append_sum(cohort = 'CD3-4-8_counts',
             persons = rslt$CD348 %>% distinct_ct())

  message("IgM measurements")
  rslt$IgM <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("IgM_mx"),
                                        is_pre = FALSE, is_closet = FALSE) 
#   rslt$IgM %>% view()
  # check units
#   rslt$IgM %>% group_by(measurement_concept_id, measurement_concept_name, unit_concept_name) %>% summarise(n = n()) %>% view()
  # everything matches to measurement_concept_id = 3028026 get units milligram per deciliter
  # for measurement_concept_id = 4164130 only take milligram per deciliter
  rslt$IgM <- rslt$IgM %>% filter(measurement_concept_id %in% c(3028026, 4164130)) %>%
                        output_tbl(name = "covar_IgM_mx")
  append_sum(cohort = 'IgM',
                 persons = rslt$IgM %>% filter(measurement_concept_id %in% c(3028026, 4164130)) %>% distinct_ct())

   # pre-transplant ferritin    
   rslt$ferritin_pre <- get_outcome_measurements(cohort = results_tbl("study_cohorts") %>% filter(aim_2a_2 == TRUE | aim_3_2 == TRUE), 
                                                                mx_codeset = load_codeset("ferritin_lab"),
                                                                is_pre = TRUE, is_closet = FALSE) 
   rslt$ferritin_post <- get_outcome_measurements(cohort = results_tbl("study_cohorts") %>% filter(aim_2b_2 == TRUE), 
                                                                mx_codeset = load_codeset("ferritin_lab"),
                                                                is_pre = FALSE, is_closet = FALSE) 
   rslt$ferritin <- rbind(rslt$ferritin_pre, rslt$ferritin_post)
#    rslt$ferritin %>% view()

   # everything matches to measurement_concept_id = 3004121, majority with units nanogram per milliliter
   # assume missing units are nanogram per milliliter
   # exclude 1 case with NA unit measurement_concept_id = 4176561
#    rslt$ferritin %>% filter(measurement_concept_id == "3015242") %>% view()
   rslt$ferritin %>% filter(measurement_concept_id != "4176561") %>% 
                mutate(ferritin_level = case_when(value_as_number < 1000 ~ "low",
                                                value_as_number >= 1000 & value_as_number <= 3000 ~ "moderate",
                                                value_as_number > 3000 ~ "high")) %>%
                output_tbl(name = "covar_ferritin_mx")                                                                
   append_sum(cohort = 'covar_ferritin_px',
                     persons = rslt$ferritin %>% distinct_ct())

   # atleast 60 days after ce_date
  rslt$creatinine <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("serum_creatinine"),
                                        is_pre = FALSE, is_closet = FALSE, min_days = 180) 
#   rslt$creatinine %>% view()
  # check units, everything matches to measurement_concept_id = 3016723, units are milligram per deciliter
#   rslt$creatinine %>% group_by(measurement_concept_id, measurement_concept_name, unit_concept_name) %>% summarise(n = n()) %>% view()
  rslt$creatinine %>% filter(!is.na(value_as_number)) %>% output_tbl(name = "covar_creatinine_mx")
  append_sum(cohort = 'serum_creatinine',
                     persons = rslt$creatinine %>% distinct_ct())
  # VOD indicators
  rslt$defibrotide <- get_outcome_drugs(cohort = results_tbl("study_cohorts"), 
                                        dx_codeset = load_codeset("defibrotide_rx"),
                                        is_pre = FALSE, is_closet = FALSE) 
#   rslt$defibrotide %>% group_by(drug_concept_name, route_concept_name) %>% 
#                 summarise(n = n()) %>% view()
  rslt$abdominal_ultrasound <- get_outcome_procedures(cohort = results_tbl("study_cohorts"),
                                                px_codeset = load_codeset("abdominal_ultrasound_px"),
                                                is_pre = FALSE, is_closet = FALSE)
  rslt$bilirubin <- get_outcome_measurements(cohort = results_tbl("study_cohorts"),
                                                mx_codeset = load_codeset("bilirubin_mx"),
                                                is_pre = FALSE, is_closet = FALSE)
   #select patients with defibrotide or (abdominal ultrasound and bilirubin > 2mg/L 3 days apart)
   rslt$defibrotide %>% output_tbl(name = "defibrotide_rx")
   rslt$abdominal_ultrasound %>% output_tbl(name = "abdominal_ultrasound_px")
   rslt$bilirubin %>% output_tbl(name = "bilirubin_mx")
  
  #select patients with defibrotide or (abdominal ultrasound and bilirubin > 2mg/L 3 days apart)
  # expecting many to many relationship because of multiple ultrasound per day and multiple bilirubin measurements per day
  rslt$bilirubin %>% filter(value_as_number> 2) %>%
                        arrange(person_id, measurement_date) %>% select(person_id, transplant_date, bilirubin = value_as_number, 
                                                bilirubin_date = measurement_date) %>% 
                                                distinct(person_id, transplant_date,  bilirubin_date) %>% 
                        inner_join(rslt$abdominal_ultrasound %>% arrange(person_id, procedure_date) %>% 
                                distinct(person_id, transplant_date, abd_ultrasound_date = procedure_date), 
                                by = c("person_id", "transplant_date"), relationship = "many-to-many") %>% 
                                filter(abs(as.numeric(difftime(abd_ultrasound_date, bilirubin_date, units = "days"))) <= 3) %>%
                                group_by(person_id, transplant_date) %>% 
                                slice_min(bilirubin_date, n = 1, with_ties = FALSE) %>%
                        full_join(rslt$defibrotide %>% 
                                select(person_id, transplant_date, defibrotide_date = drug_exposure_start_date) %>% 
                                group_by(person_id, transplant_date) %>%
                                slice_min(abs(as.numeric(difftime(defibrotide_date, transplant_date, units = "days"))), n = 1, with_ties = FALSE),
                                by = c("person_id", "transplant_date"), relationship = "many-to-many") %>% 
                                mutate(vod = TRUE) %>% 
                                output_tbl(name = "covar_vod_dx")

  append_sum(cohort = 'defibrotide_drugs',
             persons = rslt$defibrotide %>% distinct_ct())

#   rslt$blood_culture <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
#                                                 mx_codeset = load_codeset("blood_culture_mx"),
#                                                 is_pre = FALSE, is_closet = FALSE) 
#   # add organism name
#   rslt$blood_culture %>% distinct(measurement_concept_id, measurement_concept_name) %>% 
#                         output_tbl(name = "blood_culture_mx", file = TRUE, local = TRUE)
#                         summarise(n = n()) %>% view()
#   rslt$blood_culture %>% filter(grepl("detected", value_as_concept_name, ignore.case = TRUE),
#                                 !grepl("not", value_as_concept_name, ignore.case = TRUE)) %>%
#                         mutate(bactereamia = case_when(grepl("Enterobacter", value_as_concept_name, ignore.case = TRUE) ~ "enterobacteria",
#                                                         grepl("pseudomonas aeruginosa", value_as_concept_name, ignore.case = TRUE) ~ "pseudomonas aeruginosa",
#                                                         grepl("gram negative", value_as_concept_name, ignore.case = TRUE) ~ "gram negative",
#                                                         grepl("Staphylococcus", value_as_concept_name, ignore.case = TRUE) ~ "staphylococcus",)) %>% select(person_id, bactereamia) %>%
#                                 view()
#   append_sum(cohort = 'blood_culture',
#              persons = rslt$blood_culture %>% distinct_ct())
 
  # atleast 60 days after ce_date
  rslt$ALT <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                mx_codeset = load_codeset("ALT_mx"),
                                is_pre = FALSE, is_closet = FALSE, min_days = 180) 
  append_sum(cohort = 'ALT',
                 persons = rslt$ALT %>% distinct_ct())
  # check units
  rslt$ALT %>% filter(!is.na(value_as_number), measurement_concept_id != "4146380") %>% 
                output_tbl(name = "covar_ALT_mx")
                # group_by(measurement_concept_id, measurement_concept_name, unit_concept_name) %>% 
                # summarise(n = n(), min = min(value_as_number), median = median(value_as_number),
                #                    max = max(value_as_number)) %>% view()
   
   # check why didnt some patients have height and heightz? what was missing? 
   rslt$heightz_pre <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("heightz_vital"),
                                        is_pre = TRUE, is_closet = FALSE) 
   rslt$heightz_post <- rslt$heightz_pre
   rslt$heightz_post <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("heightz_vital"),
                                        is_pre = FALSE, is_closet = FALSE) 
   rslt$heightz <- rbind(rslt$heightz_pre, rslt$heightz_post) %>% output_tbl(name = "covar_heightz")
   append_sum(cohort = 'heightz',
                 persons = rslt$heightz %>% distinct_ct())
  
#    rslt$height <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
#                                         mx_codeset = load_codeset("height_vital"),
#                                         is_pre = FALSE, is_closet = FALSE) 

#    append_sum(cohort = 'height',
#                  persons = rslt$height %>% distinct_ct())
    
   rslt$gvhd <- get_outcome_conditions(cohort = results_tbl("study_cohorts"), 
                                                dx_codeset = load_codeset("gvhd_dx"),
                                                is_pre = FALSE, is_closet = FALSE) 
#    rslt$gvhd %>% group_by(condition_concept_id, condition_concept_name) %>% 
#                 summarise(n = n()) %>% view()
   rslt$gvhd %>% output_tbl(name = "covar_gvhd_dx")
   append_sum(cohort = 'gvhd',
                  persons = rslt$gvhd %>% distinct_ct())

   rslt$conditioning_dx <- get_outcome_drugs(results_tbl("study_cohorts"), 
                                        dx_codeset = load_codeset("conditioning_drugs_rx"),
                                        is_pre = TRUE, is_closet = TRUE) 

#    rslt$conditioning_dx %>% view()     
   append_sum(cohort = 'conditioning_rx',
                    persons = rslt$conditioning_dx %>% distinct_ct())

   rslt$conditioning_px <- get_outcome_procedures(results_tbl("study_cohorts"),
                                                px_codeset = load_codeset("conditioning_tbt_px"),
                                                is_pre = TRUE, is_closet = TRUE) 
#     rslt$conditioning_px %>% view()
    append_sum(cohort = 'conditioning_px',
                      persons = rslt$conditioning_px %>% distinct_ct())

    rslt$conditioning_dx %>% select(person_id, transplant_date, 
                                concpet_id = drug_concept_id, 
                                concept_name = drug_concept_name, 
                                exposure_date = drug_exposure_start_date, 
                                route = route_concept_name) %>% mutate(type = "drugs") %>%
                             union(rslt$conditioning_px %>% select(person_id, 
                                                                transplant_date,
                                                                concpet_id = procedure_concept_id,
                                                                concept_name = procedure_concept_name,
                                                                exposure_date = procedure_date) %>% mutate(type = "drugs", route = NA)) %>%
                             output_tbl("conditioning_agents_rx_px")
    # neutrophil counts
    rslt$ANC <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                        mx_codeset = load_codeset("neutrophils_mx"),
                                        is_pre = FALSE, is_closet = FALSE)   
    # check units
    rslt$ANC %>% output_tbl(name = "covar_anc_mx_original")
    rslt$ANC %>% filter(measurement_concept_id %in% c("3013650", "3017501", "3017732", "3038972")) %>% 
        mutate(value_as_number = case_when(unit_concept_name == "Kelvin per microliter" ~ value_as_number * 1000,
                                        unit_concept_name == "thousand per microliter" ~ value_as_number*1000,
                                        unit_concept_name == "thousand per cubic millimeter" ~ value_as_number*1000,
                                        unit_concept_name == "billion per liter" ~ value_as_number,
                                        unit_concept_name == "cells per microliter" ~ value_as_number,
                                        unit_concept_name == "per microliter" ~ value_as_number,
                                        unit_concept_name == "microliter" ~ value_as_number,
                                        unit_concept_name == "per cubic millimeter" ~ value_as_number,
                                        TRUE ~ NA)) %>%
        filter(!is.na(value_as_number)) %>% output_tbl(name = "covar_anc_mx")
    # calculate ANC engraftment date
    # anc_date is the first of 3 consecutive days with ANC >= 500 after nadir
    nadir <- 7 
    results_tbl("covar_anc_mx") %>% collect_new() %>% filter(!is.na(value_as_number)) %>%
                group_by(person_id, transplant_date) %>%
                arrange(person_id, transplant_date, measurement_date) %>%
                mutate(anc_day = as.numeric(difftime(measurement_date, transplant_date, units = "days")),
                       val_ct = ifelse(value_as_number >= 500, 1, 0)) %>% 
                mutate(anc_consec = ifelse((lead(anc_day, 1) - anc_day == 1) & (lead(anc_day, n =2) - anc_day== 2) &
                                          (val_ct + lead(val_ct, 1) + lead(val_ct, n =2)) == 3, TRUE, FALSE)) %>% 
                ungroup() %>% 
                filter(anc_consec == TRUE, anc_day >= nadir) %>%
                group_by(person_id, transplant_date) %>%
                slice_min(measurement_date, n = 1, with_ties = FALSE) %>%
                rename(anc_engraftment_date = measurement_date) %>% 
                select(-anc_day, -val_ct) %>%
                output_tbl("covar_anc_engraftment_mx")

   rslt$platelet <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
                                                mx_codeset = load_codeset("platelet_mx"),
                                                is_pre = FALSE, is_closet = FALSE) 
        #   rslt$platelet %>% view()   
        # check units
        # thousand per microliter = per cubic millimeter = Kelvin per microliter (this K actually means thousands)
        # for measurement_concept_id = 3007461 with different ranges, it seems like we could ignore them because 
        # some of the high ranges are not relevant. e.g no values between 300 and 750
        # exclude the entitic volume as it measures the volume rather than the counts

        #   results_tbl("covar_platelet_mx") %>% group_by(measurement_concept_id, 
        #                                 measurement_concept_name, 
        #                                 unit_concept_name,
        #                                 range_high, range_low) %>% summarise(n = n()) %>% view()
  rslt$platelet %>% output_tbl(name = "covar_platelet_mx_original")
  rslt$platelet %>% filter(measurement_concept_id %in% c( "3007461", "3024929", "4267147"),
                          !is.na(value_as_number)) %>% 
                    mutate(value_as_number = value_as_number * 1000, unit_concept_name = "per microliter") %>%
                    output_tbl(name = "covar_platelet_mx")
                                
  append_sum(cohort = 'platelet',
                 persons = rslt$platelet %>% distinct_ct())  

    # platelet transfusion
    rslt$platelet_transfusion_pre <- get_outcome_procedures(results_tbl("study_cohorts"),
                                                        px_codeset = load_codeset("platelet_transfusion_px"),
                                                        is_pre = TRUE, is_closet = FALSE)
    rslt$platelet_transfusion_pre %>% output_tbl(name = "covar_platelet_transfusion_pre_px")
    rslt$platelet_transfusion %>% output_tbl(name = "covar_platelet_transfusion_post_px")
    results_tbl("covar_platelet_transfusion_pre_px") %>% union(results_tbl("covar_platelet_transfusion_post_px")) %>% 
                collect_new() %>% filter(!is.na(procedure_date)) %>%
                output_tbl(name = "covar_platelet_transfusion_px")
    rslt$platelet_transfusion %>% view()
    # how many transplant patients who never had a platelet transfusion
    # cohort_ct = 849
    # n = 426 who never had an transfusion??? 
    no_transfusion_id <- results_tbl("study_cohorts") %>% 
                distinct(person_id, transplant_date) %>% 
                anti_join(results_tbl("covar_platelet_transfusion_px") %>% 
                            distinct(person_id, transplant_date), by = c("person_id", "transplant_date")) %>%
                            compute_new()
    # did we miss sth?
    cdm_tbl("procedure_occurrence") %>% inner_join(no_transfusion_id, by = "person_id") %>% 
                filter(procedure_date <= transplant_date + days(180), procedure_date >= transplant_date - days(180),
                        grepl("transfusion", procedure_concept_name, ignore.case = TRUE)) %>%
                # distinct(person_id, transplant_date, procedure_concept_id, procedure_concept_name, .keep_all = TRUE) %>%
                select(person_id, transplant_date, procedure_concept_id, procedure_concept_name, procedure_date) %>% 
                arrange(person_id, transplant_date, procedure_date) %>% view()
                
    # calculate platelet engraftment date
    
    # get the first of 3 consecutive days with a platelet count of 20,000/uL or higher
    # when there were multiple measurements on the same day, we take the mean
    # but then Nora suggested taking the max values
     tbl <- results_tbl("covar_platelet_mx") %>% collect_new() %>% 
                group_by(person_id, transplant_date, measurement_date) %>%
                summarise(value_as_number = max(value_as_number, na.rm = TRUE)) %>%
                mutate(platelet_day = as.numeric(difftime(measurement_date, transplant_date, units = "days")),
                       val_ct = ifelse(value_as_number >= 20000, 1, 0)) %>% 
                mutate(platelet_consec = ifelse((lead(platelet_day, 1) - platelet_day == 1) & (lead(platelet_day, n =2) - platelet_day== 2) &
                                          (val_ct + lead(val_ct, 1) + lead(val_ct, n =2)) == 3, TRUE, FALSE)) %>% 
                ungroup() %>% 
                filter(platelet_consec == TRUE) %>%
                select(-val_ct)

      # The first of 3 consecutive days with a platelet count of 20,000/uL or higher 
      # in the absence of platelet transfusion for 7 consecutive days.
#       tbl %>% inner_join(results_tbl("covar_platelet_transfusion_px") %>% 
#                                 distinct(person_id, transplant_date, procedure_date) %>% 
#                                 collect_new() %>%
#                                 # ignore the platelet transfusion that happened outside of 100 days within the transplant date
#                                 filter(abs(as.numeric(difftime(procedure_date, transplant_date, units = "days"))) <= 100),
#                         by = c("person_id", "transplant_date")) %>% 
#                         filter(as.numeric(difftime(measurement_date, procedure_date, units = "days")) >= 7) %>%
                        
#                         # slice_min(measurement_date, n = 1, with_ties = FALSE) %>%
#                         select(person_id, transplant_date, 
#                                 platelet_transfusion_date = procedure_date, 
#                                 platelet_engraftment_date = measurement_date) %>%
#                         group_by(person_id, transplant_date) %>%
#                         slice_max(platelet_transfusion_date, n = 1, with_ties = TRUE) %>% ungroup() %>%
#                         group_by(person_id, transplant_date, platelet_transfusion_date) %>%
#                         slice_min(platelet_engraftment_date, n = 1, with_ties = FALSE) %>%
#                         mutate(platelet_engraftment_day = as.numeric(difftime(platelet_engraftment_date, transplant_date, units = "days"))) %>%
#                         # filter(as.numeric(difftime(platelet_engraftment_date, transplant_date, units = "days")) >= 7) %>%
#                         arrange(person_id, transplant_date, platelet_engraftment_date, platelet_transfusion_date) %>%
#                         view()
        t1 <- tbl %>% left_join(results_tbl("covar_platelet_transfusion_px") %>% 
                                distinct(person_id, transplant_date, procedure_date) %>% 
                                collect_new() %>%
                                mutate(transfusion_days_since_transplant = as.numeric(difftime(procedure_date, transplant_date, units = "days"))) %>% 
                                # ignore the platelet transfusion that happened outside of 100 days within the transplant date
                                filter(abs(transfusion_days_since_transplant) <= 100),
                        by = c("person_id", "transplant_date")) %>% 
                        mutate(#days_since_transplant = as.numeric(difftime(platelet_engraftment_date, transplant_date, units = "days")),
                                days_since_transfusion = as.numeric(difftime(measurement_date, procedure_date, units = "days"))) %>%
                        filter(days_since_transfusion >= 7) %>%
                        group_by(person_id, transplant_date) %>%
                        slice_min(days_since_transfusion, n = 1, with_ties = FALSE) %>%
                        slice_min(platelet_day, n = 1, with_ties = FALSE) %>%
                        select(person_id, transplant_date, 
                                last_transfusion_date = procedure_date, 
                                platelet_engraftment_date = measurement_date,
                                days_since_transplant = platelet_day,
                                days_since_last_transfusion = days_since_transfusion) %>%
                        arrange(person_id, transplant_date, platelet_engraftment_date, last_transfusion_date) %>% collect()

        # case when engraftment date is the same as the transplant date
        t1 %>% filter(platelet_engraftment_date == transplant_date) %>% 
                distinct(person_id, transplant_date, .keep_all = TRUE) %>% 
                inner_join(results_tbl("covar_platelet_mx") %>% collect_new(), by = c("person_id", "transplant_date")) %>%              
                mutate(platelet_day = as.numeric(difftime(measurement_date, transplant_date, units = "days")),
                       val_ct = ifelse(value_as_number >= 20000, 1, 0)) %>% 
                mutate(platelet_consec = ifelse((lead(platelet_day, 1) - platelet_day == 1) & (lead(platelet_day, n =2) - platelet_day== 2) &
                                          (val_ct + lead(val_ct, 1) + lead(val_ct, n =2)) == 3, TRUE, FALSE)) %>% 
                ungroup() %>% filter(platelet_consec == TRUE) %>%
                select(person_id, transplant_date, measurement_date, 
                        value_as_number, platelet_day, platelet_consec, 
                        last_transfusion_date, days_since_transplant, days_since_last_transfusion, 
                        platelet_engraftment_date) %>% view()

                # output_tbl("covar_platelet_engraftment_mx")

        # need to double check the very high engraftment_day, must be due to error with max(procedure_date)
        # it look like the engraftment date for this patient is the date of the transplant
        # the reason is that this patient did not have any platelet transfusion before on immediately after the transplant
        # thus when slicing from the days since last transfusion, it will slice at the next transfusion which is 10 months later
        # 2 things could have gone wrong: wrong unit conversion for the platelet count 
        # or miss out codes from the platelet transfusion
        id = "1544646"
        tbl %>% filter(person_id == id) %>% #platelet counts
                arrange(person_id, transplant_date, measurement_date) %>%
        view()
        results_tbl("covar_platelet_transfusion_px") %>% filter(person_id == id) %>% view()
        tbl %>% left_join(results_tbl("covar_platelet_transfusion_px") %>% 
                                distinct(person_id, transplant_date, procedure_date) %>% 
                                collect_new() %>%
                                mutate(transfusion_days_since_transplant = as.numeric(difftime(procedure_date, transplant_date, units = "days"))) %>% 
                                # ignore the platelet transfusion that happened outside of 100 days within the transplant date
                                filter(abs(transfusion_days_since_transplant) <= 100),
                        by = c("person_id", "transplant_date")) %>% 
                        mutate(#days_since_transplant = as.numeric(difftime(platelet_engraftment_date, transplant_date, units = "days")),
                                days_since_transfusion = as.numeric(difftime(measurement_date, procedure_date, units = "days"))) %>%
                        filter(days_since_transfusion >= 7, person_id == id) %>%
                        group_by(person_id, transplant_date) %>%
                        slice_min(days_since_transfusion, n = 1, with_ties = FALSE) %>%
                        slice_min(platelet_day, n = 1, with_ties = FALSE) %>%
                        select(person_id, transplant_date, 
                                last_transfusion_date = procedure_date, 
                                platelet_engraftment_date = measurement_date,
                                days_since_transplant = platelet_day,
                                days_since_last_transfusion = days_since_transfusion) %>%
                        arrange(person_id, transplant_date, platelet_engraftment_date, last_transfusion_date) %>%
                        view()
    # who had a another transplant after the first one?
#     rslt$cohort %>% inner_join(results_tbl("transplant_px"), by = c("person_id")) %>% 
#                 group_by(person_id) %>% 
#                 filter(transplant_date > ce_date) %>% view()
#                 summarise(relapse = n()) %>% filter(relapse > 1) %>% view()
    # transplant counts: 
    # how many transplants did each person have?
    # find out the visit dates for each transplant
#     visit_durations <- results_tbl("transplant_px") %>% inner_join(cdm_tbl("procedure_occurrence"), 
#                                 by = c("person_id", "transplant_occurrence_id" = "procedure_occurrence_id")) %>%
#                         select(person_id, transplant_occurrence_id, transplant_date, visit_occurrence_id, transplant_concept_name) %>%
#                         inner_join(cdm_tbl("visit_occurrence") %>% select(visit_occurrence_id,
#                                                                           visit_start_date,
#                                                                           visit_end_date), by = c("visit_occurrence_id")) %>% 
#                         distinct(person_id, transplant_date, .keep_all = TRUE) %>%
#                         collect_new()
                                                
    # find all patients with at least 2 transplants, filter out case when procedures outside visit windows 
    # group by visit, if the difference between 2 transplants is 1 day, then we'll take the 2nd one as the ce_date assuming that the first one
    # was cancelled
    # n = 604
#     t <- rslt$cohort %>% mutate(transplant_date = as.Date(transplant_date, format = "%Y-%m-%d")) %>%
#                 select(-transplant_date) %>%
#                 inner_join(visit_durations %>% filter(transplant_date >= visit_start_date, 
#                         transplant_date <= visit_end_date), by = c("person_id")) %>% 
#                 group_by(person_id) %>%
#                 summarise(transplant_count = n_distinct(transplant_date))

#      # n = 549
#      single_transplants_patid <- t %>% filter(transplant_count == 1) %>% 
#                 inner_join(visit_durations %>% filter(transplant_date >= visit_start_date, 
#                         transplant_date <= visit_end_date), by = c("person_id")) %>% 
#                 distinct(person_id, transplant_date, transplant_occurrence_id)
                 
#     # when the gap between 2 transplants is a day, i think we could say that this is the same transplant and the transplant date is the later one
#     # n = 9
#     single_transplants_dup_patid <- rslt$cohort %>% mutate(transplant_date = as.Date(transplant_date, format = "%Y-%m-%d")) %>%
#                 select(-transplant_date) %>%
#                 inner_join(visit_durations %>% filter(transplant_date >= visit_start_date, 
#                         transplant_date <= visit_end_date), by = c("person_id")) %>% 
#                 group_by(person_id) %>% 
#                 mutate(next_transplant_date = lead(transplant_date, order_by = transplant_date), 
#                        next_transplant_occurrence_id = lead(transplant_occurrence_id, order_by = transplant_date),
#                         transplant_count = n_distinct(transplant_date)) %>%
#                 filter(transplant_count == 2, next_transplant_date - transplant_date <= 2) %>%
#                 select(person_id, transplant_date = next_transplant_date, transplant_occurrence_id = next_transplant_occurrence_id)
     
#      # when the gap between 2 transplants is more than a day but there were only 2 transplants, it's ok and we'll take the 1st transqplant date
#      # as ce_date and the 2nd one as the relapse date, 
#      # n = 38 
#      # n = 10 cases when both transplants happened in the same visit, already checked, not duplicated transplants
#      two_transplants_patid <- rslt$cohort %>% mutate(transplant_date = as.Date(transplant_date, format = "%Y-%m-%d")) %>%
#                 select(-transplant_date) %>%
#                 inner_join(visit_durations %>% filter(transplant_date >= visit_start_date, 
#                         transplant_date <= visit_end_date), by = c("person_id")) %>% 
#                 group_by(person_id) %>% 
#                 mutate(next_transplant_date = lead(transplant_date, order_by = transplant_date), 
#                        next_transplant_occurrence_id = lead(transplant_occurrence_id, order_by = transplant_date),
#                         transplant_count = n_distinct(transplant_date)) %>%
#                 filter(transplant_count == 2, next_transplant_date - transplant_date > 2) %>%
#                 select(person_id, transplant_date, next_transplant_date, transplant_occurrence_id, next_transplant_occurrence_id)

#      two_transplants_patid <- two_transplants_patid %>% select(person_id, transplant_date, transplant_occurrence_id) %>% 
#                 union(two_transplants_patid %>% select(person_id, transplant_date = next_transplant_date, transplant_occurrence_id = next_transplant_occurrence_id))

#       single_transplants_patid %>% union(single_transplants_dup_patid) %>% 
#                 union(two_transplants_patid) %>%
#                 output_tbl("transplant_patid", local = TRUE, file = TRUE)
#       # now onto patients with more than 2 transplants
#       # n = 8
#       # 4180183 probably had 3, the last 2 are 1 day apart, we should take the 1st one as ce_date and the 3rd one as relapse date
#       # let's just discard these patients
#       t %>% filter(transplant_count > 2) %>% select(person_id) %>% 
#                 output_tbl("multi_transplants_patid", local = TRUE, file = TRUE)

      rslt$death <- cdm_tbl("death") %>% inner_join(results_tbl("study_cohorts"), by = "person_id") %>%
            select(person_id, death_date, cause_concept_name, death_impute_concept_name) %>%
            compute_new()
      rslt$death %>% collect_new() %>% output_tbl("covar_death")
      append_sum(cohort = 'death',
                 persons = rslt$death %>% distinct_ct())

   
  # Write step summary log to CSV and/or database,
  # as determined by configuration
  output_sum(name = "outcome_covariates_attrition", local = TRUE, file = TRUE)  

  message('Done.')

  invisible(rslt)

}
