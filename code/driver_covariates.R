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
  library(ggplot2)
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
  
  # exclude the ones with percent as unit
  rslt$CD348 %>% mutate(CDtype = case_when(grepl("T4 helper", measurement_concept_name, ignore.case = TRUE) ~ "CD4",
                                grepl("CD4", measurement_concept_name, ignore.case = TRUE) ~ "CD4",
                                grepl("T8 suppressor", measurement_concept_name, ignore.case = TRUE) ~ "CD8",
                                grepl("CD8", measurement_concept_name, ignore.case = TRUE) ~ "CD8",
                                grepl("CD3 cells", measurement_concept_name, ignore.case = TRUE) ~ "CD3",
                                grepl("CD3", measurement_concept_name, ignore.case = TRUE) ~ "CD3",
                                TRUE ~ NA)) %>%
                                filter(!is.na(value_as_number), unit_concept_name != "percent") %>%
                filter(!((measurement_concept_id == "3011412" & unit_concept_name == "No matching concept" & range_high < 100) |
                        (measurement_concept_id == "3028167" & unit_concept_name == "No matching concept" & is.na(range_high) & is.na(range_low))|
                        (measurement_concept_id == "3045450" & unit_concept_name == "No matching concept" & is.na(range_high) & is.na(range_low))|
                        (measurement_concept_id == "3014037" & unit_concept_name == "No matching concept" & range_high <100) |
                        (measurement_concept_id == "3045450" & unit_concept_name == "No information"  & is.na(range_high) & is.na(range_low)))) %>%
                collect() %>% output_tbl(name = "covar_CD348_mx")
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
                                filter(!((unit_concept_name == "NA" & is.na(range_low)) | 
                                        (measurement_concept_id == "3028026"& unit_concept_name == "No matching concept" & is.na(range_low))|
                                        (unit_concept_name == "No information" & is.na(range_low)))) %>% collect() %>% 
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
   rslt$ferritin <- rbind(rslt$ferritin_pre %>% mutate(ferritin_type = "pre"), 
                        rslt$ferritin_post %>% mutate(ferritin_type = "post")) %>% 
                        pivot_longer(cols = c(aim_2a_2:aim_3_2), names_to = "aim", values_to = "which_aim") %>%
                        filter(which_aim == TRUE) %>% 
                        select(person_id, 
                        transplant_date,
                        aim,
                        ferritin_date = measurement_date,
                        ferritin = value_as_number,
                        ferritin_type, 
                        measurement_concept_id, measurement_concept_name, unit_concept_name, range_high, range_low) %>% 
                        mutate(ferritin_level = case_when(ferritin < 1000 ~ "low",
                                                ferritin >= 1000 & ferritin <= 3000 ~ "moderate",
                                                ferritin > 3000 ~ "high",
                                                TRUE ~ NA)) 
   # majority with units nanogram per milliliter
   # assume missing units are nanogram per milliliter  
   # no unit conversion done                                              
   rslt$ferritin %>% collect() %>% output_tbl(name = "covar_ferritin_mx")                                                            
   append_sum(cohort = 'covar_ferritin_px',
                     persons = rslt$ferritin %>% distinct_ct())

   # the previous ferritin table was indexed based on the transplant dates, we just need all the ferritin measurements 
   rslt$ferritin <- find_measurements(cohort = results_tbl("study_cohorts"), 
                                mx_codeset = load_codeset("ferritin_lab")) %>% 
                    select(person_id, 
                        ferritin_date = measurement_date,
                        ferritin = value_as_number,
                        measurement_concept_id, measurement_concept_name, unit_concept_name, range_high, range_low) %>% 
                        mutate(ferritin_level = case_when(ferritin < 1000 ~ "low",
                                                ferritin >= 1000 & ferritin <= 3000 ~ "moderate",
                                                ferritin > 3000 ~ "high",
                                                TRUE ~ NA)) 
    rslt$ferritin %>% collect() %>% output_tbl(name = "covar_ferritin_mx")    
   
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
  rslt$defibrotide <- find_drugs(dx_codeset = load_codeset("defibrotide_rx") %>% mutate(type = FALSE)) %>%
                        inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id")
  rslt$abdominal_ultrasound <- find_procedures(procedure_codeset_name = "abdominal_ultrasound_px") %>%
                        inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id")
  rslt$bilirubin <- find_measurements(cohort = results_tbl("study_cohorts"),
                                mx_codeset = load_codeset("bilirubin_mx"))

   #select patients with defibrotide or (abdominal ultrasound and bilirubin > 2mg/L 2 days apart)
   rslt$defibrotide %>% collect() %>% output_tbl(name = "defibrotide_rx")
   rslt$abdominal_ultrasound %>% collect() %>% output_tbl(name = "abdominal_ultrasound_px")
   rslt$bilirubin %>% collect() %>% output_tbl(name = "bilirubin_mx")
  
  # expecting many to many relationship because of multiple ultrasound per day and multiple bilirubin measurements per day
  # patients with (abdominal ultrasound and bilirubin > 2mg/L 3 days apart) within 90 days after transplant
  # dont do the filter dates now because we dont know the transplant dates yet
  bili_and_ultrasound <- results_tbl("bilirubin_mx") %>% filter(value_as_number> 2) %>%
                        select(person_id, bilirubin = value_as_number, 
                                bilirubin_date = measurement_date) %>% 
                                distinct(person_id, bilirubin_date) %>% 
                        inner_join(results_tbl("abdominal_ultrasound_px") %>% 
                                select(person_id, abd_ultrasound_date = procedure_date) %>%
                                distinct(person_id, abd_ultrasound_date), 
                                by = c("person_id")) %>% #, relationship = "many-to-many") %>% 
                        mutate(min_bili_abd_date := pmin(bilirubin_date, abd_ultrasound_date)) %>% collect() %>%
                        filter(abs(as.numeric(difftime(abd_ultrasound_date, bilirubin_date, units = "days"))) <= 2) %>%
                        distinct(person_id, min_bili_abd_date) %>% 
                        mutate(vod_cat_1 = 1)
                        # group_by(person_id, transplant_date) %>% 
                        # slice_min(bilirubin_date, n = 1, with_ties = FALSE) %>%
                        # filter(as.numeric(difftime(transplant_date, min_bili_abd_date, units = "days")) <= 90,
                        #         min_bili_abd_date >= transplant_date)
    bili_and_ultrasound %>% view()  
    
    # patients with (defibrotide AND abdominal ultrasound 2 days apart) within 90 days after transplant
    defi_and_ultrasound <- results_tbl("defibrotide_rx") %>% #filter(as.numeric(difftime(drug_exposure_start_date, transplant_date, units = "days")) <= 90, 
                                                #        drug_exposure_start_date >= transplant_date) %>%
                           distinct(person_id, drug_exposure_start_date) %>%  
                           inner_join(results_tbl("abdominal_ultrasound_px") %>% 
                                     select(person_id, abd_ultrasound_date = procedure_date) %>%
                                     distinct(), by = c("person_id")) %>%
                           mutate(min_defi_abd_date := pmin(drug_exposure_start_date, abd_ultrasound_date)) %>% collect() %>%
                           filter(abs(as.numeric(difftime(abd_ultrasound_date, drug_exposure_start_date, units = "days"))) <= 2) %>%
                           distinct(person_id, min_defi_abd_date) %>% mutate(vod_cat_2 = 2)
    defi_and_ultrasound %>% view()

    # (Abdominal ultrasound AND total bilirubin >2 mg/dL within 2 days AND significant weight gain defined as: 
    # weight gain on 3 consecutive measurements within 10 days of ultrasound OR weight gain >5% of 
    # admission weight within 10 days of ultrasound) within 90 days of transplant     thed host    
    # need the admission date
    rslt$weight <- results_tbl("cdm_measurement") %>% filter(measurement_concept_id == 3013762) %>%
        inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id") %>%
        select(weight_date = measurement_date, measurement_id, person_id, weight = value_as_number)

    bili_and_ultrasound_and_weight_gain <- rslt$weight %>% collect() %>% 
        group_by(person_id, weight_date) %>%
        slice_max(weight, n = 1, with_ties = FALSE) %>% ungroup() %>%
        inner_join(bili_and_ultrasound, by = "person_id") %>%
        distinct(person_id, min_bili_abd_date, weight_date, .keep_all = TRUE) %>% mutate(vod_cat = 3) %>% 
        group_by(person_id) %>%
        arrange(weight_date) %>% 
        mutate(weight_gain_consec = if_else((weight - lag(weight, default = first(weight), n = 1L)) > 0 & (lag(weight, default = first(weight), n = 1L) - lag(weight, default = first(weight), n = 2L)) > 0, 1, 0)) %>% 
        filter(abs(as.numeric(difftime(min_bili_abd_date, weight_date, units = "days"))) <= 10, 
                weight_gain_consec == 1) %>%
        slice_min(abs(as.numeric(difftime(min_bili_abd_date, weight_date, units = "days"))), n = 1, with_ties = FALSE) %>% 
        select(person_id, weight_date, min_bili_abd_date) %>%
        mutate(min_bili_abd_weight_date := pmin(weight_date, min_bili_abd_date)) %>% 
        distinct(person_id, min_bili_abd_weight_date) %>%
        mutate(vod_cat_3 = 3) %>% ungroup() 
        
    bili_and_ultrasound_and_weight_gain %>% view()
    # Peritoneal drain placement procedure during transplant admission
    # (this should show up as a procedure code during admission)  
#     results_tbl("cdm_procedure_occurrence") %>% inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id") %>%
#         filter(grepl("peritoneal",procedure_concept_name, ignore.case = TRUE),
#                 grepl("insertion",procedure_concept_name, ignore.case = TRUE) |
#                 grepl("catheter",procedure_concept_name, ignore.case = TRUE)) %>% 
#                 distinct(procedure_concept_name,  procedure_concept_id) %>% 
#                 view()

    # could not find any codeset for peritoneal drain placement
    # find codeset for abdominal drainage instead
#     results_tbl("cdm_procedure_occurrence") %>% inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id") %>%
#         filter(grepl("abdominal",procedure_concept_name, ignore.case = TRUE),
#                 grepl("drainage",procedure_concept_name, ignore.case = TRUE)) %>% 
#                 distinct(procedure_concept_name, procedure_concept_id) %>%
#                 distinct(procedure_concept_name,  procedure_concept_id)

    # the relevant codes for Peritoneal drain placement procedure
    # and output results to a table
    results_tbl("cdm_procedure_occurrence") %>% inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id") %>%
                filter(procedure_concept_id %in% c(4816373,  2003547, 44816426, 2614742, 2781295, 2109464)) %>%  
                collect() %>% output_tbl(name = "covar_abd_drainage_px")

    # do the date filtering later
    abd_drainage <- results_tbl("covar_abd_drainage_px") %>% 
                inner_join(results_tbl("study_cohorts") %>% distinct(person_id, transplant_date), by = "person_id") %>%
                # filter(abs(as.numeric(difftime(procedure_date, transplant_date, units = "days"))) <= 7) %>% 
                select(person_id, abd_date = procedure_date) %>%
                mutate(vod_cat_4 = 4) 

    # Veno-occlusive disease (VOD) or sinusoidal obstruction syndrome (SOS) on patient problem list --- new criterion        
    results_tbl("cdm_condition_occurrence") %>% inner_join(results_tbl("study_cohorts") %>% distinct(person_id), by = "person_id") %>%
        filter(grepl("veno-occlusive", condition_concept_name, ignore.case = TRUE)) %>% 
        # distinct(condition_concept_name) %>% view()
        output_tbl(name = "covar_venoocclusive_dx")

    vod_dx <- results_tbl("covar_venoocclusive_dx") %>% collect_new() %>%
                inner_join(cohort_covars %>% distinct(person_id, transplant_date), by = "person_id") %>%
                filter(condition_status_concept_id == 4230359 |
                        condition_type_concept_id %in% c(2000000089)) %>%
                # filter(condition_start_date >= transplant_date) %>%
                select(person_id, veno_date = condition_start_date) %>% 
                mutate(vod_cat_5 = 5)

    defi_and_ultrasound %>% view()
    defi_and_ultrasound %>% distinct_ct()

    VOD <- defi_and_ultrasound %>% full_join(bili_and_ultrasound, by = c("person_id")) %>% 
                full_join(bili_and_ultrasound_and_weight_gain, by = c("person_id")) %>%
                full_join(abd_drainage %>% collect(), by = c("person_id")) %>%
                full_join(vod_dx, by = c("person_id")) %>% 
                mutate(vod = TRUE) %>% output_tbl(name = "covar_vod_dx")

    VOD %>% view()
    

  append_sum(cohort = 'defibrotide_drugs',
             persons = rslt$defibrotide %>% distinct_ct())

  rslt$blood_culture <- find_organism(cohort=results_tbl("study_cohorts"),
                                     lab_tbl=cdm_tbl("measurement_organism")) 
  rslt$blood_culture %>% filter(specimen_source_value != "Urine|Unknown (Blood)", 
                                !grepl("yeast", organism_concept_name, ignore.case = TRUE)) %>% 
                                output_tbl(name = "covar_bacteremia_dx")

  rslt$blood_culture %>% view()
  append_sum(cohort = 'blood_culture',
             persons = rslt$blood_culture %>% distinct_ct())
 
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
    
   rslt$gvhd <- find_conditions(cohort = results_tbl("study_cohorts"), 
                                dx_codeset = load_codeset("gvhd_dx")) 
#    rslt$gvhd %>% group_by(condition_concept_id, condition_concept_name) %>% 
#                 summarise(n = n()) %>% view()
   rslt$gvhd %>% collect() %>% output_tbl(name = "covar_gvhd_dx")
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
    dx_conditioning <- load_codeset("conditioning_drugs_rx") %>% 
                        mutate(type = case_when(grepl("busulfan", concept_name, ignore.case = TRUE) ~ "busulfan",
                                                grepl("alemtuzumab", concept_name, ignore.case = TRUE) ~ "alemtuzumab",
                                                grepl("fludarabine", concept_name, ignore.case = TRUE) ~ "fludarabine",
                                                grepl("cyclophosphamide", concept_name, ignore.case = TRUE) ~ "cyclophosphamide",
                                                grepl("hydroxyurea", concept_name, ignore.case = TRUE) ~ "hydroxyurea",
                                                grepl("melphalan", concept_name, ignore.case = TRUE) ~ "melphalan",
                                                grepl("thiotepa", concept_name, ignore.case = TRUE) ~ "thiotepa",
                                                grepl("treosulfan", concept_name, ignore.case = TRUE) ~ "treosulfan",
                                                TRUE ~ "Other")) 
    find_drugs(dx_conditioning) %>% collect() %>% output_tbl("conditioning_drugs_rx")
    find_procedures(procedure_codeset_name = "conditioning_tbi_px") %>% collect() %>% output_tbl("conditioning_tbi_px")

    results_tbl("conditioning_drugs_rx") %>% select(person_id,
                                concept_id = drug_concept_id, 
                                concept_name = drug_concept_name, 
                                exposure_date = drug_exposure_start_date,
                                drug_type) %>% mutate(conditioning_type = "drugs") %>%
                             union(results_tbl("conditioning_tbi_px") %>% select(person_id, 
                                                                concept_id = procedure_concept_id,
                                                                concept_name = procedure_concept_name,
                                                                exposure_date = procedure_date) %>% mutate(conditioning_type = "tbi", drug_type = NA)) %>%
                             collect() %>% output_tbl("conditioning_agents_rx_px")

    # neutrophil counts with respect to transplant date
#     rslt$ANC <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
#                                         mx_codeset = load_codeset("neutrophils_mx"),
#                                         is_pre = FALSE, is_closet = FALSE)
    # all neutrophil counts for patients in the cohort
    rslt$ANC <- find_measurements(cohort = results_tbl("study_cohorts") %>%
                                        select(person_id, transplant_date, site), 
                                mx_codeset = load_codeset("neutrophils_mx"))  
    # check units
    rslt$ANC %>% collect() %>% output_tbl(name = "covar_anc_mx_original")
    rslt$ANC <- results_tbl("covar_anc_mx_original") %>% collect_new()

    # all the records without information on measurement units are ways after transplant, could look at this later
    rslt$ANC %>% filter(unit_concept_name == "No information", !is.na(value_as_number)) %>% 
                        select(person_id, transplant_date, measurement_date, value_as_number, unit_concept_name) %>% 
                        arrange(person_id, measurement_date) %>% view()
    # check valid concepts 
    rslt$ANC %>% distinct(measurement_concept_id, measurement_concept_name, unit_concept_name) %>% view()

    # these are the ones to use
    rslt$ANC %>% filter(measurement_concept_id %in% c("3013650", "3017501", "3017732", "3038972", "3015586", "3046321", "3018199"),
                unit_concept_name != "No information", !is.na(value_as_number)) %>%
                group_by(measurement_concept_id, 
                                    measurement_concept_name, 
                                    unit_concept_name,
                                    range_high, range_low) %>% summarise(n = n()) %>% view()
    
    rslt$ANC <- results_tbl("covar_anc_mx") %>% collect_new()

      message("last in-person visit")
      rslt$visits <- cdm_tbl("visit_occurrence") %>% inner_join(results_tbl("study_cohorts"), by = "person_id") %>%
            select(person_id, visit_occurrence_id, visit_start_date, visit_end_date, visit_concept_id) %>%
            collect()
      
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
