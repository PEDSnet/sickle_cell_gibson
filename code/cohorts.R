find_scd_conditions <- function(){
  # patients with relevant non-malignant blood disorders include: sickle cell disease, beta thalassemia major, aplastic anemia, Diamond-Blackfan anemia (congenital pure red cell aplasia)
  scd_id <- load_codeset("scd_dx") %>% mutate(type = "scd") %>% #sickle cell disease
            union(load_codeset("beta_thal_major_dx") %>% mutate(type = "bta")) %>% #beta thalassemia major
            union(load_codeset("aplastic_anemia_dx") %>% mutate(type = "aa")) %>% #aplastic anemia
            union(load_codeset("diamond_blackfan_anemia_dx") %>% mutate(type = "dba")) %>% ##diamond blackfan anemia
            compute_new(temp = TRUE, name = "scd_id")
  tbl <- cdm_tbl("condition_occurrence") %>% inner_join(scd_id, by = c("condition_concept_id" = "concept_id")) %>%
                select(scd_occurrence_id=condition_occurrence_id, 
                        person_id, 
                        scd_concept_id = condition_concept_id, 
                        scd_concept_name= condition_concept_name, 
                        scd_start_date = condition_start_date, 
                        scd_type = type, 
                        scd_start_age_in_months = condition_start_age_in_months, site) 
  return(tbl)
}

# could be used for mri, phleb
find_procedures <- function(procedure_codeset_name){
    procedure_id <- load_codeset(procedure_codeset_name) %>% compute_new()
    tbl <- cdm_tbl("procedure_occurrence") %>% inner_join(procedure_id, by = c("procedure_concept_id" = "concept_id")) %>%
                select(procedure_occurrence_id, 
                        person_id, 
                        procedure_concept_id, 
                        procedure_concept_name, 
                        procedure_date, 
                        visit_occurrence_id, 
                        site)
    return(tbl)
}

find_transplants <- function(procedure_codeset_name){
    procedure_id <- load_codeset(procedure_codeset_name) %>% compute_new()
    tbl <- cdm_tbl("procedure_occurrence") %>% inner_join(procedure_id, by = c("procedure_concept_id" = "concept_id")) %>%
                select(procedure_occurrence_id, 
                        person_id, 
                        procedure_concept_id, 
                        procedure_concept_name, 
                        procedure_date, 
                        visit_occurrence_id, 
                        site)
    
    # get visit start and end dates
    tbl <- tbl %>% inner_join(cdm_tbl("visit_occurrence") %>% select(person_id,
                                                            visit_occurrence_id, 
                                                            visit_start_date, 
                                                            visit_end_date), by = c("visit_occurrence_id", 
                                                                                  "person_id")) %>% 
                                                        filter(procedure_date >= visit_start_date,
                                                               procedure_date <= visit_end_date) %>% collect_new()

    tbl <- tbl %>% rename(transplant_concept_id = procedure_concept_id, 
                          transplant_concept_name = procedure_concept_name, 
                          transplant_date = procedure_date, 
                          transplant_occurrence_id = procedure_occurrence_id,
                          transplant_site = site) %>%
            mutate(transplant_type = case_when(grepl("autologous",transplant_concept_name) ~ "autologous",
                                                grepl("allogeneic",transplant_concept_name) ~ "allogeneic",
                                                TRUE ~ "other")) 
                                                                                 
    # find all patients with at least 2 transplants, filter out case when procedures outside visit windows 
    # group by visit, if the difference between 2 transplants is 1 day, then we'll take the 2nd one as the ce_date assuming that the first one
    # was cancelled
    t <- tbl %>% #mutate(transplant_date = as.Date(transplant_date, format = "%Y-%m-%d")) %>%
                group_by(person_id) %>%
                summarise(transplant_count = n_distinct(transplant_date))

    single_transplants_patid <- t %>% filter(transplant_count == 1) %>% 
                inner_join(tbl, by = c("person_id")) %>% 
                distinct(person_id, transplant_date, transplant_occurrence_id)
                 
    # when the gap between 2 transplants is a day, i think we could say that this is the same transplant and the transplant date is the later one
    # n = 9
    single_transplants_dup_patid <- t %>% filter(transplant_count == 2) %>% 
                inner_join(tbl, by = c("person_id")) %>% 
                group_by(person_id) %>% 
                mutate(next_transplant_date = lead(transplant_date, order_by = transplant_date), 
                       next_transplant_occurrence_id = lead(transplant_occurrence_id, order_by = transplant_date),
                        transplant_count = n_distinct(transplant_date)) %>%
                filter(next_transplant_date - transplant_date <= 2) %>%
                select(person_id, transplant_date = next_transplant_date, transplant_occurrence_id = next_transplant_occurrence_id)

    # when the gap between 2 transplants is more than a day but there were only 2 transplants, it's ok and we'll take the 1st transqplant date
     # as ce_date and the 2nd one as the relapse date, 
     # n = 38 
     # n = 10 cases when both transplants happened in the same visit, already checked, not duplicated transplants
     two_transplants_patid <- t %>% filter(transplant_count == 2) %>% 
                inner_join(tbl, by = c("person_id")) %>% 
                group_by(person_id) %>% 
                mutate(next_transplant_date = lead(transplant_date, order_by = transplant_date), 
                       next_transplant_occurrence_id = lead(transplant_occurrence_id, order_by = transplant_date),
                        transplant_count = n_distinct(transplant_date)) %>%
                filter(next_transplant_date - transplant_date > 2) %>%
                select(person_id, transplant_date, next_transplant_date, transplant_occurrence_id, next_transplant_occurrence_id)

     two_transplants_patid <- two_transplants_patid %>% select(person_id, transplant_date, transplant_occurrence_id) %>% 
                union(two_transplants_patid %>% select(person_id, transplant_date = next_transplant_date, transplant_occurrence_id = next_transplant_occurrence_id))

    valid_patid <- single_transplants_patid %>% union(single_transplants_dup_patid) %>% 
                union(two_transplants_patid) 
      # now onto patients with more than 2 transplants
      # n = 8
      # 4180183 probably had 3, the last 2 are 1 day apart, we should take the 1st one as ce_date and the 3rd one as relapse date
      # let's just discard these patients
    invalid_patid <- t %>% filter(transplant_count > 2) %>% select(person_id)

    tbl <- tbl %>% inner_join(valid_patid, by = c("person_id", "transplant_occurrence_id", "transplant_date")) %>% 
                anti_join(invalid_patid, by = c("person_id"))

    return(tbl)
}

find_phleb <- function(){
    phleb_id <- load_codeset("therapeutic_phleb_px") %>% compute_new(temp = TRUE, name = "phleb_id")
    tbl <- cdm_tbl("procedure_occurrence") %>% inner_join(phleb_id, by = c("procedure_concept_id" = "concept_id")) %>%
                select(phleb_occurrence_id= procedure_occurrence_id, 
                        person_id, 
                        phleb_concept_id = procedure_concept_id, 
                        phleb_concept_name = procedure_concept_name, 
                        phleb_date = procedure_date, 
                        phleb_site = site)
    return(tbl)
}

find_drugs <- function(dx_codeset){
    tbl <- cdm_tbl("drug_exposure") %>% inner_join(dx_codeset, by = c("drug_concept_id" = "concept_id")) %>%
                select(drug_exposure_id, 
                        person_id, 
                        drug_concept_id, 
                        drug_concept_name, 
                        drug_exposure_start_date, 
                        frequency,
                        dose_unit_concept_name,
                        effective_drug_dose,
                        quantity,
                        days_supply,
                        drug_type = type)
    return(tbl)
}

find_measurements <- function(mx_codeset, cohort = NULL){
    tbl <- cdm_tbl("measurement") %>% inner_join(mx_codeset, by = c("measurement_concept_id" = "concept_id")) %>%
                select(measurement_id, 
                        person_id, 
                        ferritin_concept_id = measurement_concept_id, 
                        ferritin_concept_name = measurement_concept_name, 
                        ferritin_date = measurement_date, 
                        operator_concept_name, 
                        value_as_number, 
                        unit_concept_name,
                        range_high,
                        range_high_operator_concept_name,
                        range_low, 
                        range_low_operator_concept_name,
                        unit_concept_name)
    return(tbl)
}

# find patients with transplant between start_date and end_date
# if multiple transplants, only keep the earliest one
# for pre-transplant measurements, get the latest one within 1 year before transplant
# for post-transplant measurements, get the earliest one within 1 year after transplant

find_scd_transplant_with_procedure <- function(scd_cohort,
                                                transplant_cohort,
                                                procedure_cohort,
                                                procedure_date,
                                                start_date = as.Date("2010-01-01"),
                                                end_date,
                                                is_pre = FALSE,
                                                cut_off = 365, min_days = 0){

    tbl <- transplant_cohort %>% filter(transplant_date >= start_date, transplant_date <= end_date) %>%
                group_by(person_id) %>%
                slice_min(transplant_date, n = 1, with_ties = FALSE) %>% ungroup() %>%
                inner_join(scd_cohort, by = "person_id") %>% collect_new() %>%
                # filter(scd_start_date <= transplant_date) %>% #condition disnosed before transplant
                group_by(person_id) %>%
                slice_min(abs(as.numeric(difftime(scd_start_date, transplant_date, units = "days"))), n = 1, with_ties = FALSE) %>% ungroup() %>%
                inner_join(procedure_cohort %>% collect_new(), by = c("person_id")) 

    # find scd patients who had a transplant between start_date and end_date
    if(is_pre){
        tbl <- tbl %>% filter({{procedure_date}} <= transplant_date)
    }else{
        tbl <- tbl %>% filter({{procedure_date}} >= transplant_date)
    }
    
    tbl <- tbl %>% filter(abs(as.numeric(difftime({{procedure_date}}, transplant_date, units = "days"))) >= min_days) %>% 
                group_by(person_id) %>%
                slice_min(abs(as.numeric(difftime({{procedure_date}}, transplant_date, units = "days"))), n = 1, with_ties = FALSE) %>% 
                ungroup() %>% collect()                      

    # get demographics
    tbl <- tbl %>% copy_to_new(df = ., name = "sdsd") %>%
                                inner_join(cdm_tbl("person") %>% 
                                mutate(gender = case_when(gender_concept_id == 8532 ~ "F",
                                                          gender_concept_id == 8507 ~ "M",
                                                          TRUE ~ NA)) %>%
                                select(person_id, birth_date, gender), by = c("person_id"), copy = TRUE) %>%
                    left_join(cdm_tbl("death") %>% select(person_id, death_date, cause_concept_id), by = c("person_id"), copy= TRUE) %>% 
                    mutate(age_at_transplant = (transplant_date - birth_date)/365.25)                                
    return(tbl)
}

# get outcome measurements within 365 days following ce_date
get_outcome_measurements <- function(cohort,
                                     lab_tbl=cdm_tbl("measurement"),
                                     mx_codeset = load_codeset(codeset_name),
                                     min_days = 0,
                                     max_days = 365,
                                     is_pre = FALSE,
                                     is_closet = TRUE) {
  tbl <- cohort %>%
    inner_join(lab_tbl %>% select(measurement_date, 
                                    person_id,
                                    measurement_concept_id,
                                    measurement_concept_name, 
                                    # measurement_type_concept_id,
                                    # operator_concept_id, 
                                    # operator_concept_name,
                                    range_high,
                                    # range_high_operator_concept_id,
                                    range_high_operator_concept_name,
                                    range_low,
                                    # range_low_operator_concept_id,
                                    range_low_operator_concept_name, 
                                    # specimen_concept_id,
                                    # unit_concept_id, 
                                    unit_concept_name, 
                                    value_as_number,
                                    # value_as_concept_id,
                                    value_as_concept_name,
                                    # measurement_order_date, 
                                    # measurement_result_date
                                    ), by = 'person_id') %>% collect_new()
  if(is_pre){
    tbl <- tbl %>% filter(as.numeric(difftime(ce_date, measurement_date, units = "days")) >= min_days,
                          as.numeric(difftime(ce_date, measurement_date, units = "days")) <= max_days)
  } else {
    tbl <- tbl %>% filter(as.numeric(difftime(measurement_date, ce_date, units = "days")) >= min_days,
                          as.numeric(difftime(measurement_date, ce_date, units = "days")) <= max_days)
  }
  if(is_closet){
    tbl <- tbl %>% group_by(person_id, transplant_date, measurement_concept_id) %>% 
                  slice_min(abs(difftime(measurement_date, ce_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup()
  }
    tbl <- tbl %>% inner_join(mx_codeset %>% select(concept_id) %>% collect_new(), 
                              by = c("measurement_concept_id" = "concept_id"))   
    return(tbl)
}

get_outcome_drugs <- function(cohort,
                              dx_codeset,
                              min_days = 0,
                              max_days = 365,
                              is_pre = FALSE,
                              is_closet = TRUE) {
  tbl <- cohort %>% inner_join(cdm_tbl("drug_exposure") %>% select(drug_exposure_id, 
                                                                person_id, 
                                                                drug_concept_id, 
                                                                drug_concept_name, 
                                                                drug_exposure_start_date, 
                                                                frequency,
                                                                dose_unit_concept_name,
                                                                effective_drug_dose,
                                                                quantity,
                                                                days_supply, 
                                                                route_concept_name), by = c("person_id")) %>% collect_new()
  if(is_pre){
    tbl <- tbl %>% filter(as.numeric(difftime(ce_date, drug_exposure_start_date, units = "days")) >= min_days,
                          as.numeric(difftime(ce_date, drug_exposure_start_date, units = "days")) <= max_days)
  } else {
    tbl <- tbl %>% filter(as.numeric(difftime(drug_exposure_start_date, ce_date, units = "days")) >= min_days,
                          as.numeric(difftime(drug_exposure_start_date, ce_date, units = "days")) <= max_days)
  }
  if(is_closet){
    tbl <- tbl %>% group_by(person_id, transplant_date, drug_concept_id) %>% 
                  slice_min(abs(difftime(drug_exposure_start_date, ce_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup()
  }
    tbl <- tbl %>% inner_join(dx_codeset %>% select(concept_id) %>% collect_new(), by = c("drug_concept_id" = "concept_id")) 
    return(tbl)
}

get_height_z_scores <- function(cohort,
                                vital_tbl = cdm_tbl("measurement_vitals"),
                                demographic_tbl = cdm_tbl("person")){
  
  # cohort_dem <- demographic_tbl %>% 
  #   inner_join(cohort, by = "patid") %>% 
  #   distinct(birth_date, patid, sex, ce_date) %>% 
  #   compute_new()
  
  # abs ht in inches
  vital <- vital_tbl %>% 
    inner_join(cohort, by = "person_id") %>% 
    compute_new()
  
  heights <- vital %>%
    filter(!is.na(ht)) %>%
    distinct(patid, ce_date, sex, birth_date, ht, height_date = measure_date) %>%
    mutate(age = (height_date - birth_date) / 365.25,
           age_ce = (ce_date - birth_date) / 365.25) %>%
    distinct(patid, ht, age, height_date, sex) %>%
    collect_new()
  
  # to get z-score, need age, sex
  heights$height_z <-
    with(heights, ifelse(
      age >= 2,
      sds(
        ht * 2.54,
        age,
        sex,
        ref = cdc.ref,
        item = "height2_20",
        type = "SDS",
        male = "M",
        female = "F"
      ),
      NA
    ))
  
  # BIV from https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm
  heights$height_z <-
    with(heights, ifelse(height_z >= -5 &
                           height_z <= 4, height_z, NA))
  
  heights <- heights %>%
    filter(!is.na(height_z))
  
}

get_outcome_conditions <- function(cohort,
                                   dx_codeset,
                                   min_days = 0,
                                   max_days = 365,
                                   is_pre = FALSE,
                                   is_closet = TRUE) {
  tbl <- cohort %>% inner_join(cdm_tbl("condition_occurrence") %>% select(condition_occurrence_id, 
                                                                person_id, 
                                                                condition_concept_id, 
                                                                condition_concept_name, 
                                                                condition_start_date), by = c("person_id")) %>% collect_new()
  if(is_pre){
    tbl <- tbl %>% filter(as.numeric(difftime(ce_date, condition_start_date, units = "days")) >= min_days,
                          as.numeric(difftime(ce_date, condition_start_date, units = "days")) <= max_days)
  } else {
    tbl <- tbl %>% filter(as.numeric(difftime(condition_start_date, ce_date, units = "days")) >= min_days,
                          as.numeric(difftime(condition_start_date, ce_date, units = "days")) <= max_days)
  }
  if(is_closet){
    tbl <- tbl %>% group_by(person_id, transplant_date, condition_concept_id) %>%
                  slice_min(abs(difftime(condition_start_date, ce_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup()
  }
    tbl <- tbl %>% inner_join(dx_codeset %>% select(concept_id) %>% collect_new(), by = c("condition_concept_id" = "concept_id"))
    return(tbl)
}

get_outcome_procedures <- function(cohort,
                                   px_codeset,
                                   min_days = 0,
                                   max_days = 365, 
                                   is_pre = FALSE,
                                   is_closet = TRUE) {
  tbl <- cohort %>% inner_join(cdm_tbl("procedure_occurrence") %>% select(procedure_occurrence_id, 
                                                                person_id, 
                                                                procedure_concept_id, 
                                                                procedure_concept_name, 
                                                                procedure_date), by = c("person_id")) %>% collect_new()
  if(is_pre){
    tbl <- tbl %>% filter(as.numeric(difftime(ce_date, procedure_date, units = "days")) >= min_days,
                          as.numeric(difftime(ce_date, procedure_date, units = "days")) <= max_days)
  } else {  
    tbl <- tbl %>% filter(as.numeric(difftime(procedure_date, ce_date, units = "days")) >= min_days,
                          as.numeric(difftime(procedure_date, ce_date, units = "days")) <= max_days)
  }
  if(is_closet){
    tbl <- tbl %>% group_by(person_id, transplant_date, procedure_concept_id) %>% 
                  slice_min(abs(difftime(procedure_date, ce_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup()
  }
    tbl <- tbl %>% inner_join(px_codeset %>% select(concept_id) %>% collect_new(), by = c("procedure_concept_id" = "concept_id"))                
    return(tbl)
}
