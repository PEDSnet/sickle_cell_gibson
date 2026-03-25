
require(tidyr)
require(knitr)
require(kableExtra)
require(stringr)
require(tibble)
require(ggplot2)
require(table1)
require(countmaskr)
require(tidyverse)

##Summary Functions

#Connect to a table outside the CDM
connect_table <- function(table_name){
  return(tbl(config('db_src'), table_name))
}

#Get Count

get_patient_count <- function(table_one){
  return(table_one %>% select(person_id) %>% distinct() %>% tally())
}

#Make percentage instead of count
make_percentage <- function(data, field = "Patient_Count"){
  data %>% mutate(Percent_Total = round(!!rlang::sym(field)*100/sum(!!rlang::sym(field)),2))
}

#Get by Gender
gender_breakdown <- function(table_one, concept_name = T){
  if(!concept_name){ 
    table_one %>% 
    left_join(vocabulary_tbl('concept') %>% select(concept_id, gender_concept_name = concept_name),
              by = c('gender_concept_id' = 'concept_id'))}
  
  table_one %>% 
    group_by(gender_concept_id, gender_concept_name) %>%
    summarize(Patient_Count = n_distinct(person_id)) %>%
    select(gender_concept_id, gender_concept_name, Patient_Count)
}

#Get by Race
race_breakdown <- function(table_one, concept_name = T){
  if(!concept_name){ 
    table_one %>% 
      left_join(vocabulary_tbl('concept') %>% select(concept_id, race_concept_name = concept_name),
                by = c('race_concept_id' = 'concept_id'))}
  
  table_one %>% 
    group_by(race_concept_id, race_concept_name) %>%
    summarize(Patient_Count = n_distinct(person_id)) %>%
    select(race_concept_id, race_concept_name, Patient_Count)
}

#Get by Ethnicity
ethnicity_breakdown <- function(table_one, concept_name = T){
  if(!concept_name){ 
    table_one %>% 
      left_join(vocabulary_tbl('concept') %>% select(concept_id, ethnicity_concept_name = concept_name),
                by = c('ethnicity_concept_id' = 'concept_id'))}
  
  table_one %>% 
    group_by(ethnicity_concept_id, ethnicity_concept_name) %>%
    summarize(Patient_Count = n_distinct(person_id)) %>%
    select(ethnicity_concept_id, ethnicity_concept_name, Patient_Count)
}

#Get by Site 
site_breakdown <- function(table_one){
  if('site' %in% colnames(table_one)){
   result <- table_one %>% 
      group_by(site) %>%
      summarize(Patient_Count = n_distinct(person_id)) %>%
      select(site, Patient_Count)
  }
  else{
  result <- table_one %>% 
    inner_join(cdm_tbl('person') %>% select(person_id, site), by = 'person_id') %>%
    group_by(site) %>%
    summarize(Patient_Count = n_distinct(person_id)) %>%
    select(site, Patient_Count)
  }
  return(result)
}

# conditions counts
concept_counts <- function(cohort, concept_name, concept_id){
  cohort_ct <- cohort %>% distinct_ct()
  if(("type" %in% colnames(cohort)) && ("site" %in% colnames(cohort))){ 
    tbl <- cohort %>% group_by(site, type, {{concept_id}}, {{concept_name}}) %>%
      summarize(persons = n_distinct(person_id), perc = round(n_distinct(person_id)/cohort_ct*100,2)) 
  } else if(("type" %in% colnames(cohort))){
    tbl <- cohort %>% group_by(type, {{concept_id}}, {{concept_name}}) %>%
      summarize(persons = n_distinct(person_id), perc = round(n_distinct(person_id)/cohort_ct*100,2)) 
  } else if(("site" %in% colnames(cohort))){
    tbl <- cohort %>% group_by(site, {{concept_id}}, {{concept_name}}) %>%
      summarize(persons = n_distinct(person_id), perc = round(n_distinct(person_id)/cohort_ct*100,2)) 
  } else{
    tbl <- cohort %>% group_by({{concept_id}}, {{concept_name}}) %>%
      summarize(persons = n_distinct(person_id), perc = round(n_distinct(person_id)/cohort_ct*100,2)) 
  }
  
  return(tbl)
}

get_codesets <- function(codeset_name) {
  read_csv(paste0('../specs/', codeset_name, '.csv'), col_types = c('c', 'c', 'c', 'c'), show_col_types = FALSE) %>%
  mutate(concept_id = as.character((concept_id)), concept_code = as.character((concept_code))) %>% collect()
}

get_results <- function(tbl_name) {
  if (data_source == 'local') {
    rslt <- read.csv(paste0('../local/', tbl_name, '.csv'))
  }
  else {
    rslt <- results_tbl(tbl_name) %>% collect()
  }
  rslt
}

prettify_kable <- function(data) {
  rslt <- data %>% kable(digits = 4, format.args = list(big.mark = ','))
  if (knitr::is_html_output()) rslt <- rslt %>%
      kable_styling(bootstrap_options = c("striped", "hover")) %>%
      column_spec(1, bold = T, border_right = T)
  rslt
}


