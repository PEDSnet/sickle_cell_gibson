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

  message('get all patients ids from study cohorts')
  message('add in patients whose transplant status is unknown')

  rslt$cohort <- results_tbl("study_cohorts") %>% select(person_id, site) #%>%
                    #  union(results_tbl("multi_transplant_cohort") %>% select(person_id, site)) %>%
                    # distinct(person_id, site)
  rslt$personal_tbls <- personal_tables() %>% filter(output_name %in% c('adt_occurrence',
                                                                        'condition_occurrence',
                                                                        'drug_exposure',
                                                                        'measurement',
                                                                        'measurement_labs',
                                                                        'measurement_organism', 
                                                                        'measurement_vitals', 
                                                                        'visit_occurrence',
                                                                        'person',
                                                                        'procedure_occurrence',
                                                                        'specimen', 
                                                                        'death'))

  rslt$impersonal_tbls <- impersonal_tables() %>% filter(output_name %in% c('condition_occurrence_concept_name',
                                                                            'drug_exposure_concept_name',
                                                                            'measurement_concept_name',
                                                                            'procedure_occurrence_concept_name'))

  rslt$mini_cdm <- build_mini_cdm(cohort=rslt$cohort,
                                    clean=FALSE,
                                    personal= NULL,  #rslt$personal_tbls,
                                    impersonal=rslt$impersonal_tbls,
                                    materialize=TRUE)
  # rslt$mini_cdm$cdm_person %>% compute_new() %>% output_tbl("cdm_person")
  # results_tbl(in_schema("test_gibson", "cdm_person")) 
  # Write step summary log to CSV and/or database,
  # as determined by configuration
  
  output_sum()

  message('Done.')

  invisible(rslt)

}
