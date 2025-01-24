.run  <- function() {

    library(tidyr)
    library(lubridate)

    library(pdftools)
    library(ggplot2)

    rslt = list()

    message("get transplant dates from chart review and exclude false positive patients")
    # remove false positives and get transplant dates from chart review 
    rslt$study_cohorts <- results_tbl("study_cohorts") %>% 
                    filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                    left_join(results_tbl("cr_data") %>% select(record_id, transplant_date_cr = transplant_date, 
                                                            transplant_type_cr = transplant_type, 
                                                            second_transplant_date_cr = second_transplant_date,
                                                            disease_relapse = graft_fail, chart_completion, eligibility) %>% distinct(), by = "record_id") %>% 
                    mutate(transplant_date = if_else(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                            transplant_type = if_else(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                            chart_completion = if_else(is.na(chart_completion), FALSE, chart_completion)) %>%
                    filter((eligibility == "Yes" & chart_completion) | !chart_completion) %>%
                    select(-transplant_date_cr, -transplant_type_cr)
    
    rslt$study_cohorts %>% view()
    message("Compute platelet engraftment dates")
    rslt$platelet <- find_measurements(cohort = rslt$study_cohorts, 
                                             mx_codeset = load_codeset("platelet_mx")) 
    rslt$platelet <- rslt$platelet %>% filter(aim_2a_1 | aim_2a_2) %>% collect()
    rslt$platelet %>% output_tbl(name = "platelet_mx_org")                                             
    rslt$platelet <- rslt$platelet %>% #filter(aim_2a_1 | aim_2a_2) %>% collect() %>% 
                     # if we have a lot of patients with missing platelet counts, we might need to check the filtered data
                    # filter(measurement_concept_id %in% c( "3007461", "3024929", "4267147"),
                    filter(!is.na(value_as_number)) %>% 
                    mutate(conversion_factor = case_when(grepl("billion per liter", unit_concept_name, ignore.case = TRUE) ~ 1000,
                                                        grepl("Kelvin per cubic millimeter", unit_concept_name, ignore.case = TRUE) ~ 1000,
                                                        grepl("thousand per cubic millimeter", unit_concept_name, ignore.case = TRUE) ~ 1000,
                                                        grepl("Kelvin per microliter", unit_concept_name, ignore.case = TRUE) ~ 1000,
                                                        grepl("thousand per microliter", unit_concept_name, ignore.case = TRUE) ~ 1000,
                                                        grepl("per cubic millimeter", unit_concept_name, ignore.case = TRUE) ~ 1,
                                                        grepl("femtoliter", unit_concept_name, ignore.case = TRUE) ~ -99,
                                                        grepl("percent", unit_concept_name, ignore.case = TRUE) ~ -99,
                                                        grepl("Meter per microliter", unit_concept_name, ignore.case = TRUE) ~ -99,
                                                        grepl("Giant platelets", measurement_concept_name, ignore.case = TRUE) ~ -99,
                                                        grepl("Platelet Ab", measurement_concept_name, ignore.case = TRUE) ~ -99, 
                                                        grepl("No information", unit_concept_name, ignore.case = TRUE) & value_as_number >= 1000 ~ 1, 
                                                        grepl("No information", unit_concept_name, ignore.case = TRUE) & value_as_number < 1000 ~ 1000,
                                                        grepl("No matching concept", unit_concept_name, ignore.case = TRUE) & value_as_number >= 1000 ~ 1, 
                                                        grepl("No matching concept", unit_concept_name, ignore.case = TRUE) & value_as_number < 1000 ~ 1000,
                                                        is.na(unit_concept_name) & value_as_number >= 1000 ~ 1, 
                                                        is.na(unit_concept_name) & value_as_number < 1000 ~ 1000,
                                                        TRUE ~ NA),
                            value_as_number = conversion_factor * value_as_number) %>%
                    filter(!is.na(conversion_factor), value_as_number >= 0)
    rslt$platelet %>% output_tbl(name = "platelet_mx")

    # check units
    # thousand per microliter = per cubic millimeter = Kelvin per microliter (this K actually means thousands)
    # the cut-off is 20K/microliter
    # for measurement_concept_id = 3007461 with different ranges, it seems like we could ignore them because 
    # some of the high ranges are not relevant. e.g no values between 300 and 750
    # exclude the entitic volume as it measures the volume rather than the counts

    rslt$platelet %>% group_by(measurement_concept_id, 
                    measurement_concept_name, 
                    unit_concept_name,
                    range_high, range_low) %>% summarise(n = n()) %>% view()   

    # platelet transfusion procedure 
    # get all the procedures within 30 days before transplant and 150 days after transplant
    # filter for transfusion and get the names
    rslt$platelet_transfusion <- rslt$study_cohorts %>% inner_join(cdm_tbl("procedure_occurrence") %>% select(procedure_occurrence_id, 
                                                                person_id, 
                                                                procedure_concept_id, 
                                                                procedure_concept_name, 
                                                                procedure_date), by = c("person_id")) %>% 
                                filter(grepl("transfusion", procedure_concept_name, ignore.case = TRUE),
                                        grepl("blood", procedure_concept_name, ignore.case = TRUE) | grepl("platelet", procedure_concept_name, ignore.case = TRUE),
                                        !grepl("stem cell", procedure_concept_name, ignore.case = TRUE),
                                        !grepl("blood bank", procedure_concept_name, ignore.case = TRUE)) %>%
                                collect_new() %>%
                                filter((as.numeric(difftime(procedure_date, transplant_date, units = "days")) <= 180 &
                                        as.numeric(difftime(procedure_date, transplant_date, units = "days")) >= 0) |
                                        (as.numeric(difftime(transplant_date, procedure_date, units = "days")) <= 30 &
                                        as.numeric(difftime(transplant_date, procedure_date, units = "days")) >= 0))
    rslt$platelet_transfusion %>% group_by(procedure_concept_name) %>% summarise(n = n_distinct(person_id)) %>% view()
    rslt$platelet_transfusion %>% output_tbl(name = "covar_platelet_transfusion_px")

    # could do either max or mean per day
    rslt$platelet <- results_tbl("platelet_mx") %>% collect_new() %>% 
                filter(aim_2a_1 | aim_2a_2) %>% 
                group_by(person_id, record_id, transplant_date, measurement_date) %>%
                summarise(value_as_number_max = max(value_as_number, na.rm = TRUE),
                        value_as_number_mean = mean(value_as_number, na.rm = TRUE)) %>%
                mutate(platelet_days_since_ce = as.numeric(difftime(measurement_date, transplant_date, units = "days")), 
                        val_ct_max = ifelse(value_as_number_max >= 20000, 1, 0),
                        val_ct_mean = ifelse(value_as_number_mean >= 20000, 1, 0)) %>% 
                ungroup() %>% group_by(person_id, record_id, transplant_date) %>%   
                mutate(platelet_consec_max = ifelse((lead(platelet_days_since_ce, 1) - platelet_days_since_ce == 1) & 
                                        (lead(platelet_days_since_ce, n =2) - platelet_days_since_ce == 2) &
                                        (val_ct_max + lead(val_ct_max, 1) + lead(val_ct_max, n =2)) == 3, TRUE, FALSE),
                        platelet_consec_mean = ifelse((lead(platelet_days_since_ce, 1) - platelet_days_since_ce == 1) & 
                                        (lead(platelet_days_since_ce, n =2) - platelet_days_since_ce == 2) &
                                        (val_ct_mean + lead(val_ct_mean, 1) + lead(val_ct_mean, n =2)) == 3, TRUE, FALSE)) %>% 
                ungroup() %>% 
                # filter(platelet_consec == TRUE) %>%
                select(person_id, record_id, transplant_date, platelet_measurement_date = measurement_date, #ce_date, birth_date, gender, 
                platelet_max_ct = value_as_number_max, platelet_mean_ct = value_as_number_mean, 
                platelet_consec_max, platelet_consec_mean, platelet_days_since_ce)

    rslt$transfusion <- results_tbl("covar_platelet_transfusion_px") %>% collect() %>%
                mutate(transfusion_days_since_ce = as.numeric(difftime(procedure_date, transplant_date, units = "days"))) %>%
                distinct(person_id, record_id, transplant_date, procedure_date, .keep_all= TRUE) %>%
                select(person_id, record_id, transplant_date, transfusion_date = procedure_date, transfusion_days_since_ce)
    rslt$transfusion %>% view()

    rslt$platelet <- rslt$platelet %>% full_join(rslt$transfusion %>% select(-transplant_date), 
                                                by = c("person_id", "platelet_days_since_ce" = "transfusion_days_since_ce", "record_id")) %>%
                    mutate( transfusion_days_since_ce = as.numeric(difftime(transfusion_date, transplant_date, units = "days")))

    # get 2nd transplant dates, for patients without chart review data, use EHR data
    # all chart reviews were done
    rslt$transplant_px <- rslt$study_cohorts %>% 
                        # filter(eligibility == "yes") %>%
                        filter(aim_2a_1 | aim_2a_2) %>% 
                        left_join(results_tbl("no_multi_transplant_px") %>% 
                        select(person_id, second_transplant_date = transplant_date), by = "person_id") %>% collect_new() %>%
                        filter(second_transplant_date >= transplant_date) %>% 
                        mutate(second_transplant_date = if_else(second_transplant_date == transplant_date, NA, second_transplant_date)) %>%
                        group_by(person_id) %>%
                        slice_min(second_transplant_date, n = 1, with_ties = FALSE, na_rm = FALSE) %>% 
                        ungroup()

    # use second transplant dates from chart review 
    rslt$transplant_px <- rslt$transplant_px %>% 
                        mutate(second_transplant_date = case_when(!is.na(second_transplant_date_cr) ~ second_transplant_date_cr, 
                                                                (is.na(second_transplant_date_cr) & chart_completion) ~ NA, #overwrite second transplant dates from EHR
                                                                TRUE ~ second_transplant_date)) %>% 
                        select(-second_transplant_date_cr)

    rslt$platelet <- rslt$transplant_px %>%
                    left_join(rslt$platelet %>% select(-transplant_date), by = c("person_id", "record_id"))

    record_ids <- rslt$platelet %>% distinct(record_id) %>% arrange(record_id) %>% pull()
    platelet_engraftment_dates <- array(NA, length(record_ids))
    nadir <- array(NA, length(record_ids))
    for (i_record in 1:length(record_ids)){
        re  <- compute_platelet_engraftment(record_ids = record_ids[i_record], 
                                            platelet_var = platelet_max_ct, #either max or mean 
                                            platelet_ct_tbl = rslt$platelet,
                                            cutoff = 20000)
        platelet_engraftment_dates[i_record] = re$platelet_engraf
        nadir[i_record] = re$nadir
    }
    rslt$platelet <- rslt$platelet %>% inner_join(data.frame(record_id = record_ids, 
                                                platelet_engraftment_date = platelet_engraftment_dates, 
                                                nadir = nadir), by = "record_id")
    rslt$platelet %>% view()
    # test case
    compute_platelet_engraftment(record_ids = "01_0021_cr", 
                                platelet_var = platelet_max_ct, #either max or mean 
                                platelet_ct_tbl = rslt$platelet,
                                cutoff = 20000)
    # Determine number of subsets/pages based on the number of unique groups
    num_plot_per_page <- 10
   # missing platelett data
    missing_ids <- c("01_0177_cr", "01_0080_cr", "01_0107_cr", "03_0022_cr", "03_0028_cr", 
                        "04_0037_cr", "05_0007_cr", "05_0008_cr", "05_0004_cr", "05_0025_cr",
                        "06_0021_cr")
    num_pages <- ceiling(rslt$platelet %>% distinct(record_id) %>% 
                        filter(!record_id %in% missing_ids) %>%
                        nrow() / num_plot_per_page)  # Assuming 8 plots per page

    # List to store individual PDF file names
    pdf_files <- list()

    # n = 422 patients
    record_ids <- rslt$platelet %>% distinct(record_id) %>% 
                filter(!record_id %in% missing_ids) %>%
                arrange(record_id) %>% pull()
   
    # Loop through pages and generate plots
    for (i in 1:num_pages) {
        pdf_files[[i]] <- generate_platelet_plots(data = rslt$platelet,
                        transplant_px = rslt$transplant_px, 
                        # transfusion_tbl = rslt$transfusion,
                        page_num = i, num_plot_per_page,
                        record_ids = record_ids)}

    # Combine all PDF files into a single one
    combined_pdf_filename <- "results/platelet_plots/combined_plots_avg_max.pdf"
    pdf_combine(pdf_files, combined_pdf_filename)

    # Remove individual PDF files
    file.remove(pdf_files)

    rslt$platelet %>% distinct(record_id, person_id, nadir, platelet_engraftment_date) %>% 
        output_tbl(name = "platelet_engraftment_dates", local = TRUE, file = TRUE)

}
