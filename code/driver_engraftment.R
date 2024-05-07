library(tidyr)
library(lubridate)

library(pdftools)
library(ggplot2)

rslt = list()
message("calculate ANC engraftment dates")
rslt$ANC <- results_tbl("covar_anc_mx") %>% collect_new()
rslt$transplant_px <- results_tbl("no_multi_transplant_px") %>% collect_new() 
rslt$transplant_px <- rslt$transplant_px %>% group_by(person_id) %>%
                    mutate(transplant_days = as.numeric(difftime(transplant_date, min(transplant_date), units = "days"))) %>%
                    summarise(transplant_1 = min(transplant_days),
                            transplant_2 = max(transplant_days),
                            transplant_date = min(transplant_date)) %>%
                    mutate(transplant_2 = if_else(transplant_2 == transplant_1, NA, transplant_2)) %>%
                    select(person_id, transplant_date, transplant_days = transplant_1, transplant_2) %>% ungroup()

# calculate ANC days since 1st transplant
rslt$ANC <- rslt$ANC %>% select(person_id, measurement_date, ANC, transplant_date, unit_concept_name) %>%
                        inner_join(rslt$transplant_px, by = c("person_id", "transplant_date"))

# For patients with available chart reviews, use chart reivew transplant dates as index dates
rslt$study_cohorts <- results_tbl("study_cohorts") %>% collect() %>%
                        left_join(results_tbl("cr_cohorts") %>% 
                                select(person_id, record_id, transplant_date, transplant_type, second_transplant_date) %>% collect(), by = "person_id")

rslt$ANC <- results_tbl("study_cohorts") %>% collect() %>%
                left_join(rslt$ANC, by = c("person_id", "transplant_date")) %>%
                left_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>% select(person_id, record_id) %>% collect(), by = "person_id")

# if there were multiple measurements on the same day, take the maximum value
rslt$ANC <- rslt$ANC %>% group_by(person_id, measurement_date) %>% 
        slice_max(ANC, with_ties = FALSE, na_rm = TRUE) %>% ungroup() %>%
        group_by(person_id) %>%
        mutate(ANC_days_since_ce = as.numeric(difftime(measurement_date, transplant_date, units = "days"))) 

# calculate ANC engraftment dates
anc_engraftment <- compute_anc_engraftment(df = rslt$ANC, threshold = 500, variable = ANC, 
                                                measurement_date, pivot_date = transplant_date, truncated_date = 80)

rslt$ANC <- rslt$ANC %>% filter(unit_concept_name != "No information") %>%
                left_join(anc_engraftment %>% select(-ANC, -measurement_date), by = c("person_id", "transplant_date")) %>%
                arrange(record_id)
rslt$ANC %>% view()
rslt$ANC %>% distinct(person_id) %>% nrow() #829
rslt$ANC %>% distinct(record_id) %>% nrow() #829

# Determine number of subsets/pages based on the number of unique groups
num_plot_per_page <- 10
num_pages <- ceiling(rslt$ANC %>% distinct(record_id) %>% nrow() / num_plot_per_page)  # Assuming 8 plots per page

# List to store individual PDF file names
pdf_files <- list()

# Loop through pages and generate plots
for (i in 1:num_pages) {
        pdf_files[[i]] <- generate_ANC_plots(data=rslt$ANC,
                                transplant_px = rslt$transplant_px, 
                                page_num = i, num_plot_per_page,
                                person_ids = rslt$ANC %>% distinct(person_id) %>% pull())
}

# Combine all PDF files into a single one
combined_pdf_filename <- "results/ANC_plots/combined_plots.pdf"
pdf_combine(pdf_files, combined_pdf_filename)

# Remove individual PDF files
file.remove(pdf_files)

# list of patients with no ANC engraftment
no_anc_ids <- c("07_0042_cr", "07_0044_cr", "07_0019_cr", "07_0048_cr", "07_0023_cr", "07_0054_cr", 
                "07_0053_cr", "07_0055_cr", "06_0021_cr", "03_0017_cr", "02_0034_cr", "02_0117_cr",
                "02_0049_cr", "02_0120_cr", "02_0043_cr", "02_0043_cr", "05_0016_cr","01_0201_cr",
                )
wrong_anc_nadir_ids <- c("07_0050_cr", "07_0052_cr", "05_0023_cr", "05_0029_cr", "05_0038_cr",
                        "06_0018_cr", "06_0012_cr", "03_0005_cr", "03_0074_cr", "03_0014_cr", "03_0021_cr",
                        "03_0032_cr", "03_0090_cr", "03_0034_cr", "02_0086_cr", "02_0090_cr", "02_0109_cr",
                        "02_0026_cr", "02_0112_cr", "02_0110_cr", "02_0029_cr", "02_0038_cr", "02_0116_cr",
                        "02_0039_cr", "02_0044_cr", "02_0048_cr", "02_0122_cr", "05_0044_cr", "07_0060_cr",
                        "07_0082_cr", "01_0087_cr", "01_0154_cr", "01_0191_cr", "01_0045_cr", "01_0140_cr",
                        "01_0009_cr")
no_data_ids <- c("07_0055_cr", "05_0006_cr", "03_0088_cr", "05_0043_cr", "07_0058_cr", "07_0079_cr")
verify_with_nora <- c("05_0012_cr", "03_0002_cr", "NA", "03_0033_cr", "03_0040_cr", "03_0036_cr", "03_0039_cr", "03_0096_cr", 
                        "03_0043_cr", "03_0047_cr", "03_0102_cr", "02_0093_cr", "02_0096_cr", "02_0020_cr", "02_0024_cr", 
                        "02_0033_cr", "02_0032_cr", "02_0057_cr", "03_0056_cr", "02_0145_cr", "02_0074_cr")
verify_transplant_status <- c("02_0131_cr")

# let's fix these first before proceeding with the rest



message("Compute platelet engraftment dates")
# rslt$platelet <- get_outcome_measurements(cohort = results_tbl("study_cohorts"), 
#                                                 mx_codeset = load_codeset("platelet_mx"),
#                                                 is_pre = FALSE, is_closet = FALSE) 

# check units
# thousand per microliter = per cubic millimeter = Kelvin per microliter (this K actually means thousands)
# the cut-off is 20K/microliter
# for measurement_concept_id = 3007461 with different ranges, it seems like we could ignore them because 
# some of the high ranges are not relevant. e.g no values between 300 and 750
# exclude the entitic volume as it measures the volume rather than the counts

platelet %>% group_by(measurement_concept_id, 
                                measurement_concept_name, 
                                unit_concept_name,
                                range_high, range_low) %>% summarise(n = n()) %>% view()

# rslt$platelet %>% output_tbl(name = "covar_platelet_mx_original")
# only look at aim 2a_1 and 2a_2
platelet <- results_tbl("covar_platelet_mx_original") %>% 
                filter(aim_2a_1 | aim_2a_2) %>% collect()
platelet %>% view()   

# if we have a lot of patients with missing platelet counts, we might need to check the filtered data
platelet %>% filter(measurement_concept_id %in% c( "3007461", "3024929", "4267147"),
                            !is.na(value_as_number)) %>% 
                        mutate(value_as_number = value_as_number, unit_concept_name = "per microliter") %>%
                        output_tbl(name = "covar_platelet_mx")
                                    

# platelet transfusion procedure 
rslt$platelet_transfusion_pre <- get_outcome_procedures(results_tbl("study_cohorts"),
                                                px_codeset = load_codeset("platelet_transfusion_px"),
                                                is_pre = TRUE, is_closet = FALSE)
rslt$platelet_transfusion_pre %>% output_tbl(name = "covar_platelet_transfusion_pre_px")
rslt$platelet_transfusion %>% output_tbl(name = "covar_platelet_transfusion_post_px")
results_tbl("covar_platelet_transfusion_pre_px") %>% union(results_tbl("covar_platelet_transfusion_post_px")) %>% 
        collect_new() %>% filter(!is.na(procedure_date)) %>%
        output_tbl(name = "covar_platelet_transfusion_px")

rslt$platelet_transfusion %>% view()
        
# calculate platelet engraftment date
tbl_platelet <- results_tbl("covar_platelet_mx") %>% collect_new() %>% 
        filter(aim_2a_1 | aim_2a_2) %>% 
        group_by(person_id, transplant_date, measurement_date) %>%
        # summarise(value_as_number = max(value_as_number, na.rm = TRUE)) %>%
        mutate(platelet_days_since_ce = as.numeric(difftime(measurement_date, transplant_date, units = "days")),
                val_ct = ifelse(value_as_number >= 20000, 1, 0)) %>% 
        mutate(platelet_consec = ifelse((lead(platelet_days_since_ce, 1) - platelet_days_since_ce == 1) & 
                                        (lead(platelet_days_since_ce, n =2) - platelet_days_since_ce == 2) &
                                        (val_ct + lead(val_ct, 1) + lead(val_ct, n =2)) == 3, TRUE, FALSE)) %>% 
        ungroup() %>% 
        # filter(platelet_consec == TRUE) %>%
        select(person_id, transplant_date, platelet_measurement_date = measurement_date, #ce_date, birth_date, gender, 
                platelet_ct = value_as_number, platelet_consec, platelet_days_since_ce)

# tbl_transfusion <- results_tbl("covar_platelet_transfusion_px") %>% collect_new() %>% 
#         filter(aim_2a_1 | aim_2a_2) %>% 
#         mutate(transfusion_days_since_ce = as.numeric(difftime(procedure_date, transplant_date, units = "days"))) %>%
#         distinct(person_id, transplant_date, procedure_date, .keep_all= TRUE) %>%
#         select(person_id, transplant_date, transfusion_date = procedure_date, transfusion_days_since_ce)
tbl_transfusion <- results_tbl("cdm_procedure_occurrence") %>% 
        inner_join(results_tbl("study_cohorts") %>% 
                        filter(aim_2a_1 | aim_2a_2) %>% 
                        distinct(person_id, transplant_date), by = "person_id") %>% 
        filter(procedure_date <= transplant_date + days(360), 
                procedure_date >= transplant_date - days(180),
                grepl("transfusion", procedure_concept_name, ignore.case = TRUE),
                grepl("blood", procedure_concept_name, ignore.case = TRUE) | 
                grepl("platelet", procedure_concept_name, ignore.case = TRUE)) %>% collect() %>%
                mutate(transfusion_days_since_ce = as.numeric(difftime(procedure_date, transplant_date, units = "days"))) %>%
                distinct(person_id, transplant_date, procedure_date, .keep_all= TRUE) %>%
                select(person_id, transplant_date, transfusion_date = procedure_date, transfusion_days_since_ce)
tbl_transfusion %>% view()
rslt$ANC %>% select(person_id, ce_date, transplant_date, birth_date, gender, ANC_measurement_date = measurement_date,
                ANC, transplant_days, transplant_2, record_id, ANC_days_since_ce, anc_engraftment_date, nadir)
rslt$platelet <- tbl_platelet %>%
                left_join(tbl_transfusion, by = c("person_id", "transplant_date")) %>% 
                left_join(results_tbl("master_xwalk_ids_ALL_AIMS") %>% select(person_id, record_id) %>% collect(), by = "person_id")

# Determine number of subsets/pages based on the number of unique groups
num_plot_per_page <- 10
num_pages <- ceiling(rslt$platelet %>% distinct(record_id) %>% nrow() / num_plot_per_page)  # Assuming 8 plots per page

# List to store individual PDF file names
pdf_files <- list()

# n = 723 patients
person_ids <- rslt$platelet %>% distinct(person_id) %>% arrange(person_id) %>% pull()

# Loop through pages and generate plots
for (i in 1:num_pages) {
        pdf_files[[i]] <- generate_platelet_plots(data=rslt$platelet,
                                transplant_px = rslt$transplant_px, 
                                page_num = i, num_plot_per_page,
                                person_ids = person_ids)
}

# Combine all PDF files into a single one
combined_pdf_filename <- "results/platelet_plots/combined_plots.pdf"
pdf_combine(pdf_files, combined_pdf_filename)

# Remove individual PDF files
file.remove(pdf_files)

# how many transplant patients who never had a platelet transfusion
# 723 patients in aim 2a_1 and 2a_2
# n = 83 patients in aim 2a_1 and 2a_2 never had platelet transfusion???
no_transfusion_id <- results_tbl("study_cohorts") %>% 
        filter(aim_2a_1 | aim_2a_2) %>% 
        distinct(person_id, transplant_date) %>% 
        anti_join(results_tbl("covar_platelet_transfusion_px") %>% 
                        distinct(person_id, transplant_date), by = c("person_id", "transplant_date")) 
# did we miss sth?


# 266
rslt$platelet %>% distinct_ct() 
rslt$platelet %>% distinct(procedure_concept_name) %>% view()   

        # distinct(person_id, transplant_date, procedure_concept_id, procedure_concept_name, .keep_all = TRUE) %>%
        select(person_id, transplant_date, procedure_concept_id, procedure_concept_name, procedure_date) %>% 
        arrange(person_id, transplant_date, procedure_date) %>% view()
