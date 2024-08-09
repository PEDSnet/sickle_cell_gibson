

.run  <- function() {

        library(tidyr)
        library(lubridate)

        library(pdftools)
        library(ggplot2)

        rslt = list()

        message("get transplant dates from chart review and exclude false positive patients")
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

        message("calculate ANC engraftment dates")
        rslt$ANC <- results_tbl("covar_anc_mx") %>% collect_new()

        # get dates of 2nd transplants (if any)
        # For patients with available chart reviews, use chart reivew transplant dates as index dates
        # Filter out false positives from chart reviews
        rslt$transplant_px <- results_tbl("study_cohorts") %>%
                        filter(!is.na(record_id)) %>% # 21 patients were missed out from chart review, might be included later
                        left_join(results_tbl("cr_data") %>% 
                        select(record_id, record_id, transplant_date_cr = transplant_date, 
                                transplant_type_cr = transplant_type, 
                                second_transplant_date_cr = second_transplant_date,
                                chart_completion, eligibility) %>% distinct(), by = "record_id") %>% 
                        mutate(transplant_date = if_else(!is.na(transplant_date_cr), transplant_date_cr, transplant_date),
                                transplant_type = if_else(!is.na(transplant_type_cr), transplant_type_cr, transplant_type),
                                chart_completion = if_else(is.na(chart_completion), FALSE, chart_completion)) %>%
                        filter((eligibility == "Yes" & chart_completion) | !chart_completion) %>%
                        select(-transplant_date_cr, -transplant_type_cr)
        rslt$transplant_px %>% view()

        rslt$transplant_px <- rslt$transplant_px %>% inner_join(results_tbl("no_multi_transplant_px") %>% 
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
        rslt$transplant_px %>% view()

        # calculate ANC days since 1st transplant
        rslt$ANC <- rslt$transplant_px %>% filter(aim_2a_1 | aim_2a_2) %>%
                left_join(rslt$ANC %>% select(-transplant_date, -site), by = "person_id") %>%
                filter(ANC_date <= transplant_date + days(365), ANC_date >= transplant_date) %>%
                mutate(ANC_days_since_ce = as.numeric(difftime(ANC_date, transplant_date, units = "days")), 
                        second_transplant_days_since_ce = as.numeric(difftime(second_transplant_date, transplant_date, units = "days")))

        # if there were multiple measurements on the same day, take the maximum value
        rslt$ANC <- rslt$ANC %>% group_by(person_id, ANC_date) %>% 
                                slice_max(ANC, n = 1, with_ties = FALSE, na_rm = TRUE) %>% ungroup() 

        # calculate ANC engraftment dates
        anc_engraftment <- compute_anc_engraftment(df = rslt$ANC, 
                                                threshold = 500, variable = ANC, 
                                                measurement_date = ANC_date, 
                                                pivot_date = transplant_date, 
                                                truncated_date = 60,
                                                max_ANC_gaps = 7)
        anc_exception <- compute_anc_engraftment_exception(df = rslt$ANC, 
                                                anc_engraftment = anc_engraftment,
                                                threshold = 500, variable = ANC, 
                                                measurement_date = ANC_date, 
                                                pivot_date = transplant_date, 
                                                truncated_date = 60,
                                                max_ANC_gaps = 7)
        anc_engraftment <- anc_engraftment %>% anti_join(anc_exception, by = c("person_id")) %>% 
                                union(anc_exception)
        # these exceptions would be coded manually: 
        anc_exception2 = data.frame(record_id = c("01_0179_cr", "01_0204_cr", "02_0020_cr", "02_0024_cr", "02_0030_cr", "02_0032_cr", "02_0040_cr",
                                                "02_0076_cr", "02_0078_cr", "02_0127_cr", "02_0138_cr", "03_0032_cr", "03_0043_cr", "03_0056_cr", "03_0071_cr", 
                                                "03_0083_cr", "04_0018_cr", "04_0037_cr", "05_0012_cr", "05_0023_cr", "05_0038_cr", "05_0029_cr", "05_0044_cr", "06_0018_cr",
                                                "06_0021_cr", "07_0002_cr", "07_0003_cr", "07_0004_cr", "07_0031_cr", "07_0086_cr", "07_0090_cr", "08_0002_cr", "08_0005_cr",
                                                "08_0006_cr", "08_0017_cr", "08_0020_cr", "08_0025_cr", "08_0030_cr", "08_0022_cr", "08_0026_cr", "08_0028_cr", "08_0039_cr",
                                                "08_0041_cr", "08_0042_cr", "08_0047_cr", "08_0057_cr", "08_0056_cr", "08_0058_cr", "08_0059_cr", "08_0061_cr"), 
                                   anc_engraftment_date = c(27, 67, 25, 13, 18, 33, 16, 90, 9, 27, 27, 19, 24, 41, 23, 21, 29, 38, 22, 11, 36, 
                                                        23, 22, 6, 0, 16, 23, 25, 15, 30, 34, 11, 11, 13, 10, 15, 12, 11, 19, 12, 28, 15, 13, 36, 22, 10, 13, 19, 11, 22))
        anc_engraftment <-  anc_engraftment %>% left_join(anc_exception2, by = "record_id") %>% 
                            mutate(anc_engraftment_date = if_else(!is.na(anc_engraftment_date.y), anc_engraftment_date.y, anc_engraftment_date.x),
                                anc_engraftment_date = if_else(record_id %in% c("01_0105_cr", "01_0113_cr", "01_0142_cr", 
                                                                                "01_0161_cr", "01_0158_cr", "07_0015_cr", 
                                                                                "07_0061_cr", "07_0058_cr", "07_0087_cr", "01_0052_cr"), NA, anc_engraftment_date),
                                nadir = if_else(is.na(anc_engraftment_date), NA, nadir)) %>% # if a patient doesnt engraft, we wont need their nadir
                            select(-anc_engraftment_date.x, -anc_engraftment_date.y)
        rslt$ANC <- rslt$ANC %>% filter(unit_concept_name != "No information") %>%
                        left_join(anc_engraftment %>% select(-ANC, -ANC_date), by = c("person_id", "transplant_date", "record_id")) %>%
                        arrange(site_id, record_id)

        rslt$ANC %>% view()
        rslt$ANC %>% distinct(person_id) %>% nrow() #659 only aim_2a

        # Determine number of subsets/pages based on the number of unique groups
        # Only for aim_2a
        num_plot_per_page <- 10
        num_pages <- ceiling(rslt$ANC %>% distinct(record_id) %>% nrow() / num_plot_per_page)  

        # List to store individual PDF file names
        pdf_files <- list()
        record_ids <- rslt$ANC %>% distinct(record_id) %>% pull()
        # Loop through pages and generate plots
        # rslt$ANC <- rslt$ANC %>% mutate(anc_engraftment_date = if_else(record_id == "01_0052_cr", NA, anc_engraftment_date),
        #                                 nadir = if_else(record_id == "01_0052_cr", NA, nadir))
        for (i in 1:num_pages) {
                pdf_files[[i]] <- generate_ANC_plots(data=rslt$ANC,
                                page_num = i, num_plot_per_page,
                                record_ids = record_ids)}

        # Combine all PDF files into a single one
        combined_pdf_filename <- "results/ANC_plots/ANC_combined_plots.pdf"
        pdf_combine(pdf_files, combined_pdf_filename)

        # check ANC engraftment dates max

        # Remove individual PDF files
        file.remove(pdf_files)

        # list of patients with no ANC engraftment
        # the easy fix is to do it manually, seriously 
        # "01_0105_cr", "01_0113_cr", "01_0142_cr", "01_0161_cr", "01_0158_cr", "07_0015_cr", "07_0061_cr", "07_0058_cr", "07_0087_cr" not engrafted
        # check units: "01_0161_cr", "01_0158_cr", "01_0188_cr", "01_0201_cr"
        # ask Nora: "01_0191_cr"
        wrong_anc_nadir_ids <- c("01_0179_cr", "01_0204_cr", "02_0020_cr", "02_0024_cr", "02_0030_cr", "02_0032_cr", "02_0040_cr",
                                "02_0076_cr", "02_0078_cr", "02_0127_cr", "02_0138_cr", "03_0032_cr", "03_0043_cr", "03_0056_cr", "03_0071_cr")
        vals = c(27, 67, 25, 13, 18, 33, 16, 90, 9, 27, 27, 19, 24, 41, 23)
        
        # let's fix these first before proceeding with the rest
        anc_engraftment <- compute_anc_engraftment(df = rslt$ANC, 
                                                threshold = 500, variable = ANC, 
                                                measurement_date = ANC_date, 
                                                pivot_date = transplant_date, 
                                                truncated_date = 60,
                                                max_ANC_gaps = 7)
        anc_exception <- compute_anc_engraftment_exception(df = rslt$ANC, 
                                                anc_engraftment = anc_engraftment,
                                                threshold = 500, variable = ANC, 
                                                measurement_date = ANC_date, 
                                                pivot_date = transplant_date, 
                                                truncated_date = 60,
                                                max_ANC_gaps = 7)
        anc_exception %>% arrange(record_id)                                         
        anc_engraftment %>% filter(record_id %in% wrong_anc_nadir_ids) %>% 
                        arrange(record_id) 


        

}