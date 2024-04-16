
library(tidyr)
library(lubridate)

# message("calculate ANC engraftment dates")
# rslt$ANC <- results_tbl("covar_anc_mx") %>% collect_new()
# rslt$transplant_px <- results_tbl("no_multi_transplant_px") %>% collect_new()

# # Determine number of subsets/pages based on the number of unique groups
# num_pages <- ceiling((rslt$ANC %>% distinct_ct()) / 8)  # Assuming 8 plots per page

# # Function to generate and save plots for each page
# generate_plots <- function(data, transplant_px, page_num) {
# data <- data %>% mutate(days = as.numeric(difftime(measurement_date, transplant_date, units = "days")))
# start_index <- (page_num - 1) * 8 + 1
# end_index <- min(page_num * 8, (rslt$ANC %>% distinct_ct()))
# person_ids <- data %>% distinct(person_id) %>% pull()
# subset_data <- subset(data, person_id %in% person_ids[start_index:end_index])
# subset_transplant <- subset(transplant_px, person_id %in% person_ids[start_index:end_index]) %>%
#                     group_by(person_id) %>%
#                     mutate(transplant_days = as.numeric(difftime(transplant_date, min(transplant_date), units = "days"))) %>%
#                     summarise(transplant_1 = min(transplant_days),
#                             transplant_2 = max(transplant_days)) %>%
#                     select(person_id, transplant_1, transplant_2) %>% ungroup()

# subset_data <- subset_data %>% inner_join(subset_transplant, by = c("person_id"))

# p <- subset_data %>% ggplot(aes(x = days, y = ANC)) + geom_point() + 
#     geom_hline(yintercept = 500, linetype = "dashed", color = "red") +  
#     # geom_vline(xintercept = transplant_1, linetype = "dashed", color = "green") +  
#     # geom_vline(xintercept = transplant_2, linetype = "dashed", color = "green") +  
#     geom_point(aes(x = transplant_1, y = 1000), color = "red", shape = 16, alpha = 0.5) + 
#     geom_point(aes(x = transplant_2, y = 1000), color = "red", shape = 16, alpha = 0.5) + 
#     facet_wrap(~person_id, ncol = 2, scales = "free") +
#     coord_cartesian(xlim = c(0, 50), ylim = c(0, 5000))

# plot_filename <- paste0("results/ANC_plots/plots_page_", page_num, ".pdf")
# pdf(plot_filename)
# print(p)
# dev.off()
# return(plot_filename)
# }

# library(pdftools)

# # List to store individual PDF file names
# pdf_files <- list()

# # Loop through pages and generate plots
# for (i in 1:num_pages) {
# pdf_files[[i]] <- generate_plots(data=rslt$ANC %>% filter(unit_concept_name != "No information"), 
#                                 transplant_px = rslt$transplant_px, i)
# }

# # Combine all PDF files into a single one
# combined_pdf_filename <- "results/ANC_plots/combined_plots.pdf"
# pdf_combine(pdf_files, combined_pdf_filename)

# # Remove individual PDF files
# file.remove(pdf_files)

# get iron overload status
# majority with units nanogram per milliliter
# assume missing units are nanogram per milliliter, no unit conversion done  
cohort_covars <- results_tbl("covar_ferritin_mx") %>% 
            filter((aim == "aim_2a_2" & ferritin_type == "pre") |
                    ((aim == "aim_2b_2" | aim == "aim_3_1" | aim == "aim_3_2") & ferritin_type == "post")|
                     aim == "aim_2a_1" | aim == "aim_2b_1",
                    !is.na(ferritin_level)) %>% collect() %>%
            group_by(person_id, transplant_date, aim) %>% 
            slice_min(abs(difftime(ferritin_date, transplant_date, units = "days")), n = 1, with_ties = FALSE) %>% ungroup() %>% 
            left_join(results_tbl("no_multi_transplant_px") %>% 
                        group_by(person_id) %>%
                        summarise(no_transplants = n_distinct(transplant_date),
                                disease_relapse_date = max(transplant_date, na.rm = TRUE)) %>%
                        filter(no_transplants > 1) %>% 
                        ungroup() %>% mutate(disease_relapse = TRUE) %>% collect(), by = "person_id") %>%
            mutate(disease_relapse = if_else(is.na(disease_relapse), FALSE, disease_relapse))
cohort_covars %>% view()
# a patient could be included in multiple aims
cohort_covars %>% distinct(person_id, aim) %>% count() # n = 1159
cohort_covars %>% distinct_ct() # n = 723
    
# overall survival
# there are 2 patients with more than 1 death causes
# for patients with more than 1 transplant, death days since ce is the number of days since 
# the second transplants
cohort_covars <- cohort_covars %>% copy_to_new(df = ., name = "sdsd") %>%
                        left_join(cdm_tbl("death") %>% 
                                select(person_id, death_date), by = "person_id") %>% collect() %>%
                        mutate(death = if_else(!is.na(death_date), TRUE, FALSE), 
                                death_days_since_ce = if_else(!is.na(death_date), 
                                                    as.numeric(difftime(death_date, transplant_date, units = "days")), NA)) %>%
                        distinct(person_id, aim, death_date, .keep_all = TRUE)
cohort_covars <- cohort_covars %>% mutate(death_days_since_last_transplant = if_else(disease_relapse, 
                                        as.numeric(difftime(death_date, disease_relapse_date, units = "days")),
                                        as.numeric(difftime(death_date, transplant_date, units = "days"))))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, death_date) %>% count() # n = 1159

# GVHD status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_gvhd_dx") %>% collect() %>%
                                mutate(gvhd_days_since_ce = as.numeric(difftime(condition_start_date, transplant_date, units = "days"))) %>%
                                group_by(person_id, transplant_date) %>%
                                filter(gvhd_days_since_ce >= 0) %>%
                                slice_min(abs(gvhd_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                ungroup() %>% select(person_id, gvhd_start_date = condition_start_date, transplant_date, gvhd_days_since_ce), by = c("person_id", "transplant_date")) %>%
                        mutate(GVHD = if_else(!is.na(gvhd_start_date), TRUE, FALSE)) 
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, gvhd_days_since_ce) %>% count() # n = 1159

# VOD status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_vod_dx") %>%
                                    mutate(vod_date = pmax(bilirubin_date, abd_ultrasound_date, defibrotide_date, na.rm = TRUE)) %>%
                                    collect() %>%
                                mutate(vod_days_since_ce = as.numeric(difftime(vod_date, transplant_date, units = "days"))) %>%
                                group_by(person_id, transplant_date) %>%
                                filter(vod_days_since_ce >= 0) %>%
                                slice_min(abs(vod_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                ungroup() %>% 
                                select(person_id, vod_date, transplant_date, vod_days_since_ce), by = c("person_id", "transplant_date")) %>%
                        mutate(VOD = if_else(!is.na(vod_date), TRUE, FALSE))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, vod_date) %>% count() # n = 1159

# bacteremia status
cohort_covars <- cohort_covars %>%
                        left_join(results_tbl("covar_bacteremia_dx") %>% collect() %>%
                                    mutate(bacteremia_days_since_ce = as.numeric(difftime(bacteremia_date, transplant_date, units = "days"))) %>%
                                    group_by(person_id, transplant_date) %>%
                                    filter(bacteremia_days_since_ce >= 0) %>%
                                    slice_min(abs(bacteremia_days_since_ce), n = 1, with_ties = FALSE) %>% 
                                    ungroup() %>% 
                                    select(person_id, bacteremia_date, transplant_date, bacteremia_days_since_ce, bacteremia), by = c("person_id", "transplant_date"))
cohort_covars %>% distinct_ct() # n = 723
cohort_covars %>% distinct(person_id, aim, bacteremia_date) %>% count() # n = 1159

# Immune reconstitution
# results_tbl("covar_CD348_mx") %>% union(results_tbl("covar_IgM_mx") %>% 
#                                         mutate(CDtype = "IgM")) %>% 
#                                 select(person_id, transplant_date, 
#                                         immune_date = measurement_date,
#                                         immune_type = CDtype,
#                                         measurement_concept_id, measurement_concept_name, unit_concept_name, range_high, range_low,
#                                         immune_val = value_as_number)

# cases when patients received both phlebotomy and chelation
cohort_covars <- cohort_covars %>% mutate(val = TRUE) %>% pivot_wider(names_from = aim, values_from = val, values_fill = FALSE) %>%
                mutate(IRT = case_when( aim_3_1 & aim_3_2 ~ "phlebotomy & chelation",
                                        aim_3_1 ~ "phlebotomy", 
                                        aim_3_2 ~ "chelation", 
                                       !aim_3_1 & !aim_3_2 ~ "No IRT", 
                                       TRUE ~ NA)) %>% 
                pivot_longer(cols = aim_2a_2:aim_3_1, names_to = "aim", values_to = "val") %>%
                filter(val) %>% select(-val)
cohort_covars %>%view()
# cohort_covars %>% group_by(person_id) %>% summarise(n = n_distinct(IRT)) %>% filter(n >1) %>%
#                 inner_join(cohort_covars, by = "person_id") %>% 
#                 select(person_id, IRT, aim) %>% 
#                 arrange(person_id, IRT) %>% view()
# cohort_covars %>% filter(!(duplicated(person_id) & IRT == "No IRT")) %>% 
#                 select(person_id, transplant_date, IRT, aim) %>%  arrange(person_id, IRT) %>%  view()
                # filter(!(duplicated(person_id, transplant_date) & (IRT == "phlebotomy" | IRT == "chelation"))) %>% view()
