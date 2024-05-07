# write a function that works for all cases
# flag the ones that did not work and investigate


# Create a function to generate facet_wrap plots for each patient ID
create_plots <- function(person_id) {
  # Subset data for the current patient ID
  subset_data <- subset(data, person_id == person_id)
  
  # Create facet plot
  plot <- ggplot(subset_data, aes(x = x_variable, y = y_variable)) +
    geom_point() +  # Example plot; replace with your desired plot type
    ggtitle(paste("Patient ID:", patient_id)) +
    theme_bw()
  
  return(plot)
}


compute_anc_engraftment <- function(df, variable, measurement_date, threshold, pivot_date, truncated_date = 100) {
    # if there were multiple measurements on the same day, take the maximum value
    # df <- df %>% group_by(person_id, {{measurement_date}}) %>% 
    #         slice_max({{variable}}, with_ties = FALSE, na_rm = TRUE) %>% ungroup() %>%
    #         group_by(person_id) %>%
    #         mutate(ANC_days_since_ce = as.numeric(difftime({{measurement_date}}, {{pivot_date}}, units = "days"))) 
    
    # find transplant nadir
    nadir_tbl <- df %>% filter({{variable}} < threshold, 
                            ANC_days_since_ce < transplant_2 | is.na(transplant_2), 
                            ANC_days_since_ce <= truncated_date) %>% 
            group_by(person_id) %>% 
            # filter(ANC_days_since_ce == min(ANC_days_since_ce, na.rm = TRUE) | ANC_days_since_ce == max(ANC_days_since_ce, na.rm = TRUE)) %>%
            summarise(nadir = (min(ANC_days_since_ce, na.rm = TRUE) + max(ANC_days_since_ce, na.rm = TRUE))/2, #min(ANC_days_since_ce, na.rm = TRUE),
                  # nadir_index = which.min(ANC_days_since_ce), 
                  nadir_2 = max(ANC_days_since_ce, na.rm = TRUE)) %>% ungroup()
    # nadir should be defined as the average of the 2 min days

  # Iterate over cell counts starting from the day after the nadir
  # only for the first transplant
    df <- df %>% left_join(nadir_tbl, by = "person_id") %>%
                filter({{variable}} >= threshold, ANC_days_since_ce >= nadir, 
                          ANC_days_since_ce < transplant_2 | is.na(transplant_2)) %>%
                          group_by(person_id) %>%
                          mutate(next_ANC = lead(ANC_days_since_ce, order = {{measurement_date}}) - ANC_days_since_ce, 
                          last_ANC = ANC_days_since_ce - lag(ANC_days_since_ce, order = {{measurement_date}}),
                          anc_engraftment_date = ANC_days_since_ce - last_ANC) %>% 
                          filter(next_ANC >= 1, last_ANC >= 1) %>% #, anc_engraftment_date > nadir_2) %>%
                          slice_min(anc_engraftment_date, na_rm = TRUE, with_ties = TRUE) %>% 
                          select(person_id, transplant_date, anc_engraftment_date , {{measurement_date}}, nadir, ANC) %>% ungroup()
    
    return(df)  # If no such day is found, return NA
}

# Function to generate and save plots for each page
generate_ANC_plots <- function(data, transplant_px, page_num, num_plot_per_page = 8, person_ids) {
    # data <- data %>% mutate(days = as.numeric(difftime(measurement_date, transplant_date, units = "days")))
    start_index <- (page_num - 1) * num_plot_per_page + 1
    end_index <- min(page_num * num_plot_per_page, data %>% distinct(person_id) %>% nrow())
    subset_data <- subset(data, person_id %in% person_ids[start_index:end_index]) %>%
                        mutate(ANC = if_else(ANC >= 5000, 5000, ANC))
    # subset_transplant <- subset(transplant_px, person_id %in% person_ids[start_index:end_index]) %>%
    #                 group_by(person_id) %>%
    #                 mutate(transplant_days = as.numeric(difftime(transplant_date, min(transplant_date), units = "days"))) %>%
    #                 summarise(transplant_1 = min(transplant_days),
    #                         transplant_2 = max(transplant_days)) %>%
    #                 select(person_id, transplant_1, transplant_2) %>% ungroup()

    # subset_data <- subset_data %>% inner_join(subset_transplant, by = c("person_id"))
    label_size <- 2
    point_size <- 0.5
    p <- subset_data %>% ggplot(aes(x = ANC_days_since_ce, y = ANC)) + geom_point(size = point_size) + 
    geom_hline(yintercept = 500, linetype = "dashed", color = "red") +  
    geom_label(aes(x = transplant_2, y = 4550, label = "T2"), size = label_size, colour = "green") +
    geom_segment(aes(x = transplant_2, xend = transplant_2, y = 0, yend = 4500), linetype = "dashed", colour = "green") +  
    geom_label(aes(x = transplant_days, y = 4550, label = "T1"), size = label_size, colour = "red") +
    geom_segment(aes(x = transplant_days, xend = transplant_days, y = 0, yend = 4500), linetype = "dashed", colour = "red") +  
    geom_label(aes(x = nadir, y = 4050, label = sprintf("Nadir %d", round(nadir))), size = label_size, colour = "blue") +
    geom_segment(aes(x = nadir, xend = nadir, y = 0, yend = 4000), linetype = "dashed", colour = "blue") +  
    geom_label(aes(x = anc_engraftment_date, y = 3250, label = sprintf("ANC %d", anc_engraftment_date)), size = label_size, colour = "magenta") +
    geom_segment(aes(x = anc_engraftment_date, xend = anc_engraftment_date, y = 0, yend = 3200), linetype = "dashed", colour = "magenta") +  
    facet_wrap(~record_id, ncol = 2, scales = "free") +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 5000))

    plot_filename <- paste0("results/ANC_plots/plots_page_", page_num, ".pdf")
    pdf(plot_filename)
    print(p)
    dev.off()
    return(plot_filename)
}

generate_platelet_plots <- function(data, transplant_px, page_num, num_plot_per_page = 8, person_ids) {
    # data <- data %>% mutate(days = as.numeric(difftime(measurement_date, transplant_date, units = "days")))
    start_index <- (page_num - 1) * num_plot_per_page + 1
    end_index <- min(page_num * num_plot_per_page, data %>% distinct(person_id) %>% nrow())
    subset_data <- subset(data, person_id %in% person_ids[start_index:end_index]) %>%
                        mutate(platelet_ct = if_else(platelet_ct >= 100000, 100000, platelet_ct))
            

    label_size <- 2
    point_size <- 0.5
    cutoff <- 20000
    p <- subset_data %>% ggplot(aes(x = platelet_days_since_ce, y = platelet_ct)) + geom_point(size = point_size) + 
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +  
    # geom_label(aes(x = transplant_2, y = 90000, label = "T2"), size = label_size, colour = "green") +
    # geom_segment(aes(x = transplant_2, xend = transplant_2, y = 0, yend = 90000), linetype = "dashed", colour = "green") +  
    # geom_label(aes(x = transplant_days, y = 90000, label = "T1"), size = label_size, colour = "red") +
    # geom_segment(aes(x = transplant_days, xend = transplant_days, y = 0, yend = 90000), linetype = "dashed", colour = "red") +  
    geom_label(aes(x = transfusion_days_since_ce, y = 50000, label = "BT"), size = label_size, colour = "red") +
    # geom_label(aes(x = nadir, y = 4050, label = sprintf("Nadir %d", round(nadir))), size = label_size, colour = "blue") +
    # geom_segment(aes(x = nadir, xend = nadir, y = 0, yend = 4000), linetype = "dashed", colour = "blue") +  
    # geom_label(aes(x = anc_engraftment_date, y = 3250, label = sprintf("ANC %d", anc_engraftment_date)), size = label_size, colour = "magenta") +
    # geom_segment(aes(x = anc_engraftment_date, xend = anc_engraftment_date, y = 0, yend = 3200), linetype = "dashed", colour = "magenta") +  
    facet_wrap(~record_id, ncol = 2, scales = "free") +
    coord_cartesian(xlim = c(0, 300), ylim = c(0, 100000))

    plot_filename <- paste0("results/platelet_plots/plots_page_", page_num, ".pdf")
    pdf(plot_filename)
    print(p)
    dev.off()
    return(plot_filename)
}
