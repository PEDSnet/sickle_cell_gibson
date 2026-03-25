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


compute_anc_engraftment <- function(df, variable, measurement_date, threshold, pivot_date, 
                                    truncated_date = 100, max_ANC_gaps = 7) {
    
    # record_ids <- df %>% distinct(record_id) %>% pull()
        # single_record = "01_0029_cr"
        # person_data <- df %>% filter(record_id %in% c("01_0029_cr", "01_0007_cr", "01_0035_cr"))

        # # # compute nadir
        # nadir_tbl <- person_data %>% filter(ANC < threshold, 
        #                     ANC_days_since_ce < second_transplant_days_since_ce | is.na(second_transplant_days_since_ce), 
        #                     ANC_days_since_ce <= truncated_date) %>%
        #           group_by(person_id, record_id) %>% 
        #           mutate(gaps = lead(ANC_days_since_ce, order = ANC_days_since_ce) - ANC_days_since_ce) %>%
        #           filter(gaps < max_ANC_gaps) %>% 
        #           summarise(nadir = (min(ANC_days_since_ce, na.rm = TRUE) + max(ANC_days_since_ce, na.rm = TRUE))/2, 
        #                     nadir_2 = max(ANC_days_since_ce, na.rm = TRUE)) %>% ungroup() %>%
        #                     view()

        # person_data %>% left_join(nadir_tbl, by = c("record_id", "person_id")) %>%
        #         filter(ANC >= threshold, ANC_days_since_ce >= nadir, 
        #                   ANC_days_since_ce < second_transplant_days_since_ce | is.na(second_transplant_days_since_ce)) %>%
        #                   group_by(person_id) %>%
        #                   mutate(next_ANC = lead(ANC_days_since_ce, order = ANC_date) - ANC_days_since_ce, 
        #                   last_ANC = ANC_days_since_ce - lag(ANC_days_since_ce, order = ANC_date),
        #                   anc_engraftment_date = ANC_days_since_ce - last_ANC) %>% 
        #                   filter(next_ANC >= 1, last_ANC >= 1) %>% #, anc_engraftment_date > nadir_2) %>%
        #                   slice_min(anc_engraftment_date, na_rm = TRUE, with_ties = TRUE) %>% 
        #                   select(person_id, record_id, transplant_date, anc_engraftment_date , ANC_date, nadir, ANC) %>% view()

    # find transplant nadir
    nadir_tbl <- df %>% filter({{variable}} < threshold, 
                            ANC_days_since_ce < second_transplant_days_since_ce | is.na(second_transplant_days_since_ce), 
                            ANC_days_since_ce <= truncated_date) %>% 
                  group_by(person_id, record_id) %>% 
                  mutate(gaps = lead(ANC_days_since_ce, order = ANC_days_since_ce) - ANC_days_since_ce) %>%
                  filter(gaps < max_ANC_gaps) %>% 
                  group_by(person_id, record_id) %>% 
                  summarise(nadir = (min(ANC_days_since_ce, na.rm = TRUE) + max(ANC_days_since_ce, na.rm = TRUE))/2, 
                            nadir_2 = max(ANC_days_since_ce, na.rm = TRUE)) %>% ungroup()

    # nadir should be defined as the average of the 2 min days

  # Iterate over cell counts starting from the day after the nadir
  # only for the first transplant
    df <- df %>% left_join(nadir_tbl, by = c("record_id", "person_id")) %>%
                filter({{variable}} >= threshold, ANC_days_since_ce >= nadir, 
                          ANC_days_since_ce < second_transplant_days_since_ce | is.na(second_transplant_days_since_ce)) %>%
                group_by(person_id, record_id) %>%
                mutate(next_ANC = lead(ANC_days_since_ce, order = {{measurement_date}}) - ANC_days_since_ce, 
                      last_ANC = ANC_days_since_ce - lag(ANC_days_since_ce, order = {{measurement_date}}),
                      anc_engraftment_date = ANC_days_since_ce - last_ANC) %>% 
                filter(next_ANC >= 1, last_ANC >= 1) %>% #, anc_engraftment_date > nadir_2) %>%
                slice_min(anc_engraftment_date, na_rm = TRUE, with_ties = TRUE) %>% 
                select(person_id, record_id, transplant_date, anc_engraftment_date , {{measurement_date}}, nadir, ANC, nadir_2) %>% ungroup()
    
    return(df)  # If no such day is found, return NA
}

compute_anc_engraftment_exception <- function(df, anc_engraftment, variable, measurement_date, threshold, pivot_date, 
                                    truncated_date = 100, max_ANC_gaps = 7) {
    
    exceptions <- anc_engraftment %>% filter(anc_engraftment_date < nadir_2) %>% distinct(record_id, person_id, anc_engraftment_date)

    # find transplant nadir
    df <- df %>% inner_join(exceptions, by = c("record_id", "person_id")) %>%
                filter(ANC_days_since_ce >= anc_engraftment_date) %>%
                compute_anc_engraftment(df = ., 
                                        threshold = threshold, variable = {{variable}}, 
                                        measurement_date = {{measurement_date}}, 
                                        pivot_date = {{pivot_date}}, 
                                        truncated_date = truncated_date,
                                        max_ANC_gaps = max_ANC_gaps)
    return(df)  # If no such day is found, return NA
}

# Define custom function to generate plot titles
custom_labeller <- function(variable, value) {
  group_data <- df[df$group == value, ]
  return(paste(value, "  cr_completed = ", group_data$chart_completion))
}

# Function to generate and save plots for each page
generate_ANC_plots <- function(data, page_num, num_plot_per_page = 8, record_ids) {
    # data <- data %>% mutate(days = as.numeric(difftime(measurement_date, transplant_date, units = "days")))
    start_index <- (page_num - 1) * num_plot_per_page + 1
    end_index <- min(page_num * num_plot_per_page, data %>% distinct(record_id) %>% nrow())
    subset_data <- subset(data, record_id %in% record_ids[start_index:end_index]) %>%
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
    # red dashed line for threshold
    geom_hline(yintercept = 500, linetype = "dashed", color = "red") +  
    # label second transplant
    geom_label(aes(x = second_transplant_days_since_ce, y = 4550, label = "T2"), size = label_size, colour = "green") +
    geom_segment(aes(x = second_transplant_days_since_ce, xend = second_transplant_days_since_ce, y = 0, yend = 4500), linetype = "dashed", colour = "green") +  
    # label first transplant
    annotate("text", x = 0, y = 4550, label = "T1", size = label_size, colour = "red") +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 4500, linetype = "dashed", colour = "red") +
    # geom_label(aes(x = 0, y = 4550, label = "T1"), size = label_size, colour = "red") +
    # geom_segment(aes(x = 0, xend = 0, y = 0, yend = 4500), linetype = "dashed", colour = "red") +  
    # label nadir
    geom_label(aes(x = nadir, y = 4050, label = sprintf("Nadir %d", round(nadir))), size = label_size, colour = "blue") +
    geom_segment(aes(x = nadir, xend = nadir, y = 0, yend = 4000), linetype = "dashed", colour = "blue") +  
    # label ANC engraftment
    geom_label(aes(x = anc_engraftment_date, y = 3250, label = sprintf("ANC %d", anc_engraftment_date)), size = label_size, colour = "magenta") +
    geom_segment(aes(x = anc_engraftment_date, xend = anc_engraftment_date, y = 0, yend = 3200), linetype = "dashed", colour = "magenta") +  
    facet_wrap(~record_id, ncol = 2, scales = "free", labeller = labeller(group = custom_labeller)) +
    # labs(title = paste0(record_id, "cr_completion = ", chart_completion)) + 
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 5000))

    plot_filename <- paste0("results/ANC_plots/plots_page_", page_num, ".pdf")
    pdf(plot_filename)
    print(p)
    dev.off()
    return(plot_filename)
}

generate_platelet_plots <- function(data, 
                                    transplant_px, 
                                    transfusion_tbl,
                                    page_num, num_plot_per_page = 8, record_ids) {
    
    start_index <- (page_num - 1) * num_plot_per_page + 1
    end_index <- min(page_num * num_plot_per_page, data %>% distinct(record_id) %>% nrow())
    subset_data <- subset(data, record_id %in% record_ids[start_index:end_index]) %>%
                        mutate(platelet_max_ct = if_else(platelet_max_ct >= 100000, 100000, platelet_max_ct),
                              platelet_mean_ct = if_else(platelet_mean_ct >= 100000, 100000, platelet_mean_ct),
                              second_transplant_days_since_ce = as.numeric(difftime(second_transplant_date, transplant_date, units = "days")))
                            
            

    label_size <- 2
    point_size <- 0.5
    cutoff <- 20000
    
    ann_text <- subset_data %>% distinct(record_id, nadir, platelet_engraftment_date, second_transplant_days_since_ce)

    p <- subset_data %>% ggplot(aes(x = platelet_days_since_ce)) + 
    geom_point(aes(y = platelet_mean_ct), size = point_size, colour = "red") + 
    # geom_point(aes(y = platelet_max_ct), size = point_size, colour = "black") + 
    geom_point(aes(y = 50000, x = transfusion_days_since_ce), size = point_size, colour = "black") + 
    # geom_label(aes(x = transfusion_days_since_ce, y = 50000, label = "BT"), size = label_size, colour = "red") +
    facet_wrap(~record_id, ncol = 2, scales = "free") 

    p <- p + 
        geom_hline(yintercept = cutoff, linetype = "dashed", color = "red") +  
        geom_label(data = ann_text, aes(x = second_transplant_days_since_ce, y = 80000, label = "T2"), size = label_size, colour = "green") +
        geom_segment(data = ann_text, aes(x = second_transplant_days_since_ce, xend = second_transplant_days_since_ce, y = 0, yend = 80000), linetype = "dashed", colour = "green") +  
        # geom_label(aes(x = transplant_days, y = 90000, label = "T1"), size = label_size, colour = "red") +
        # geom_segment(aes(x = transplant_days, xend = transplant_days, y = 0, yend = 90000), linetype = "dashed", colour = "red") +  
        geom_label(data = ann_text, aes(x = nadir, y = 80000, label = sprintf("Nadir %d", round(nadir))), size = label_size, colour = "blue") +
        geom_segment(data = ann_text, aes(x = nadir, xend = nadir, y = 0, yend = 80000), linetype = "dashed", colour = "blue") +  
        geom_label(data = ann_text, aes(x = as.numeric(platelet_engraftment_date), y = 70000, label = sprintf("platelet %s", platelet_engraftment_date)), size = label_size, colour = "magenta") +
        geom_segment(data = ann_text, aes(x = as.numeric(platelet_engraftment_date), xend = as.numeric(platelet_engraftment_date), y = 0, yend = 70000), linetype = "dashed", colour = "magenta") +
        # scale_x_continuous(minor_breaks = seq(-5, 100, by = 5),  # Minor grid every 5
        # breaks = seq(0, 100, by = 20) +           # Major grid every 2
        # theme(panel.grid.minor = element_line(color = "gray", linetype = "dotted"))) +
        coord_cartesian(xlim = c(-5, 100), ylim = c(0, 100000))
    
    plot_filename <- paste0("results/platelet_plots/plots_page_", page_num, ".pdf")
    pdf(plot_filename)
    print(p)
    dev.off()
    return(plot_filename)
}

#' Determine initial platelet engraftment day
#'
#' Takes a tibble with labs, a patient MRN, a tibble with platelet transfusions,
#' and a transplant date, and returns the day post transplant on which Plt20
#' was reached.
#'
#' Rules for Platelet engraftment:
#'
#' 1. First of 3 consecutive Plt > 20 values obtained on separate days
#'
#' 2. Don't start until nadir post-HCT (local minimum with Plt < 20 within
#'    2 weeks post transplant)
#'
#' 3. Ignore results within 7d post platelet transfusion (until including day +6)
#'
#' 4. Plt50 is the first of 3 consecutive Plt > 50 on or after Plt20 day
#'
#' @param all_infusion_dates A vector of all transplant dates for this patient,
#' so if engraftment date is after subsequent transplant, count as "never
#' engrafted".
#' @param platelet_transfusion_dates Vector of platelet transfusion dates. Only necessary for \code{plt50_day} to feed a vector from the plt50 functional environment so \code{first_transfusion_day} passes in plt20 code. Defaults to NA.
#'
#' @import rlang lubridate purrr
#' @export
compute_platelet_engraftment <- function(record_ids, 
                                        platelet_var, #either max or mean
                                        cutoff = 20000,
                                        platelet_ct_tbl
) {

  platelet_ct_tbl <- subset(platelet_ct_tbl, record_id %in% record_ids) %>% 
                              filter(platelet_days_since_ce >= -7)
  t <- platelet_ct_tbl %>% filter(platelet_days_since_ce < 0, !is.na(transfusion_days_since_ce)) %>% 
            pull(transfusion_days_since_ce)  

  if(is_empty(t)){
    # no transfusion before transplant
    platelet_ct_tbl <- platelet_ct_tbl %>% filter(platelet_days_since_ce >=0)
  } else{
    platelet_ct_tbl <- platelet_ct_tbl %>% filter(platelet_days_since_ce >= (min(t, na.rm = TRUE) + 7))
  }

  platelet_transfusion_dates <- subset(platelet_ct_tbl, !is.na(transfusion_days_since_ce)) %>% 
                          arrange(transfusion_days_since_ce) %>% 
                          pull(transfusion_days_since_ce)
  
  platelet_ct_tbl <- platelet_ct_tbl %>%
      mutate(
        invalid = map_lgl(
          platelet_days_since_ce,
          ~ if (
            any(
              seq(
                from = .x - 6,
                to = .x,
                by = 1
              ) %in% platelet_transfusion_dates
            )
          ) {
            TRUE
          } else {
            FALSE
          }
        )
      )
  # platelet_ct_tbl %>% view()
  #
  # Find first local minimum < 20
  # Fix for #29: allow for nadir to be on day 0
  #

  platelet_ct_tbl <- platelet_ct_tbl %>%
    mutate(
      nadir = ifelse(
        ({{platelet_var}} <= lag({{platelet_var}}) | platelet_days_since_ce == 0) &
          {{platelet_var}} <= lead({{platelet_var}}) &
          {{platelet_var}} < cutoff,
        # (platelet_max_ct<= lag(platelet_max_ct) | platelet_days_since_ce == 0) &
        #   platelet_max_ct <= lead(platelet_max_ct) &
        #   platelet_max_ct < cutoff,
        TRUE,
        FALSE
      )
    )

  # Handle edge case with no labs
  if (nrow(platelet_ct_tbl) == 0) {
    # warning(glue("plt20_day(): no platelet counts for {patient_mrn} found between {start_date} and {end_date}\n"))
    # return(NA_real_)
    return(list(nadir = NA, platelet_engraft = "No labs"))
  }

  nadir_day <- platelet_ct_tbl %>% filter(nadir == TRUE) %>% slice(1) %>%
    pull(platelet_days_since_ce)

  # Handle edge case where first value is lowest value and <20
  temp = platelet_ct_tbl %>% pull({{platelet_var}}) 
  # temp = platelet_ct_tbl %>% pull(platelet_max_ct) 
  # platelet_days_since_ce_0 <- platelet_ct_tbl %>% filter(row_number(platelet_days_since_ce) == which.min(platelet_max_ct)) %>% 
  platelet_days_since_ce_0 <- platelet_ct_tbl %>% filter(row_number(platelet_days_since_ce) == which.min({{platelet_var}})) %>% 
              pull(platelet_days_since_ce)
  platelet_ct_0 <- platelet_ct_tbl %>% slice(1) %>% pull({{platelet_var}})
  if ( platelet_ct_0 == min(temp, na.rm = TRUE) & platelet_ct_0 < cutoff) {
    nadir_day <- platelet_ct_tbl %>% filter({{platelet_var}} == platelet_ct_0) %>% slice(1) %>% pull(platelet_days_since_ce)
  }

  # Handle platelet transfusions after transplant
  # by setting day of first transfusion after txp
  # as the smaller of nadir_day and transfusion day
  first_transfusion_day <- platelet_ct_tbl %>% filter(platelet_days_since_ce %in% platelet_transfusion_dates) %>%
                                                slice(1) %>% pull(platelet_days_since_ce)

  # if (length(first_transfusion_day) > 0) {
  #   nadir_day <- max(nadir_day, first_transfusion_day) #????
  # }

  # Handle no nadir found or nadir too late (>3 wk post txp)
  if (length(nadir_day) == 0) {  #|| nadir_day > 21
    if(length(first_transfusion_day) > 0 && first_transfusion_day == 0){
      # platelet count never below cutoff and transfusion on day 0
      return(list(nadir = 0, platelet_engraft = "7"))
    } else {
      return(list(nadir = NA, platelet_engraft = "0")) #"Never below"
    }
  }
  

  nadir_end <- nadir_day
  check_engraftment <- TRUE
  while(check_engraftment){
    # platelet_engraf_date <- platelet_ct_tbl %>% filter(platelet_days_since_ce > nadir_end) %>%
    # # Are there 3 consecutive VALID days with PLT >= 20?
    # mutate(plt20 = ifelse(platelet_max_ct >= cutoff & invalid == FALSE &
    #     lead(n = 1, platelet_max_ct) >= cutoff &
    #     lead(n = 1, invalid) == FALSE &
    #     lead(n = 2, platelet_max_ct) >= cutoff &
    #     lead(n = 2, invalid) == FALSE,
    #   TRUE,
    #   FALSE
    # )) %>% filter(plt20 == TRUE) %>% slice(1) %>%
    # pull(platelet_days_since_ce)
    platelet_engraf_date <- platelet_ct_tbl %>% filter(platelet_days_since_ce > nadir_end) %>%
    # Are there 3 consecutive VALID days with PLT >= 20?
    mutate(plt20 = ifelse({{platelet_var}} >= cutoff & invalid == FALSE &
        lead(n = 1, {{platelet_var}}) >= cutoff &
        lead(n = 1, invalid) == FALSE &
        lead(n = 2, {{platelet_var}}) >= cutoff &
        lead(n = 2, invalid) == FALSE,
      TRUE,
      FALSE
    )) %>% filter(plt20 == TRUE) %>% slice(1) %>%
    pull(platelet_days_since_ce)

    #after engraftment does the value dip again with 5 days? If yes, shift nadir_end to the engraftment date
    if(length(platelet_engraf_date) == 0){
      return(list(nadir = nadir_day, platelet_engraft = "Never engrafted"))
    } else {
      t1 <- platelet_ct_tbl %>% filter(platelet_days_since_ce >= platelet_engraf_date, 
                                    platelet_days_since_ce <= (platelet_engraf_date + 5),
                                    {{platelet_var}} < cutoff)
      if(nrow(t1)==0){
        check_engraftment <- FALSE
      } else{
        nadir_end <- platelet_engraf_date
        # print(nadir_end)
      }
    }
    
  }
  

  # Fix for Issue #1: if engraftment date is after subsequent transplant,
  # count as "never engrafted"
  # all_infusion_dates vextor of transplant dates

  if (is.na(distinct(select(platelet_ct_tbl, second_transplant_date)))) {
    next_transplant_day <- platelet_ct_tbl %>% filter(!is.na(second_transplant_date)) %>%
                          pull(platelet_days_since_ce)

    if (length(next_transplant_day) == 0) {
      next_transplant_day <- NA_real_
    } else {
      next_transplant_day <- min(next_transplant_day, na.rm = TRUE)
    }

    if (!is.na(next_transplant_day) &&
      length(platelet_engraf_date) != 0 && platelet_engraf_date > next_transplant_day) {
      # return(NA_real_)
      return(list(nadir = nadir_day, platelet_engraf = platelet_engraf_date))
    }
  }

  if (is.null(platelet_engraf_date) || length(platelet_engraf_date) == 0) {
    # NA_real_
    return(list(nadir = nadir_day, platelet_engraft = "Never engrafted"))
  } else {
    return(list(nadir = nadir_day, platelet_engraft = as.character(platelet_engraf_date)))
  }
  
}
