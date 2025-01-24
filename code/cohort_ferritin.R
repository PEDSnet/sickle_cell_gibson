ferritin_classification <- function(cohort, 
                                ferritin_cutoff1 = 1000, ferritin_cutoff2 = 2500,
                                ferritin_type, no_ferritin = 1, keep_tie = FALSE, 
                                cutoff_days = 365.25,
                                slice_by = "most_recent") {

        cohort <- cohort %>% collect() %>% group_by(person_id, transplant_date) %>%
                        mutate(ferritin_days_since_transplant = as.numeric(difftime(ferritin_date, transplant_date, units = "days"))) %>%
                        filter(abs(ferritin_days_since_transplant) <= cutoff_days) %>% ungroup()

        if (ferritin_type == "pre") { # most recent values 
                cohort <- cohort %>% filter(ferritin_days_since_transplant <= 0) %>% 
                        group_by(person_id, transplant_date) %>% 
                        slice_min(as.numeric(difftime(transplant_date, ferritin_date, units = "days")), n = no_ferritin, with_ties = keep_tie) %>% 
                        summarise(ferritin = mean(ferritin, na.rm = TRUE),
                        ferritin_days_since_transplant = mean(ferritin_days_since_transplant, na.rm = TRUE)) %>% ungroup() %>%
                        mutate(ferritin_type = "pre")                                        
        } else if (ferritin_type == "post") { # most recent values 
                cohort <- cohort %>% filter(ferritin_days_since_transplant > 0) %>% 
                        group_by(person_id, transplant_date) %>% 
                        slice_min(ferritin_days_since_transplant, n = no_ferritin, with_ties = keep_tie) %>% 
                        summarise(ferritin = mean(ferritin, na.rm = TRUE), 
                                 ferritin_days_since_transplant = mean(ferritin_days_since_transplant, na.rm = TRUE)) %>% 
                        ungroup() %>% mutate(ferritin_type = "post")
        } else if( ferritin_type == "pre_max") {
                cohort <- cohort %>% filter(ferritin_days_since_transplant <= 0) %>% 
                        group_by(person_id, transplant_date) %>% 
                        slice_max(ferritin, n = no_ferritin, with_ties = keep_tie) %>% 
                        summarise(ferritin = mean(ferritin, na.rm = TRUE),
                        ferritin_days_since_transplant = mean(ferritin_days_since_transplant, na.rm = TRUE)) %>% 
                        ungroup() %>% mutate(ferritin_type = "pre_max")
        } else if( ferritin_type == "post_max") {
                cohort <- cohort %>% filter(ferritin_days_since_transplant > 0) %>% 
                        group_by(person_id, transplant_date) %>% 
                        slice_max(ferritin, n = no_ferritin, with_ties = keep_tie) %>% 
                        summarise(ferritin = mean(ferritin, na.rm = TRUE),
                        ferritin_days_since_transplant = mean(ferritin_days_since_transplant, na.rm = TRUE)) %>% 
                        ungroup() %>% mutate(ferritin_type = "post_max")
        }
        cohort <- cohort %>% mutate(ferritin_level = case_when(ferritin < ferritin_cutoff1 ~ "low",
                                                        ferritin >= ferritin_cutoff1 & ferritin < ferritin_cutoff2 ~ "moderate",
                                                        ferritin >= ferritin_cutoff2 ~ "high",
                                                        TRUE ~ NA)) 
        return(cohort)
        
}

val_classification <- function(cohort, val, label, cutoffs){
    
    if(length(cutoffs) == 1){
        # only 2 levels 
        cohort <- cohort %>% mutate({{label}} := if_else({{val}} < cutoffs[1], "low", "high"))
    } else if(length(cutoffs) == 2){
        cohort <- cohort %>% mutate({{label}} := case_when({{val}} < cutoffs[1] ~ "low",
                                                        {{val}} >= cutoffs[1] & {{val}} <= cutoffs[2] ~ "moderate",
                                                        {{val}} > cutoffs[2] ~ "high",
                                                        TRUE ~ NA)) 
    }
    cohort <- cohort %>% mutate({{label}} := factor({{label}}, levels = c("low", "moderate", "high")))
    return(cohort)
}

val_extraction <- function(cohort, no_value = 1, 
                           grouping_id = person_id, 
                           value_name = LIC, 
                           value_date = LIC_date, 
                           value_type = "LIC", 
                           keep_tie = FALSE, 
                           cutoff_days = 365.25, cutoffs, 
                           slice_by = "most_recent") {
                                
        cohort <- cohort %>% collect() %>% 
                        mutate(!!sym(paste0(value_type, "_days_since_transplant")) := as.numeric(difftime({{value_date}}, transplant_date, units = "days"))) %>%
                        filter(abs(!!sym(paste0(value_type, "_days_since_transplant"))) <= cutoff_days) %>% 
                        ungroup()

        if (slice_by == "most_recent") { # most recent values 
                cohort <- cohort %>% 
                        group_by({{grouping_id}}, !!sym(paste0(value_type, "_type"))) %>% 
                        slice_min(abs(!!sym(paste0(value_type, "_days_since_transplant"))), n = no_value, with_ties = keep_tie) %>% 
                        mutate(!!sym(paste0(value_type, "_", as.character(no_value))) := mean({{value_name}}, na.rm = TRUE)) %>% 
                        mutate(!!sym(paste0(value_type, "_days_since_transplant_", as.character(no_value))) := mean(!!sym(paste0(value_type, "_days_since_transplant")), na.rm = TRUE)) %>% 
                        slice_min(abs(!!sym(paste0(value_type, "_days_since_transplant")) ), n = 1, with_ties = FALSE) %>% ungroup()   
                cohort <- cohort %>% val_classification(cohort = ., 
                                                        val = {{value_name}}, 
                                                        label = !!sym(paste0(value_type, "_level")), 
                                                        cutoffs = cutoffs) %>%
                                    val_classification(cohort = ., 
                                                        val = !!sym(paste0(value_type, "_", as.character(no_value))), 
                                                        label = !!sym(paste0(value_type, "_", as.character(no_value), "_level")), 
                                                        cutoffs = cutoffs)
        } else if( slice_by == "max") {
                cohort <- cohort %>% 
                        group_by({{grouping_id}}, !!sym(paste0(value_type, "_type"))) %>% 
                        slice_max({{value_name}}, n = no_value, with_ties = keep_tie) %>% 
                        mutate(!!sym(paste0(value_type, "_", as.character(no_value), "_max")) := mean({{value_name}}, na.rm = TRUE)) %>% 
                        mutate(!!sym(paste0(value_type, "_days_since_transplant_", as.character(no_value), "_max")) := mean(!!sym(paste0(value_type, "_days_since_transplant")), na.rm = TRUE)) %>% 
                        slice_max({{value_name}}, n = 1, with_ties = FALSE) %>% ungroup() %>%
                        rename(!!sym(paste0(value_type, "_max")) := {{value_name}}, 
                            !!sym(paste0(value_type, "_days_since_transplant_max")) := !!sym(paste0(value_type, "_days_since_transplant")),
                            !!sym(paste0(value_type, "_type_max")) := !!sym(paste0(value_type, "_type")))   
                cohort <- cohort %>% val_classification(cohort = ., 
                                                        val = !!sym(paste0(value_type, "_", as.character(no_value), "_max")), 
                                                        label = !!sym(paste0(value_type, "_", as.character(no_value), "_max_level")), 
                                                        cutoffs = cutoffs) %>%
                                    val_classification(cohort = ., 
                                                        val = !!sym(paste0(value_type, "_max")), 
                                                        label = !!sym(paste0(value_type, "_max_level")), 
                                                        cutoffs = cutoffs)
        } 

        return(cohort)
        
}



# flows data
# data <- ferritin %>% group_by(ferritin_pre, ferritin_pre_max) %>% 
#             summarise(n = n_distinct(person_id)) %>% ungroup() %>%
#             filter(!is.na(ferritin_pre), !is.na(ferritin_pre_max)) 

# # From these flows we need to create a node data frame: it lists every entities involved in the flow
# nodes <- data.frame(
#   name=c(as.character(data$ferritin_pre), 
#         as.character(data$ferritin_pre_max)) %>% unique()
# )

# # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
# data$IDsource <- match(data$ferritin_pre, nodes$name)-1 
# data$IDtarget <- match(data$ferritin_pre_max, nodes$name)-1

# # Make the Network
# p <- sankeyNetwork(Links = data, Nodes = nodes,
#               Source = "IDsource", Target = "IDtarget",
#               Value = "n", NodeID = "name", 
#               sinksRight=FALSE)
# print(p)



# Linear Interpolation Function to get ferritin values at a specific time point
# given the ferritin values measured on ferritin_dates
# compute the ferritin value at start_date and end_date for treatment

linear_interpolate_per_patient <- function(df, df_interpolated, day_window = 60) {
  # Ensure the data frame has the correct column names
  if (!all(c("ferritin_date", "ferritin") %in% colnames(df))) {
    stop("The data frame must have columns named 'ferritin_date' and ''ferritin'.")
  }
  
  # Sort the data frame by time (just in case it's not sorted)
  df <- df %>% group_by(person_id) %>% arrange(ferritin_date) %>% ungroup()

  # Perform interpolation for each patient and query time
  result <- df_interpolated
  result$interpolated_ferritin <- NA  # Add a column for the results
  result$ferritin0 <- NA
  result$ferritin1 <- NA

  result$t0 <- NA
  result$t1 <- NA

  for (i in 1:dim(df_interpolated)[1]) {
    person_id <- df_interpolated$person_id[i]
    t <- df_interpolated$interpolated_t[i]
    
    # Filter the data for the current patient
    patient_data <- df[df$person_id == person_id, ]
    
    # Check if query_time is within the range of the patient's time column
    if (nrow(patient_data) == 0 || 
        t < min(patient_data$ferritin_date) || 
        t > max(patient_data$ferritin_date)) {
      warning(paste("Query time for person_id", person_id, "is outside the range or data is missing."))
      # if the value is within 30 days, use it
        if(t >= min(patient_data$ferritin_date) - day_window && t <= min(patient_data$ferritin_date)){
          lower = min(patient_data$ferritin_date)
          result$interpolated_ferritin[i] <- patient_data$ferritin[patient_data$ferritin_date == lower]
          result$ferritin0[i] <- result$interpolated_ferritin[i]
          result$ferritin1[i] <- result$interpolated_ferritin[i]
          result$t0[i] <- lower
          result$t1[i] <- lower
        } else if (t <= max(patient_data$ferritin_date) + day_window && t >= max(patient_data$ferritin_date)){
        upper = max(patient_data$ferritin_date)
        result$interpolated_ferritin[i] <- patient_data$ferritin[patient_data$ferritin_date == upper]
        result$ferritin0[i] <- result$interpolated_ferritin[i]
        result$ferritin1[i] <- result$interpolated_ferritin[i]
        result$t0[i] <- upper
        result$t1[i] <- upper
        }
#       result[i, c("interpolated_ferritin", "ferritin0", "t0", "t1", "ferritin1")] <- rep(NA, 5)
      next
    }
    
    # Find the two nearest time points for interpolation
    lower <- max(patient_data$ferritin_date[patient_data$ferritin_date <= t])
    upper <- min(patient_data$ferritin_date[patient_data$ferritin_date >= t])
    
    result$t0[i] <- lower
    result$t1[i] <- upper

    # If the query time exactly matches a time point, return the corresponding cell count
    if (lower == upper) {
      result$interpolated_ferritin[i] <- patient_data$ferritin[patient_data$ferritin_date == t]
      result$ferritin0[i] <- result$interpolated_ferritin[i]
      result$ferritin1[i] <- result$interpolated_ferritin[i]
      next
    }
    
    # Get the corresponding cell counts for lower and upper times
    lower_count <- patient_data$ferritin[patient_data$ferritin_date == lower]
    upper_count <- patient_data$ferritin[patient_data$ferritin_date == upper]
    
    # Perform linear interpolation
    interpolated_ferritin <- lower_count + 
      (upper_count - lower_count) * (t - lower) / (upper - lower)
    
    result$interpolated_ferritin[i] <- interpolated_ferritin
    result$ferritin0[i] <- lower_count
    result$ferritin1[i] <- upper_count
  }
  
  return(result)

}

