# survival analysis
# source("code/EDA.R")
# install.packages(c("survival", "survminer"))
library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
require(ggplot2)
require(lubridate)
# install.packages("ggsurvfit")
# the wasy survival is coded here is different from table 1
# this is overall survival, death of any causes 

cohort_covars <- results_tbl("analytics_dataset") %>% collect() %>% 
                    mutate(age_at_ce = as.numeric(difftime(transplant_date, birth_date, units = "days"))/365.25, 
                            ethnicity = factor(ethnicity, levels = c("Hispanic or Latino", "Not Hispanic or Latino", "Other", "Unknown"))) %>%
          mutate(survival = if_else(survival, 1, 0), 
                survival_1yr = if_else(death_days_since_ce >= 365.25 | is.na(death_days_since_ce), 0, 1),
                survival_3yr = if_else(death_days_since_ce >= (365.25*3) | is.na(death_days_since_ce), 0, 1),
                survival_5yr = if_else(death_days_since_ce >= (365.25*5) | is.na(death_days_since_ce), 0, 1),
                graft_fail_days_since_ce = as.numeric(difftime(graft_fail_date, transplant_date, units = "days"))) %>%
          mutate(age_at_ce = round(as.numeric(difftime(transplant_date, birth_date, unit = "days"))/365.25), 
         donor_relation = if_else(is.na(donor_relation), "Unknown", donor_relation),
         match_status = if_else(is.na(match_status), "Unknown", match_status),
         second_transplant = if_else(is.na(second_transplant_date), "No", "Yes"),
         gender = if_else(gender == "FEMALE", "Female", "Male"),
         scd_type = recode(scd_type, "aa" = "aplastic anemia", "dba" = "Diamond-Blackfan anemia", "bta" = "beta thalassemia major", "scd" = "sickle cell disease")) %>%
        mutate(CD3_reconstitution = case_when(CD3_reconstitute_3mon ~ "0-3 months",
                                        CD3_reconstitute_6mon ~ "3-6 months",
                                        CD3_reconstitute_9mon ~ "6-9 months",
                                        CD3_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         CD4_reconstitution = case_when(CD4_reconstitute_3mon ~ "0-3 months",
                                        CD4_reconstitute_6mon ~ "3-6 months",
                                        CD4_reconstitute_9mon ~ "6-9 months",
                                        CD4_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         CD8_reconstitution = case_when(CD8_reconstitute_3mon ~ "0-3 months",
                                        CD8_reconstitute_6mon ~ "3-6 months",
                                        CD8_reconstitute_9mon ~ "6-9 months",
                                        CD8_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months"),
         IgM_reconstitution = case_when(IgM_reconstitute_3mon~ "0-3 months",
                                        IgM_reconstitute_6mon ~ "3-6 months",
                                        IgM_reconstitute_3mon ~ "6-9 months",
                                        IgM_reconstitute_12mon ~ "9-12 months", TRUE ~ "Never/after 12 months")) %>%
         mutate(across(ends_with("reconstitution"), ~factor(.x, levels = c("0-3 months", "3-6 months", "6-9 months", "9-12 months", "Never/after 12 months")))) %>% 
         mutate(vod_subgroup = if_else(vod, scd_type, "No VOD")) %>% 
         mutate(vod_subgroup = factor(vod_subgroup, levels = c("beta thalassemia major", "sickle cell disease", "No VOD")))%>%
         mutate(busulfan_and_vod = case_when(has_busulfan & vod ~ "busulfan & VOD",
                                            has_busulfan & !vod ~ "busulfan & No VOD", 
                                            TRUE ~ "Other categories")) %>%
         mutate(across(c("has_busulfan"), ~if_else(is.na(.x) | !(.x), "No", "Yes")))

cohort_covars <- cohort_covars %>% 
  mutate(#survival = if_else(survival, "Alive", "Death"), 
         #censoring status 0=censored, 1=dead
         survival_time_1yr_censored = case_when(death_days_since_ce <= 365.25 ~ death_days_since_ce, 
                                                death_days_since_ce > 365.25 ~ 365.25,
                                                is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25, na.rm = TRUE),
                                                TRUE ~ NA_real_),
         survival_time_3yr_censored = case_when(death_days_since_ce <= 365.25*3 ~ death_days_since_ce, 
                                                death_days_since_ce > 365.25*3 ~ 365.25*3,
                                                is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*3, na.rm = TRUE),
                                                TRUE ~ NA_real_),
        survival_time_5yr_censored = case_when(death_days_since_ce <= 365.25*5 ~ death_days_since_ce, 
                                               death_days_since_ce > 365.25*5 ~ 365.25*5, 
                                               is.na(death_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25*5, na.rm = TRUE),
                                               TRUE ~ NA_real_),
         graft_fail_days_since_ce_5yr_censored = case_when(graft_fail_days_since_ce > 365.25*5 ~ 365.25*5, #censored at 5 years
                                                           graft_fail_days_since_ce <= 365.25*5 ~ graft_fail_days_since_ce, #graft failure
                                                           death_days_since_ce <= 365.25*5 ~ death_days_since_ce, #death: competing event
                                                           last_in_person_visit_since_ce <= 365.25*5 & is.na(graft_fail_days_since_ce) ~ last_in_person_visit_since_ce, #censored at last visit
                                                           last_in_person_visit_since_ce > 365.25*5 & is.na(graft_fail_days_since_ce) ~ 365.25*5, #censored at 5 years follow-up
                                                            TRUE ~ NA),
         # graft failure is a competing event: 0: censored, 1: graft failure, 2:death                                       
         graft_fail_5yr = case_when(graft_fail_days_since_ce <= 365.25*5 ~ "1",
                                    death_days_since_ce <= 365.25*5 ~ "2", 
                                    graft_fail_days_since_ce > 365.25*5 | !is.na(last_in_person_visit_since_ce) ~ "0",
                                    TRUE ~ NA))

cohort_covars <- cohort_covars %>% #mutate(ferritin_pre = factor(ferritin_pre, levels = c("low", "moderate", "high")),
         #ferritin_pre_max = factor(ferritin_pre_max, levels = c("low", "moderate", "high"))) %>%
  mutate(bacteremia_time_1yr_censored = case_when(bacteremia_days_since_ce > 365.25 ~ 365.25, 
                                                  bacteremia_days_since_ce >= 0 & bacteremia_days_since_ce <= 365.25 ~ bacteremia_days_since_ce,death_days_since_ce <= 365.25 ~ death_days_since_ce, #if die before 5 years, mark as censored 
                                                  is.na(bacteremia_days_since_ce) & !is.na(last_in_person_visit_since_ce) ~ pmin(last_in_person_visit_since_ce, 365.25, na.rm = TRUE),
                                               TRUE ~ NA_real_)) %>%
  mutate(bacteremia_1yr = case_when(is.na(bacteremia_days_since_ce) & bacteremia_time_1yr_censored < 365.25 ~ 0, 
                                bacteremia_time_1yr_censored >= 0 & bacteremia_time_1yr_censored < 365.25 & !is.na(bacteremia_days_since_ce)~ 1, 
                                bacteremia_time_1yr_censored >= (365.25) ~0)) %>%
  mutate(graft_manip = if_else(is.na(graft_manip), "Unknown", graft_manip)) %>%
  mutate(disease_relapse = if_else(is.na(graft_fail_date), "No", "Yes")) 

cohort_covars %>% output_tbl("survival_analysis_dataset")

# results_tbl("cr_data") %>% collect() %>% output_tbl("cr_data", local = TRUE, file = TRUE)
                              
# stat_dataset %>% filter(chart_completion, eligibility == "Yes") %>%
#       group_by(survival_5yr, transplant_type_combined) %>% summarise(n = n_distinct(record_id)) %>% view()

# now get ferritin time series for aim 3
# there are patients who only had 1 ferritin measurement after transplant
# rslt$potential_aim3 <- results_tbl("analytics_dataset") %>%
#       filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date, IRT) %>%
#       left_join(results_tbl("covar_ferritin_mx"), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>% 
#       filter(transplant_date - ferritin_date <= 365) %>%
#       group_by(person_id) %>% 
#       mutate(n = n()) %>% filter(n >=1) %>% ungroup()

# # we want the base ferritin levels to be closet to the start date of treatment
# rslt$phleb_px <- results_tbl("analytics_dataset") %>%
#       # filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date) %>%
#       find_procedures(procedure_codeset_name = "therapeutic_phleb_px") %>%
#       filter(procedure_date >= transplant_date) %>%
#       filter(procedure_date - transplant_date <= 365) %>%
#       group_by(person_id) %>%
#       collect_new() %>%
#       mutate(start_date = min(procedure_date, na.rm = TRUE), 
#              end_date = max(procedure_date, na.rm = TRUE),
#              duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration) %>% 
#       mutate(route = NA) %>%
#       mutate(IRT_type = "phlebotomy") 

# # patients who got chelation
# dx_codeset<- load_codeset("defibrotide_rx") %>% mutate(type = "defibrotide") %>% #defibrotide
#                 union(load_codeset("deferoxamine_rx") %>% mutate(type = "deferoxamine")) %>% #deferoxamine
#                 union(load_codeset("deferasirox_rx") %>% mutate(type = "deferasirox")) %>% #deferasirox
#                 compute_new(temp = TRUE, name = "drug_id")

# rslt$chelation_rx <- results_tbl("analytics_dataset") %>%
#       # filter(chart_completion, eligibility == "Yes") %>%
#       select(record_id, person_id, transplant_date) %>%
#       find_drugs(dx_codeset) %>%
#       filter(drug_exposure_start_date >= transplant_date) %>%
#       filter(drug_exposure_start_date - transplant_date <= 365) %>% 
#       distinct(person_id, drug_exposure_start_date, IRT_type = type, .keep_all = TRUE) 

# # for intravenous drugs, assume a single one per day
# rslt$chelation_iv <- rslt$chelation_rx %>% filter(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous")) %>%
#       group_by(person_id) %>%
#       collect_new() %>%
#       mutate(start_date = min(drug_exposure_start_date, na.rm = TRUE), 
#              end_date = max(drug_exposure_start_date, na.rm = TRUE),
#              duration = as.numeric(end_date - start_date)) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration, .keep_all = TRUE) %>% 
#       mutate(route = "IV") %>%
#       select(person_id, record_id, start_date, end_date, transplant_date, duration, IRT_type, route)

# # oral drugs
# rslt$chelation_oral <- rslt$chelation_rx %>% filter(!(route_source_value %in% c("Intravenous", "IV", "Injection", "Subcutaneous"))) %>%
#       filter(!is.na(quantity)) %>%
#       collect_new() %>%
#       group_by(person_id) %>%
#       mutate(start_date = min(drug_exposure_start_date, na.rm = TRUE), 
#              end_date = drug_exposure_start_date + days(quantity),
#              duration = quantity) %>% ungroup() %>%
#       distinct(person_id, record_id, start_date, end_date, transplant_date, duration, .keep_all = TRUE) %>% 
#       select(person_id, record_id, start_date, end_date, transplant_date, duration, IRT_type, route = route_source_value) 

# # Only took the duration of the first prescription, some prescriptions have up to 11 refills 

# # append the chelation and IRT patients
# rslt$irt <- rslt$phleb_px %>%
#           union(rslt$chelation_iv) %>%
#           union(rslt$chelation_oral) 

# # who got more than 1 treatment
# rslt$irt_exclude <- rslt$irt %>% group_by(record_id) %>% 
#       summarise(n = n()) %>% 
#       filter(n> 1) %>% ungroup() %>% view()
      
# rslt$irt_exclude %>% 
#       inner_join(rslt$irt, by = "record_id") %>% 
#       arrange(person_id) %>% view()
      
# # exclude patients with more than 1 treatment for now: 
# rslt$irt <- rslt$irt %>% anti_join(rslt$irt_exclude, by = "record_id")

# rslt$irt %>% distinct_ct()

# # ferritin_0 is the one closest to the start date of treatment
# rslt$ferritin_start <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>%
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       # filter(transplant_date - ferritin_date <= 365) %>%
#       group_by(person_id) %>% 
#       filter(start_date <= ferritin_date) %>%
#       slice_max(start_date, n =1, with_ties = FALSE) %>% ungroup() %>%
#       select(person_id, record_id, transplant_date, start_date, ferritin0_date = ferritin_date, ferritin0 = ferritin) 
      

# rslt$ferritin_end <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>%
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>% 
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       filter(end_date <= ferritin_date) %>%
#       group_by(person_id) %>% 
#       slice_min(end_date, n =1, with_ties = FALSE) %>% ungroup() %>%
#       select(person_id, record_id, transplant_date, end_date, ferritinf_date = ferritin_date, ferritinf = ferritin)

# rslt$ferritin_all <- rslt$irt %>% copy_to_new(df = ., name = "irt") %>% 
#       inner_join(results_tbl("covar_ferritin_mx") %>%
#                   select(person_id, ferritin_date, ferritin), by = "person_id") %>%
#       filter(transplant_date <= ferritin_date) %>%
#       filter(!is.na(ferritin)) %>%
#       filter(transplant_date - ferritin_date <= 365) %>% collect_new()

# print(rslt$ferritin_all %>% 
#       mutate(ferritin_date = as.numeric(ferritin_date - transplant_date)) %>%
#       mutate(start_date = as.numeric(start_date - transplant_date)) %>%
#       mutate(end_date = as.numeric(end_date - transplant_date)) %>% 
#       filter(duration > 0) %>%
#       ggplot(aes(x = ferritin_date, y = ferritin, group = record_id)) + 
#       geom_point() + 
#       geom_line() +
#       geom_vline(aes(xintercept = start_date), linetype = "dashed", color = "red") +
#       geom_vline(aes(xintercept = end_date), linetype = "dashed", color = "red") +
#       facet_wrap(~record_id, ncol = 3)) 

# # special cases: patients with only 1 treatment: duration == 0
# rslt$ferritin_start %>% inner_join(rslt$ferritin_end, by = c("record_id", "person_id", "transplant_date")) %>% 
#             mutate(duration = end_date - start_date) %>% view()


# graft failure competing risk model
pre_ferritin_cox_2 <- pre_ferritin_cox %>% 
                          filter(transplant_type_combined != "Unknown") %>% 
                          # filter(transplant_type_combined != "autologous/genetherapy") %>% 
                          filter(!is.na(graft_fail_5yr)) %>% 
                          # mutate(graft_fail_5yr = recode(graft_fail_5yr, "0" = "event-free", "1" = "graft failure", "2" = "death")) %>%
                          # mutate(graft_fail_5yr = factor(graft_fail_5yr,
                                                  # levels = c("event-free", "graft failure", "death"))) %>% 
                          filter(graft_fail_days_since_ce_5yr_censored >= 0) %>%
                          mutate(transplant_type_combined = relevel(factor(transplant_type_combined), ref = "matched related"))

# Prepare transition matrix
tmat <- trans.comprisk(2, names = c("event-free", "graft failure", "death"))

# Run msprep
pre_ferritin_cox_2$graft_fail <- as.numeric(pre_ferritin_cox_2$graft_fail_5yr == 1)
pre_ferritin_cox_2$death <- as.numeric(pre_ferritin_cox_2$graft_fail_5yr == 2)

pre_ferritin_cox_2_long <- msprep(
  time = c(NA, "graft_fail_days_since_ce_5yr_censored", "graft_fail_days_since_ce_5yr_censored"), 
  status = c(NA, "graft_fail", "death"), 
  data = pre_ferritin_cox_2, 
  keep = c("scd_type", "has_busulfan", "transplant_type_combined", "age_at_ce", "gender", "ferritin_pre"), 
  trans = tmat)

events(pre_ferritin_cox_2_long)

pre_ferritin_cox_2_long %>% view()

# Run cox model
pre_ferritin_cox_2_expanded <- expand.covs(pre_ferritin_cox_2_long, "ferritin_pre")
colnames(pre_ferritin_cox_2_expanded)
pre_ferritin_cox_2_expanded %>% view()
c1 <- coxph(Surv(time, status) ~ ferritin_premoderate.1 + ferritin_prehigh.2 +
                                  ferritin_prehigh.1 + ferritin_prehigh.2 + 
                                  strata(trans),
            data = pre_ferritin_cox_2_expanded)

# Visualising cumulative baseline hazards using plot.msfit()
# Data to predict
ferritin_low <- data.frame(
  ferritin_premoderate.1 = c(0, 0, 0, 0),
  ferritin_premoderate.2 = c(0, 0, 0, 0),
  ferritin_prehigh.1 = c(0, 0, 0, 0),
  ferritin_prehigh.2 = c(0, 0, 0, 0),
  trans = c(1, 2), 
  strata = c(1, 2)
)

# Make msfit object
msf.ferritin_low <- msfit(
  object = c1, 
  newdata = ferritin_low, 
  trans = tmat)

# Plot cumulative baseline hazards
# different line colors for different strata
plot(msf.ferritin_low, 
    type = "single",
    xlim = c(0, 5*365.25), 
    col = c("red", "blue"), 
    lty = c("dashed", "solid"), 
    legend.pos = "top",
    lwd = 2,
    use.ggplot = TRUE,
    # conf.type = "log", 
  # conf.int = 0.95
  ) + # Add title and center
  ggtitle("Cumulative baseline hazards") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot cumulative baseline hazards using facet wrap for event types

plot(msf.ferritin_low, 
    xlim = c(0, 5*365.25), 
    col = c("red", "blue"), 
    lty = c("dashed", "solid"), 
    legend.pos = "top",
    lwd = 2,
    use.ggplot = TRUE,
    type = "separate", #for 2 boxes
    scale_type = "fixed"

  ) + # Add title and center
  ggtitle("Cumulative baseline hazards") +
  theme(plot.title = element_text(hjust = 0.5))

# Filled/stacked plots
# patient-specific transition probabilities using probtrans(). 
# we would like to predict from the beginning of follow-up (predt = 0).
pt.ferritin <- probtrans(msf.ferritin_low, predt = 0)

# Example predict at different times
summary(pt.ferritin, times = c(1, 5, 10))

# stacked plot 
plot(pt.ferritin, from = 1)

# reordering the event types to have 1 competing event at the bottom and 1 at the top

plot(pt.ferritin, 
    use.ggplot = TRUE, # return a ggplot object and put legend outside
    ord = c(2, 1, 3))

# each event type in one subplot
par(mfrow = c(1, 3))
plot(pt.ferritin, 
    type = "separate",
    conf.int = 0.95, # 95% level
  conf.type = "log")

# comparing different ferritin levels
# 1. Prepare patient data 
ferritin_low <- data.frame(
  ferritin_premoderate.1 = c(0, 0, 0, 0),
  ferritin_premoderate.2 = c(0, 0, 0, 0),
  ferritin_prehigh.1 = c(0, 0, 0, 0),
  ferritin_prehigh.2 = c(0, 0, 0, 0),
  trans = c(1, 2), 
  strata = c(1, 2)
)

ferritin_med <- data.frame(
  ferritin_premoderate.1 = c(1, 0, 0, 0),
  ferritin_premoderate.2 = c(0, 1, 0, 0),
  ferritin_prehigh.1 = c(0, 0, 0, 0),
  ferritin_prehigh.2 = c(0, 0, 0, 0),
  trans = c(1, 2), 
  strata = c(1, 2)
)

ferritin_high <- data.frame(
  ferritin_premoderate.1 = c(0, 0, 0, 0),
  ferritin_premoderate.2 = c(0, 0, 0, 0),
  ferritin_prehigh.1 = c(1, 0, 0, 0),
  ferritin_prehigh.2 = c(0, 1, 0, 0),
  trans = c(1, 2), 
  strata = c(1, 2)
)

# 2. Make msfit objects
msf.low <- msfit(c1, ferritin_low, trans = tmat)
msf.med <- msfit(c1, ferritin_med, trans = tmat)
msf.high <- msfit(c1, ferritin_high, trans = tmat)


# 3. Make probtrans objects
pt.low <- probtrans(msf.low, predt = 0)
pt.med <- probtrans(msf.med, predt = 0)
pt.high <- probtrans(msf.high, predt = 0)


vis.multiple.pt(
  x = list(pt.low, pt.med, pt.high), 
  from = 1,
  to = 2, 
  # conf.type = "log",
  cols = c("red", "green", "blue"),
  labels = c("low", "moderate", "high"),
  legend.title = "Pre-transplant ferritin",
  ylab = "Probability of graft failure"
)

vis.multiple.pt(
  x = list(pt.low, pt.med, pt.high), 
  from = 1,
  to = 3, 
  # conf.type = "log",
  cols = c("red", "green", "blue"),
  labels = c("low", "moderate", "high"),
  legend.title = "Pre-transplant ferritin",
  ylab = "Probability of death without graft failure"
)

## another way
pt.f_low <- probtrans(msf.f_low, predt = 0)[[1]]
idx1 <- (pt.f_low$time < 5)

pt.f_moderate <- probtrans(msf.f_moderate, predt = 0)[[1]]
idx2 <- (pt.f_moderate$time < 5)

plot(c(0, pt.f_low$time[idx1]), c(0, pt.f_low$pstate2[idx1]), type = "s",
        # type = "single",
        ylim = c(0, 0.5), xlab = "Years from transplant", ylab = "Probability of graft failure",
        lwd = 2, use.ggplot = TRUE) 
lines(c(0, pt.f_moderate$time[idx2]), c(0, pt.f_moderate$pstate2[idx2]), type = "s", 
      lwd = 2, col = 2)
title(main = "AIDS")