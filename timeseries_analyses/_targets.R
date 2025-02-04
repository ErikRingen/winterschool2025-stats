library(targets)
source("R/functions.R")

# Dependencies for this pipeline
tar_option_set(packages = c(
  "tidyverse",
  "rethinking",
  "rstan",
  "patchwork",
  "ggsci",
  "RColorBrewer",
  "cubelyr",
  "proxy",
  "foreach",
  "doParallel",
  "ggdist"
))

list(
  # Extract raw data
  tar_target(rdata_file,
             "data/input-to-models.RData", format = "file"),
  tar_target(footage_file,
             "data/Study1_footage_details.xlsx"),
  tar_target(infants_file,
             "data/Infant_age.xlsx"),
  tar_target(d_raw_time,
             extract_raw_data(rdata_file, footage_file, infants_file, "time")),
  
  ## Generate data for Shiny dashboard #####
  tar_target(d_dash,
             dashboard_data(d_raw_time)),

  # Data for time series model 
  tar_target(d_time,
             time_data(d_raw_time)),
  
  ## Generate data for time-averaged model comparison
  tar_target(d_avg,
             time_data_avg(d_time)),
  
  ###### Time Averaged Analyses #########################
  tar_target(model_compare,
             avg_model_comparison(d_avg)),
  
  ###### Time Series Analyses ##########################
  # Fit time series model
  tar_target(fit_time,
             fit_time_series(d_time, n_iter = 1000, n_chains = 8)),
  
  # Make predictions for time series
  tar_target(pred_time,
             pred_species_pid_time(fit_time, d_time)),
  
  # Make plots + export summaries of posterior predictions to .csv
  # time series plots
  tar_target(p_species_time_agent,
             plot_species_time_AOI(pred_time, "Agent")),
  
  tar_target(p_species_time_patient,
             plot_species_time_AOI(pred_time, "Patient")),
  
  tar_target(p_species_time_other,
             plot_species_time_AOI(pred_time, "Other")),
  
  tar_target(p_species_time_OR,
             plot_species_time_OR(pred_time)),
  
  tar_target(p_species_time_APO,
             plot_species_time_OR_APO(pred_time)),
  
  tar_target(p_species_time_AO,
             plot_species_time_OR_AO(pred_time)),
  
  # AOI size plots
  tar_target(p_species_AOI_size,
             plot_species_AOI_size(fit_time, d_time)),
  
  # AP movement diff
  tar_target(p_species_AP_md,
             plot_species_AP_md(fit_time, d_time))
  
)
