# First, source all files in dir to get bigger named functions
functions <- list.files(
  path = "R",
  pattern = "*.R",
  recursive = TRUE)

sapply(paste0("R/", functions[functions != "functions.R"]), source, .GlobalEnv)
rm(functions)

# Define col pal
species_colors <- ggsci::pal_aaas("default")(5)
names(species_colors) <- c("Human (adult)", "Chimpanzee", "Gorilla", "Orangutan", "Human (infant)")
col_scale <- ggplot2::scale_color_manual(name = "Species", values = species_colors)
fill_scale <- ggplot2::scale_fill_manual(name = "Species", values = species_colors)

AOI_col_scale <- ggplot2::scale_color_manual(name = "resp", values = c(Other = "#107D81", Patient = "#F9A630", Agent = "#b22222"), breaks = c('Agent', 'Patient', 'Other')) 
AOI_fill_scale <- ggplot2::scale_fill_manual(name = "resp", values = c(Other = "#107D81", Patient = "#F9A630", Agent = "#b22222"), breaks = c('Agent', 'Patient', 'Other')) 

# Get relevant dataframe from RData file
extract_raw_data <- function(RData_file, footage_file, infants_file, outcome, bin_size = 5) {
  load(RData_file,  temp_env <- new.env())
  
  d_footage <- readxl::read_xlsx(footage_file) %>% 
    mutate(
      food_clip = case_when(
      is.na(`Food clip`) ~ "not food",
      `Food clip` == "Y" ~ "food"
    ),
      tool_use = case_when(
        `Tool use` == "NA" ~ "not tool use",
        `Tool use` == "N" ~ "not tool use",
        `Tool use` == "Y" ~ "tool use"
      ),
    camera_facing = case_when(
      `Camera facing` == "A" ~ "agent only",
      `Camera facing` == "both" ~ "agent and patient",
      `Camera facing` == "N" ~ "neither",
      `Camera facing` == "P" ~ "patient only"
     ),
    gaze_at_camera = case_when(
      `Gaze at camera` == "A" ~ "agent only",
      `Gaze at camera` == "N" ~ "neither",
      `Gaze at camera` == "P" ~ "patient only"
    )
    ) %>% 
    rename(stimulus = `Export file name`) %>% 
    select(stimulus, food_clip, tool_use, camera_facing, gaze_at_camera)
  
  d_infants <- readxl::read_xlsx(infants_file) %>% 
    mutate(Participant.name = paste0("Human baby_Participant", `Participant Number`)) %>% 
    rename(age_days = `Age in Days`) %>% 
    select(Participant.name, age_days) %>% 
    group_by(Participant.name) %>% 
    summarise(age_days = mean(age_days))
  
  env_list <- as.list(temp_env)
  
  if (outcome == "time") { 
    
    if (bin_size == 2.5) d <- env_list$data_timecourse_normalized_2.5percent
    if (bin_size == 5) d <- env_list$data_timecourse_normalized_5percent
    
    d$Stimulus <- sapply( str_split(d$Stimulus, "_"), "[", 1 )
    
    # Link up to footage data and infant participant data
    d <- d %>% 
      left_join(d_footage, by = c("Stimulus" = "stimulus")) %>% 
      left_join(d_infants)
    
    d <- d %>% 
      mutate(log_agent_size_z = scale(log(AOI_SIZE_AGENT)),
             log_patient_size_z = scale(log(AOI_SIZE_PATIENT)),
             log_AP_size = log(AOI_SIZE_AGENT / AOI_SIZE_PATIENT),
             AP_movement_diff = `Agent-Patient movement difference`,
             time_c = round(TIME_normalized - A_Action_StartTime_normalized, 3)
    ) %>% 
  mutate(stimulus_cat = case_when(
    food_clip == "not food" & Category == "inanimate" ~ 1,
    food_clip == "food" & Category == "inanimate" ~ 2,
    Category == "social" ~ 3
  ))
    
    # Drop observations before action starts, only relevant for a subset of videos
    d <- d %>% 
      filter(time_c >= 0)
  }
  
  return(d)
}

# Time-averaged data #####
time_data_avg <- function(data) {
  
  d <- data$d_time
  
  d_avg <- d %>% 
    group_by(TrialID) %>% 
    summarise(
      Agent_Sum = sum(Agent_Sum),
      Patient_Sum = sum(Patient_Sum),
      Other_Sum = sum(Other_Sum),
      Participant.name = unique(Participant.name),
      Species = unique(Species),
      FootageSpecies = unique(FootageSpecies),
      Stimulus = unique(Stimulus),
      log_AP_size = mean(log_AP_size),
      AP_md = mean(AP_md),
      stimulus_cat = unique(stimulus_cat)
    )
  
  return(d_avg)
  
}


# Data to export for dashboard #####
dashboard_data <- function(data){
  d <- data
  
  d <- d %>% 
    mutate(
      prop_agent = d$Agent_Sum / (d$Agent_Sum + d$Patient_Sum + d$Other_Sum),
      prop_patient = d$Patient_Sum / (d$Agent_Sum + d$Patient_Sum + d$Other_Sum),
      prop_other = d$Other_Sum / (d$Agent_Sum + d$Patient_Sum + d$Other_Sum),
      inanimate_social = d$Category
      ) %>% 
    select(TrialID, Participant.name, Species, Stimulus, FootageSpecies, inanimate_social, food_clip, tool_use, camera_facing, gaze_at_camera, age_days, log_AP_size, AP_movement_diff, prop_agent, prop_patient, prop_other, time_c)
  
  return(d)
}




# Additional processing for time series model, stan indices etc.
time_data <- function(data) {
  
  d <- data
  
  d$y <- with(d, cbind(Agent_Sum, Patient_Sum, Other_Sum))
  d$N_trials <- with(d, Agent_Sum + Patient_Sum + Other_Sum)
  d$y_prop <- d$y / d$N_trials
  
  d$AP_md <- ifelse(d$`Agent-Patient movement difference` == "A", 1, 0)
  
  N_obs <- nrow(d)
  
  ### Species indices #######
  species_id <- match(d$Species, unique(d$Species))
  N_species <- max(species_id)
  
  ### Footage species indices ######
  fs_id <- match(d$FootageSpecies, unique(d$FootageSpecies))
  N_fs <- max(fs_id)
  
  ### PID indices #########
  pid <- match(d$Participant.name, unique(d$Participant.name))
  N_pid <- max(pid)
  
  ### Time indices ########
  unique_times <- sort(unique(d$time_c))
  
  N_times <- length(unique_times)
  
  time_dist_mat <- as.matrix(dist(unique_times))
  time_id <- match(d$time_c, unique_times)
  
  ### Stimulus indices ####
  stimulus_id <- match(d$Stimulus, unique(d$Stimulus))
  N_stimuli <- max(stimulus_id)
  
  ### Trial indices ####
  trial_id <- match(d$TrialID, unique(d$TrialID))
  N_trials <- max(trial_id)
  
  ###########################
  data_list <- list(
    # indices and housekeeping for Stan
    y_prop = d$y_prop,
    N_obs = N_obs,
    N_species = N_species,
    N_fs = N_fs,
    N_pid = N_pid,
    N_stimuli = N_stimuli,
    N_trials = N_trials,
    K = 3, # number of categories
    species_id = species_id,
    fs_id = fs_id,
    pid = pid,
    stimulus_id = stimulus_id,
    trial_id = trial_id,
    
    # time indices
    N_times = N_times,
    time_dist_mat = time_dist_mat,
    time_id = time_id,
    time = d$time_c,
    
    # predictors
    stimulus_cat = d$stimulus_cat,
    AP_md = d$AP_md,
    log_AP_size = d$log_AP_size
  )
  
  return(list(d_time = d, stan_data_time = data_list))
}

fit_time_series <- function(data, n_iter, n_chains) {
    
    m_time <- rstan::stan_model("stan/m_time_GP.stan")
    data_list <- data$stan_data_time
  
    fit <- rstan::sampling(
      m_time,
      data = data_list,
      iter = n_iter,
      chains = n_chains,
      cores = n_chains,
      # exclude these parameters from returned fit
      pars = c("trial_b0_z", "time_trial_z", "time_trial_v", "lp__"),
      include = F,
      refresh = 1
    )

    return(fit)
}


GP_pred <- function(x_train, y_train, x_test, alpha) {
  # x_train: training data inputs
  # y_train: training data outputs
  # x_test: test data inputs
  
  # calculate kernel matrix for training data
  K <- (exp(-alpha * as.matrix(proxy::dist(x_train))) + (diag(length(x_train)) * 0.00001))
  K_inv <- chol2inv(chol(K))
  
  # calculate means for test data
  k_star <- exp(-alpha * as.matrix(proxy::dist(x_test, x_train)))
  mu_star <- k_star %*% K_inv %*% y_train
  
  # return the mean of the predicted values
  return(mu_star)
}


timeseries_pred <- function(
    post,
    time_train,
    time_test,
    n_samps = length(post$lp__),
    species_RE = F,
    pid_RE = F,
    stimulus_RE = F,
    fs_RE = F,
    species_id = 1,
    pid = 1,
    stimulus_id = 1,
    fs_id = 1,
    stimulus_cat = 3,
    log_AP_size = 0,
    AP_md = 0
) {
  
  post2 <- post
  
  # Generate predictions for new time smooths, add to post2
  post2$time_pred <- array(0, dim = c(n_samps, 2, length(time_test)))
  post2$time_species_pred <- array(0, dim = c(n_samps, 2, length(time_test)))
  post2$time_pid_pred <- array(0, dim = c(n_samps, 2, length(time_test)))
  post2$time_stimulus_pred <- array(0, dim = c(n_samps, 2, length(time_test)))
  
  for (k in 1:2) {
    for (i in 1:n_samps) {
      ## Make prediction for time Gaussian process ####
      post2$time_pred[i, k, ] <- GP_pred(x_train = time_train, y_train = post2$time_v[i, stimulus_cat, ,k], x_test = time_test, alpha = post2$alpha_time[i, stimulus_cat, k])[,1]
      
      post2$time_species_pred[i, k, ] <- GP_pred(x_train = time_train, y_train = post$time_species_v[i, species_id, stimulus_cat, ,k], x_test = time_test, alpha = post2$alpha_time_species[i, stimulus_cat, k])[,1]
      
      post2$time_pid_pred[i, k, ] <- GP_pred(x_train = time_train, y_train = post2$time_pid_v[i, pid, stimulus_cat, ,k], x_test = time_test, alpha = post2$alpha_time_pid[i, stimulus_cat, k])[,1]
      
      post2$time_stimulus_pred[i, k, ] <- GP_pred(x_train = time_train, y_train = post2$time_stimulus_v[i, stimulus_id, ,k], x_test = time_test, alpha = post2$alpha_time_stimulus[i, k])[,1]
    }
  }
  
  preds <- array(0, dim = c(n_samps, 3, length(time_test)))
  
  for (time in 1:length(time_test)) {
    theta <- matrix(NA, nrow = n_samps, ncol = 3)
    
    for (k in 1:2) {

      intercept_stem = with(post2, 
                            b0[,stimulus_cat, k] + 
                              species_RE * species_b0_z[,species_id, stimulus_cat, k] * sqrt(sigmasq_b0[,stimulus_cat, k] * phi_b0[,stimulus_cat, k, 1]) + 
        pid_RE * pid_b0_z[, pid, stimulus_cat, k] * sqrt(sigmasq_b0[, stimulus_cat, k] * phi_b0[, stimulus_cat, k, 2]) + 
        stimulus_RE * stimulus_b0_z[, stimulus_id,k] * sqrt(sigmasq_b0[, stimulus_cat, k] * phi_b0[, stimulus_cat, k, 3]) + 
        fs_RE * fs_b0_z[, fs_id, k] * sqrt(sigmasq_b0[, stimulus_cat, k] * phi_b0[, stimulus_cat, k, 5]) + 
        fs_RE * species_RE * fs_species_b0_z[, fs_id, species_id, k] * sqrt(sigmasq_b0[, stimulus_cat, k] * phi_b0[, stimulus_cat, k, 6]) +
        
      (b_AP_size[, stimulus_cat, k] + 
         species_RE * species_AP_size_z[, species_id, stimulus_cat, k] * sqrt(sigmasq_AOI_size[, stimulus_cat, k] * phi_AOI_size[, stimulus_cat, k, 1]) + 
          pid_RE * pid_AP_size_z[, pid, stimulus_cat, k] * sqrt(sigmasq_AOI_size[, stimulus_cat, k] * phi_AOI_size[, stimulus_cat, k, 2])) * log_AP_size +
        
        (b_AP_movement[, stimulus_cat, k] + 
           species_RE * species_AP_movement_z[, species_id, stimulus_cat, k] * sqrt(sigmasq_AP_md[, stimulus_cat, k] * phi_AP_md[, stimulus_cat, k, 1]) + 
           pid_RE * pid_AP_movement_z[, pid, stimulus_cat, k] * sqrt(sigmasq_AP_md[, stimulus_cat, k] * phi_AP_md[, stimulus_cat, k, 2])) * AP_md
      )
      
      # GP smooths for time
      time_stem = with(post2, time_pred[, k, time] * sqrt(sigmasq_time_mu[ ,stimulus_cat, k]) + 
        species_RE * time_species_pred[,k, time] * sqrt(sigmasq_time[ ,stimulus_cat, k] * phi_time[ , stimulus_cat, k, 1]) + 
        pid_RE * time_pid_pred[,k, time] * sqrt(sigmasq_time[ ,stimulus_cat, k] * phi_time[ ,stimulus_cat, k, 2]) + 
        stimulus_RE * time_stimulus_pred[,k, time] * sqrt(sigmasq_time[ ,stimulus_cat, k] * phi_time[ ,stimulus_cat, k, 3])
      )
      
      theta[,k] = intercept_stem + time_stem;
    }
    theta[,3] = 0; # ref category (other)
  
    # Apply softmax to get probabilities
    post_epreds <- t(apply(theta, 1, softmax))
    preds[ , , time] <- post_epreds
}

return(preds)
}
