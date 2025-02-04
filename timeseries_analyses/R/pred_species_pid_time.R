pred_species_pid_time <- function(fit, data){
  
  d <- data$d_time
  data_list <- data$stan_data_time
  
  time_train <- sort(unique(d$time_c)) # input training points to GP
  time_test <- seq(from = min(d$time_c), to = max(d$time_c), length.out = 20) # time values to predict
  
  post <- extract.samples(fit)
  n_samps <- length(post$lp__)
  
  #### Make species average predictions #########
  d$species_id <- data_list$species_id
  
  # get labels to match indices
  d_species <- d %>% 
    group_by(species_id) %>% 
    summarise(Species = unique(Species))
  
  newdata_time <- expand.grid(
    species_id = 1:data_list$N_species,
    stimulus_cat = 1:3
  )
  
  newdata_time <- newdata_time %>% 
    mutate(
      condition = 1:n()
    )
  
  n_conditions <- max(newdata_time$condition)
  
  registerDoParallel(cores = 8)
  
  pred_time <- array(
    # loop over conditions
    foreach(j = 1:n_conditions, .combine = 'c') %dopar% {
      timeseries_pred(post = post,
                      time_train = time_train,
                      time_test = time_test,
                      n_samps = n_samps,
                      species_id = newdata_time$species_id[j],
                      species_RE = T,
                      stimulus_cat = newdata_time$stimulus_cat[j],
      )},
    # structure the array
    dim = c(n_samps, 3, length(time_test), n_conditions), dimnames = list(samps = 1:n_samps, resp = c("Agent", "Patient", "Other"), time_c = time_test, condition = 1:n_conditions))
  
  # Pivot array long
  pred_time_long <- pred_time %>% 
    cubelyr::as.tbl_cube(met_name = "est") %>% 
    as_tibble %>% 
    left_join(newdata_time) %>% 
    left_join(d_species)
  
  ### PID ###############################
  # get labels to match indices
  d$pid <- data_list$pid
  
  d_pid <- d %>% 
    group_by(pid) %>% 
    summarise(Participant.name = unique(Participant.name), Species = unique(Species), species_id = unique(species_id))
  
  newdata_time_pid <- expand.grid(
    pid = 1:data_list$N_pid,
    stimulus_cat = 1:3
  )
  
  newdata_time_pid <- left_join(newdata_time_pid, d_pid)
  
  newdata_time_pid <- newdata_time_pid %>% 
    mutate(
      condition = 1:n()
    )
  
  n_conditions <- max(newdata_time_pid$condition)
  
  registerDoParallel(cores = 1)
  
  pred_pid_time <- array(
    # loop over conditions
    foreach(j = 1:n_conditions, .combine = 'c') %dopar% {
      timeseries_pred(post = post,
                      time_train = time_train,
                      time_test = time_test,
                      pid = newdata_time_pid$pid[j],
                      species_id = newdata_time_pid$species_id[j],
                      species_RE = T,
                      pid_RE = T,
                      stimulus_cat = newdata_time_pid$stimulus_cat[j]
      )},
    # structure the array
    dim = c(n_samps, 3, length(time_test), n_conditions), dimnames = list(samps = 1:n_samps, resp = c("Agent", "Patient", "Other"), time_c = time_test, condition = 1:n_conditions))
  
  # Pivot array long
  pred_pid_time_long <- pred_pid_time %>% 
    cubelyr::as.tbl_cube(met_name = "est") %>% 
    as_tibble %>% 
    left_join(newdata_time_pid) %>% 
    left_join(d_pid)
  
  ######
  
  ## Change factor level order ###########
  pred_time_long$stimulus_cat <- case_when(pred_time_long$stimulus_cat == 1 ~ "inanimate (not food)",
                                           pred_time_long$stimulus_cat == 2 ~ "inanimate (food)",
                                           pred_time_long$stimulus_cat == 3 ~ "social",
                                           TRUE ~ as.character(pred_time_long$stimulus_cat))
  
  pred_pid_time_long$stimulus_cat <- case_when(pred_pid_time_long$stimulus_cat == 1 ~ "inanimate (not food)",
                                               pred_pid_time_long$stimulus_cat == 2 ~ "inanimate (food)",
                                               pred_pid_time_long$stimulus_cat == 3 ~ "social",
                                               TRUE ~ as.character(pred_pid_time_long$stimulus_cat))
  
  pred_time_long$resp <- factor(pred_time_long$resp, levels = (c("Agent", "Patient", "Other")))
  pred_pid_time_long$resp <- factor(pred_pid_time_long$resp, levels = (c("Agent", "Patient", "Other")))
  
  return(list(pred_time_long = pred_time_long, pred_pid_time_long = pred_pid_time_long))
}
