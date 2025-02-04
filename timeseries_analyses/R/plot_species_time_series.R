plot_species_time_series <- function(fit, data, bin_size = 5){
  
  d <- data$d_time
  data_list <- data$stan_data_time
  
  time_train <- sort(unique(d$time_c)) # input training points to GP
  time_test <- seq(from = min(d$time_c), to = max(d$time_c), length.out = 30) # time values to predict
  
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
  
  pred_time_long_summary <- pred_time_long %>% 
    group_by(condition, time_c, resp) %>% 
    summarise(med = median(est), Species = unique(Species), stimulus_cat = unique(stimulus_cat)) %>% 
    mutate(stimulus_cat = case_when(stimulus_cat == 1 ~ "inanimate (not food)",
                                 stimulus_cat == 2 ~ "inanimate (food)",
                                 stimulus_cat == 3 ~ "social"))
  
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
  
  pred_pid_time_summary <- pred_pid_time_long %>% 
    group_by(condition, time_c, resp) %>% 
    summarise(med = median(est), Species = unique(Species), Participant.name = unique(Participant.name), stimulus_cat = unique(stimulus_cat)) %>% 
    mutate(stimulus_cat = case_when(stimulus_cat == 1 ~ "inanimate (not food)",
                                 stimulus_cat == 2 ~ "inanimate (food)",
                                 stimulus_cat == 3 ~ "social"))
  
  ######
  
  ## Change factor level order ###########
  pred_time_long_summary$resp <- factor(pred_time_long_summary$resp, levels = (c("Agent", "Patient", "Other")))
  pred_pid_time_summary$resp <- factor(pred_pid_time_summary$resp, levels = (c("Agent", "Patient", "Other")))
  
  ## Export estimates #####################
  write_csv(pred_time_long_summary %>% ungroup() %>% dplyr::select(-c(condition)), file = "pred_time_series.csv")
  
  #########################################
  # Now we can plot! ######################
  
  p_adult <- ggplot(filter(pred_time_long_summary, Species == "Human"), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Human"), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("") +
    ylab("") +
    AOI_col_scale +
    AOI_fill_scale +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_infant <- ggplot(filter(pred_time_long_summary, Species == "Human baby"), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Human baby"), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("") +
    ylab("") +
    AOI_col_scale +
    AOI_fill_scale +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_chimp <- ggplot(filter(pred_time_long_summary, Species == "Chimpanzee"), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Chimpanzee"), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("") +
    ylab("") +
    AOI_col_scale +
    AOI_fill_scale +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_gorilla <- ggplot(filter(pred_time_long_summary, Species == "Gorilla"), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Gorilla"), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("") +
    ylab("") +
    AOI_col_scale +
    AOI_fill_scale +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_orang <- ggplot(filter(pred_time_long_summary, Species == "Orangutan"), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Orangutan"), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("Time since A -> P action") +
    ylab("") +
    AOI_col_scale +
    AOI_fill_scale +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  #########################################
  p <- p_adult / p_infant / p_chimp / p_gorilla / p_orang
  
  ggsave(paste("time_course", bin_size, "new.pdf", sep = "_"), plot = p, height = 14, width = 13)
  
  return(p)
}
