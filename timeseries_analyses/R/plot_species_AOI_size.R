plot_species_AOI_size <- function(fit, data){
  
  d <- data$d_time
  data_list <- data$stan_data_time
  
  time_train <- sort(unique(d$time_c)) # input training points to GP
  time_test <- 0.5 # time values to predict
  
  # Size of AP movement diff to predict
  AOI_seq <- seq(from = min(d$log_AP_size), to = max(d$log_AP_size), length.out = 20)
  
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
    stimulus_cat = 1:3,
    log_AP_size = AOI_seq
  )
  
  newdata_time <- newdata_time %>% 
    mutate(
      condition = 1:n()
    )
  
  n_conditions <- max(newdata_time$condition)
  
  pred_time <- array(NA, dim = c(n_samps, 3, length(time_test), n_conditions), dimnames = list(samps = 1:n_samps, resp = c("Agent", "Patient", "Other"), time_c = time_test, condition = 1:n_conditions))
  
  for (j in 1:n_conditions) {
    pred_time[, , , j] = timeseries_pred(post = post,
                                         time_train = time_train,
                                         time_test = time_test,
                                         n_samps = n_samps,
                                         species_id = newdata_time$species_id[j],
                                         species_RE = T,
                                         stimulus_cat = newdata_time$stimulus_cat[j],
                                         log_AP_size = newdata_time$log_AP_size[j]
    )
  }
  
  # Pivot array long
  pred_time_long <- pred_time %>% 
    cubelyr::as.tbl_cube(met_name = "est") %>% 
    as_tibble %>% 
    left_join(newdata_time) %>% 
    left_join(d_species)
  
  # Summarize, averaging over time
  pred_AOI_summary <- pred_time_long %>% 
    group_by(resp, stimulus_cat, Species, log_AP_size) %>% 
    summarise(med = median(est),
              lower_50 = PI((est), prob = 0.5)[1],
              upper_50 = PI((est), prob = 0.5)[2],
              lower_90 = PI((est), prob = 0.9)[1],
              upper_90 = PI((est), prob = 0.9)[2]) 
  
  d_stim_cat <- d %>% 
    group_by(stimulus_cat) %>% 
    summarise(min_cat = min(log_AP_size), max_cat = max(log_AP_size))
  
  # Constrain the range of the predictor depending on stimulus category
  pred_AOI_summary <- pred_AOI_summary %>% 
    left_join(d_stim_cat) %>% 
    filter(log_AP_size >= min_cat & log_AP_size <= max_cat) %>% 
    mutate(stimulus_cat = case_when(stimulus_cat == 1 ~ "inanimate (not food)",
                                    stimulus_cat == 2 ~ "inanimate (food)",
                                    stimulus_cat == 3 ~ "social"))
  
  # Export estimates to .csv
  write_csv(pred_AOI_summary, file = "pred_AOI_size.csv")
  
  #########################################
  # Now we can plot! ######################
  p_adult <- ggplot(filter(pred_AOI_summary, Species == "Human"), aes(x = log_AP_size, y = med, color = resp, fill = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_x") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    ylab("") +
    xlab("") +
    labs(title = "Human (adult)") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_infant <- ggplot(filter(pred_AOI_summary, Species == "Human baby"), aes(x = log_AP_size, y = med, color = resp, fill = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_x") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    ylab("") +
    xlab("") +
    labs(title = "Human (infant)") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(), strip.text.x = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_chimp <- ggplot(filter(pred_AOI_summary, Species == "Chimpanzee"), aes(x = log_AP_size, y = med, color = resp, fill = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_x") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    xlab("") +
    ylab("Pr(AOI Fixation)") +
    labs(title = "Chimpanzee") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(), strip.text.x = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_gorilla <- ggplot(filter(pred_AOI_summary, Species == "Gorilla"), aes(x = log_AP_size, y = med, color = resp, fill = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_x") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    ylab("") +
    xlab("") +
    labs(title = "Gorilla") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(), strip.text.x = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_orang <- ggplot(filter(pred_AOI_summary, Species == "Orangutan"), aes(x = log_AP_size, y = med, color = resp, fill = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_x") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = log_AP_size), alpha = 0.15, color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    ylab("") +
    xlab("log(Agent Size / Patient Size)") +
    labs(title = "Orangutan") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(), strip.text.x = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  #########################################
  p <- p_adult / p_infant / p_chimp / p_gorilla / p_orang
  
  ggsave("AOI_size_species.pdf", plot = p, height = 14, width = 13)
  
  return(p)
}
