plot_species_AP_md <- function(fit, data){
  
  d <- data$d_time
  data_list <- data$stan_data_time
  
  time_train <- sort(unique(d$time_c)) # input training points to GP
  time_test <- 0.5 # time values to predict
  
  # AP movement diff to predict
  AP_md <- c(0,1)
  
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
    AP_md = AP_md
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
                                         AP_md = newdata_time$AP_md[j]
    )
  }
  
  # Pivot array long
  pred_time_long <- pred_time %>% 
    cubelyr::as.tbl_cube(met_name = "est") %>% 
    as_tibble %>% 
    left_join(newdata_time) %>% 
    left_join(d_species)
  
  # Summarize, averaging over time
  pred_md_summary <- pred_time_long %>% 
    group_by(resp, stimulus_cat, Species, AP_md) %>% 
    summarise(med = median(est),
              lower_50 = PI((est), prob = 0.5)[1],
              upper_50 = PI((est), prob = 0.5)[2],
              lower_90 = PI((est), prob = 0.9)[1],
              upper_90 = PI((est), prob = 0.9)[2]) %>% 
    mutate(stimulus_cat = case_when(stimulus_cat == 1 ~ "inanimate (not food)",
                                    stimulus_cat == 2 ~ "inanimate (food)",
                                    stimulus_cat == 3 ~ "social"))
  
  pred_md_summary$AP_movement <- ifelse(pred_md_summary$AP_md == 0, "A/P", "A")
  pred_md_summary$AP_movement <- factor(pred_md_summary$AP_movement, levels = c("A/P", "A"))
  
  pred_md_summary$resp <- factor(pred_md_summary$resp, levels = (c("Agent", "Patient", "Other")))
  
  # Export estimates to .csv
  write_csv(pred_md_summary, file = "pred_AOI_MD.csv")
  
  #########################################
  # Now we can plot! ######################
  p_adult <- ggplot(filter(pred_md_summary, Species == "Human"), aes(x = AP_movement, y = med, color = resp, fill = resp)) + 
    facet_grid(resp ~ stimulus_cat) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, x = AP_movement), alpha = 0.4, width = 0) +
    geom_errorbar(aes(ymin = lower_50, ymax = upper_50, x = AP_movement), alpha = 0.4, width = 0, lwd = 1.4) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(resp, stimulus_cat))) +
    ylab("") +
    xlab("") +
    labs(title = "Human (adult)") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing.y = unit(3, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(), strip.text.y = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_infant <- ggplot(filter(pred_md_summary, Species == "Human baby"), aes(x = AP_movement, y = med, color = resp, fill = resp)) + 
    facet_grid(resp ~ stimulus_cat) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, x = AP_movement), alpha = 0.4, width = 0) +
    geom_errorbar(aes(ymin = lower_50, ymax = upper_50, x = AP_movement), alpha = 0.4, width = 0, lwd = 1.4) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(resp, stimulus_cat))) +
    ylab("") +
    xlab("") +
    labs(title = "Human (infant)") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing.y = unit(3, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_chimp <- ggplot(filter(pred_md_summary, Species == "Chimpanzee"), aes(x = AP_movement, y = med, color = resp, fill = resp)) + 
    facet_grid(resp ~ stimulus_cat) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, x = AP_movement), alpha = 0.4, width = 0) +
    geom_errorbar(aes(ymin = lower_50, ymax = upper_50, x = AP_movement), alpha = 0.4, width = 0, lwd = 1.4) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(resp, stimulus_cat))) +
    xlab("") +
    ylab("Pr(AOI Fixation)") +
    labs(title = "Chimpanzee") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing.y = unit(3, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(),strip.text.x = element_blank(), strip.text.y = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_gorilla <- ggplot(filter(pred_md_summary, Species == "Gorilla"), aes(x = AP_movement, y = med, color = resp, fill = resp)) + 
    facet_grid(resp ~ stimulus_cat) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, x = AP_movement), alpha = 0.4, width = 0) +
    geom_errorbar(aes(ymin = lower_50, ymax = upper_50, x = AP_movement), alpha = 0.4, width = 0, lwd = 1.4) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(resp, stimulus_cat))) +
    ylab("") +
    xlab("") +
    labs(title = "Gorilla") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing.y = unit(3, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(), strip.text.x = element_blank(), axis.text.y = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  p_orang <- ggplot(filter(pred_md_summary, Species == "Orangutan"), aes(x = AP_movement, y = med, color = resp, fill = resp)) + 
    facet_grid(resp ~ stimulus_cat) +
    geom_errorbar(aes(ymin = lower_90, ymax = upper_90, x = AP_movement), alpha = 0.4, width = 0) +
    geom_errorbar(aes(ymin = lower_50, ymax = upper_50, x = AP_movement), alpha = 0.4, width = 0, lwd = 1.4) +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(resp, stimulus_cat))) +
    ylab("") +
    xlab("") +
    labs(title = "Orangutan") +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    AOI_col_scale +
    AOI_fill_scale +
    theme_minimal(base_size = 15) +
    theme(panel.spacing.y = unit(3, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", legend.title = element_blank(), strip.text.x = element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  #########################################
  p <- (p_adult + p_infant) / (p_chimp + p_gorilla) / (p_orang + plot_spacer())
  
  ggsave("AP_md_species.pdf", plot = p, height = 14, width = 13)
  
  return(p)
}
