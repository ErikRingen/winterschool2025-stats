plot_species_time_OR_APO <- function(preds) {
  
  pred_time_long <- preds$pred_time_long
  pred_pid_time_long <- preds$pred_pid_time_long
  
  ###########################################################
  ### Odds ratio (agent / patient) ##########################
  pred_time_long_OR_summary <- pred_time_long %>% 
    pivot_wider(names_from = resp, values_from = est) %>% 
    mutate(est = logit(Agent + Patient)) %>% 
    group_by(condition, time_c) %>% 
    summarise(med = median(est),
              lower_50 = PI((est), prob = 0.5)[1],
              upper_50 = PI((est), prob = 0.5)[2],
              lower_90 = PI((est), prob = 0.9)[1],
              upper_90 = PI((est), prob = 0.9)[2],
              Species = unique(Species), stimulus_cat = unique(stimulus_cat)) %>% 
    mutate(
      
      lower_50_1 = case_when(
        (lower_50 < 0 & upper_50 < 0) ~ lower_50,
        (lower_50 > 0 & upper_50 > 0) ~ -9999,
        (lower_50 < 0 & upper_50 > 0) ~ lower_50
      ),
      lower_50_2 = case_when(
        lower_50 < 0 & upper_50 < 0 ~ -9999,
        lower_50 > 0 & upper_50 > 0 ~ lower_50,
        lower_50 < 0 & upper_50 > 0 ~ 0
      ),
      upper_50_1 = case_when(
        lower_50 < 0 & upper_50 < 0 ~ upper_50,
        lower_50 > 0 & upper_50 > 0 ~ -9999,
        lower_50 < 0 & upper_50 > 0 ~ 0
      ),
      upper_50_2 = case_when(
        lower_50 < 0 & upper_50 < 0 ~ -9999,
        lower_50 > 0 & upper_50 > 0 ~ upper_50,
        lower_50 < 0 & upper_50 > 0 ~ upper_50
      ),
      lower_90_1 = case_when(
        (lower_90 < 0 & upper_90 < 0) ~ lower_90,
        (lower_90 > 0 & upper_90 > 0) ~ -9999,
        (lower_90 < 0 & upper_90 > 0) ~ lower_90
      ),
      lower_90_2 = case_when(
        lower_90 < 0 & upper_90 < 0 ~ -9999,
        lower_90 > 0 & upper_90 > 0 ~ lower_90,
        lower_90 < 0 & upper_90 > 0 ~ 0
      ),
      upper_90_1 = case_when(
        lower_90 < 0 & upper_90 < 0 ~ upper_90,
        lower_90 > 0 & upper_90 > 0 ~ -9999,
        lower_90 < 0 & upper_90 > 0 ~ 0
      ),
      upper_90_2 = case_when(
        lower_90 < 0 & upper_90 < 0 ~ -9999,
        lower_90 > 0 & upper_90 > 0 ~ upper_90,
        lower_90 < 0 & upper_90 > 0 ~ upper_90
      )
    ) %>% 
    mutate(across(where(is.numeric), ~na_if(., -9999)))
  
  ## Export estimates #####################
  write_csv(pred_time_long_OR_summary %>% ungroup() %>% dplyr::select(-c(condition, upper_50_1, upper_50_2, upper_90_1, upper_90_2, lower_50_1, lower_50_2, lower_90_1, lower_90_2)), file = "pred_time_series_OR_APO.csv")
  
  ## odds ratios for PIDs
  pred_pid_time_summary <- pred_pid_time_long %>% 
    pivot_wider(names_from = resp, values_from = est) %>% 
    mutate(est = logit(Agent + Patient)) %>% 
    group_by(pid, condition, time_c) %>% 
    summarise(med = median(est), Species = unique(Species), Participant.name = unique(Participant.name), stimulus_cat = unique(stimulus_cat))
  
  #########################################
  # Now we can plot! ######################
  
  p_adult <-  ggplot(filter(pred_time_long_OR_summary, Species == "Human"), aes(x = time_c, y = med)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower_50_1, ymax = upper_50_1, x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = lower_50_2, ymax = upper_50_2, x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_ribbon(aes(ymin = lower_90_1, ymax = upper_90_1, x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = lower_90_2, ymax = upper_90_2, x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Human"), aes(x = time_c, y = med, group = Participant.name), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    ylab("") +
    xlab("") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous() +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5)) 
  
  
  p_infant <-  ggplot(filter(pred_time_long_OR_summary, Species == "Human baby"), aes(x = time_c, y = med)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower_50_1, ymax = upper_50_1, x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = lower_50_2, ymax = upper_50_2, x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_ribbon(aes(ymin = lower_90_1, ymax = upper_90_1, x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = lower_90_2, ymax = upper_90_2, x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Human baby"), aes(x = time_c, y = med, group = Participant.name), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    ylab("") +
    xlab("") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous() +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_chimp <-  ggplot(filter(pred_time_long_OR_summary, Species == "Chimpanzee"), aes(x = time_c, y = med)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, NA, lower_50), ymax = ifelse(upper_50 < 0, upper_50, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, lower_50, 0), ymax = ifelse(upper_50 > 0, upper_50, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, NA, lower_90), ymax = ifelse(upper_90 < 0, upper_90, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, lower_90, 0), ymax = ifelse(upper_90 > 0, upper_90, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Chimpanzee"), aes(x = time_c, y = med, group = Participant.name), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    ylab("Log Odds of Agent/Patient vs. Other Fixation") +
    xlab("") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous() +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_gorilla <-  ggplot(filter(pred_time_long_OR_summary, Species == "Gorilla"), aes(x = time_c, y = med)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, NA, lower_50), ymax = ifelse(upper_50 < 0, upper_50, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, lower_50, 0), ymax = ifelse(upper_50 > 0, upper_50, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, NA, lower_90), ymax = ifelse(upper_90 < 0, upper_90, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, lower_90, 0), ymax = ifelse(upper_90 > 0, upper_90, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Gorilla"), aes(x = time_c, y = med, group = Participant.name), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    ylab("") +
    xlab("") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous() +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", axis.text.x = element_blank(), legend.title = element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  p_orang <- ggplot(filter(pred_time_long_OR_summary, Species == "Orangutan"), aes(x = time_c, y = med)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, NA, lower_50), ymax = ifelse(upper_50 < 0, upper_50, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_50 > 0, lower_50, 0), ymax = ifelse(upper_50 > 0, upper_50, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, NA, lower_90), ymax = ifelse(upper_90 < 0, upper_90, 0), x = time_c), alpha = 0.4, fill = "#107D81", color = NA) +
    geom_ribbon(aes(ymin = ifelse(lower_90 > 0, lower_90, 0), ymax = ifelse(upper_90 > 0, upper_90, NA), x = time_c), alpha = 0.4, fill = "#d25d28", color = NA) +
    geom_line(lwd = 1.3, alpha = 0.75) +
    geom_line(data = filter(pred_pid_time_summary, Species == "Orangutan"), aes(x = time_c, y = med, group = Participant.name), lwd = 0.5, alpha = 0.3, linetype = "solid") +
    ylab("") +
    xlab("Proportion time since action starts") +
    scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
    scale_y_continuous() +
    theme_minimal(base_size = 22) +
    theme(panel.spacing = unit(1.5, "lines"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", legend.title = element_blank(),
          plot.title = 
            element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text.x = element_blank()) 
  
  #########################################
  p <- p_adult / p_infant / p_chimp / p_gorilla / p_orang
  
  ggsave("time_course_OR_APO.pdf", plot = p, height = 14, width = 13)
  
  return(p)
}
