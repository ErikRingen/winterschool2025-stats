plot_species_time_AOI <- function(preds, AOI) {
  
  pred_time_long <- preds$pred_time_long
  pred_pid_time_long <- preds$pred_pid_time_long
  
  #########################################
  ## Summary df for the prob of all three AOIs
  pred_time_long_summary <- pred_time_long %>% 
    group_by(condition, time_c, resp) %>% 
    summarise(med = median(est), lower_50 = PI((est), prob = 0.5)[1],
              upper_50 = PI((est), prob = 0.5)[2],
              lower_90 = PI((est), prob = 0.9)[1],
              upper_90 = PI((est), prob = 0.9)[2], Species = unique(Species), stimulus_cat = unique(stimulus_cat))
  
  ## Export estimates #####################
  write_csv(filter(pred_time_long_summary, resp == AOI) %>% ungroup() %>% dplyr::select(-c(condition)), file = paste0("pred_time_series_", AOI, ".csv"))
  
  # Summary df for all three AOIs, participant-specific
  pred_pid_time_summary <- pred_pid_time_long %>% 
    group_by(condition, time_c, resp) %>% 
    summarise(med = median(est), Species = unique(Species), Participant.name = unique(Participant.name), stimulus_cat = unique(stimulus_cat))
  
  p_adult <- ggplot(filter(pred_time_long_summary, Species == "Human" & resp == AOI), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = time_c), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = time_c), alpha = 0.4, color = NA) +
    geom_line(aes(x = time_c, y = med), lwd = 1.3, alpha = 0.75) + 
    geom_line(data = filter(pred_pid_time_summary, Species == "Human" & resp == AOI), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), color = "black",lwd = 0.5, alpha = 0.3, linetype = "solid") +
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
            element_text(hjust = 0.5)) 
  
  p_infant <- ggplot(filter(pred_time_long_summary, Species == "Human baby" & resp == AOI), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = time_c), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = time_c), alpha = 0.4, color = NA) +
    geom_line(aes(x = time_c, y = med), lwd = 1.3, alpha = 0.75) + 
    geom_line(data = filter(pred_pid_time_summary, Species == "Human baby" & resp == AOI), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), color = "black", lwd = 0.5, alpha = 0.3, linetype = "solid") +
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
  
  p_chimp <- ggplot(filter(pred_time_long_summary, Species == "Chimpanzee" & resp == AOI), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = time_c), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = time_c), alpha = 0.4, color = NA) +
    geom_line(aes(x = time_c, y = med), lwd = 1.3, alpha = 0.75) + 
    geom_line(data = filter(pred_pid_time_summary, Species == "Chimpanzee" & resp == AOI), aes(x = time_c, y = med, group = interaction(Participant.name, resp)),color = "black", lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("") +
    ylab(paste0("Pr(", AOI, " Fixation)")) +
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
  
  p_gorilla <- ggplot(filter(pred_time_long_summary, Species == "Gorilla" & resp == AOI), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = time_c), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = time_c), alpha = 0.4, color = NA) +
    geom_line(aes(x = time_c, y = med), lwd = 1.3, alpha = 0.75) + 
    geom_line(data = filter(pred_pid_time_summary, Species == "Gorilla" & resp == AOI), aes(x = time_c, y = med, group = interaction(Participant.name, resp)),color = "black", lwd = 0.5, alpha = 0.3, linetype = "solid") +
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
  
  p_orang <- ggplot(filter(pred_time_long_summary, Species == "Orangutan" & resp == AOI), aes(x = time_c, y = med, fill = resp, color = resp)) + 
    facet_grid(~ stimulus_cat, scales = "free_y", space = "free_y") +
    geom_ribbon(aes(ymin = lower_50, ymax = upper_50, x = time_c), alpha = 0.4, color = NA) +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90, x = time_c), alpha = 0.4, color = NA) +
    geom_line(aes(x = time_c, y = med), lwd = 1.3, alpha = 0.75) + 
    geom_line(data = filter(pred_pid_time_summary, Species == "Orangutan" & resp == AOI), aes(x = time_c, y = med, group = interaction(Participant.name, resp)), color = "black", lwd = 0.5, alpha = 0.3, linetype = "solid") +
    xlab("Proportion time since action starts") +
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
  
  #########################################
  p <- p_adult / p_infant / p_chimp / p_gorilla / p_orang
  
  ggsave(paste0("time_course_", AOI, ".pdf"), plot = p, height = 14, width = 13)
  
  return(p)
}
