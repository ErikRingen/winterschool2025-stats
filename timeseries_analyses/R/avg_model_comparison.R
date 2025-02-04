avg_model_comparison <- function(d) {
  
  d$stimulus_category <- case_when(
    d$stimulus_cat == 1 ~ "inanimate, not food",
    d$stimulus_cat == 2 ~ "inanimate, food",
    d$stimulus_cat == 3 ~ "social"
  )
  
  d$size <- with(d, Agent_Sum + Patient_Sum + Other_Sum)
  d$y <- with(d, cbind(Agent_Sum, Patient_Sum, Other_Sum))
  
  d$y_prop <- d$y / d$size
  
  adjust_proportions <- function(proportions, epsilon=1e-4) {
    
    # Adjust the proportions based on the conditions
    adjusted <- ifelse(proportions == 0, proportions + epsilon, 
                       ifelse(proportions == 1, proportions - epsilon, proportions))
    
    # Normalize so they sum to 1
    adjusted <- adjusted / sum(adjusted)
    
    return(adjusted)
  }
  
  V <- Vectorize(adjust_proportions)
  d$y_prop2 <- t(apply(d$y_prop, 1, adjust_proportions))
  
  ###############################
  # create indices
  species_id <- match(d$Species, unique(d$Species))
  pid <- match(d$Participant.name, unique(d$Participant.name))
  stimulus_id <- match(d$Stimulus, unique(d$Stimulus))
  stimulus_cat <- match(d$stimulus_cat, unique(d$stimulus_cat))
  fs_id <- match(d$FootageSpecies, unique(d$FootageSpecies))
  fs_species_id <- match(paste(d$FootageSpecies, d$Species), unique(paste(d$FootageSpecies, d$Species)))
  
  N_species <- max(species_id)
  N_pid <- max(pid)
  N_stimuli <- max(stimulus_id)
  N_fs <- max(fs_id)
  
  data_list <- list(
    y = d$y_prop2,
    K = 3,
    N = nrow(d),
    N_species = N_species,
    N_stimuli = N_stimuli,
    N_pid = N_pid,
    N_fs = N_fs,
    species_id = species_id,
    pid = pid,
    stimulus_id = stimulus_id,
    fs_id = fs_id,
    fs_species_id = fs_species_id,
    stimulus_cat = stimulus_cat,
    log_AP_size = d$log_AP_size,
    AP_md = d$AP_md
  )
  
  mod_null <- stan_model(file = "stan/m_dirichlet_pred_null.stan")
  mod <- stan_model(file = "stan/m_dirichlet_pred.stan") 
  
  n_iter <- 2000
  fit_null <- sampling(mod_null, data = data_list, chains = 8, cores = 8, init = "0", iter = n_iter)
  fit <- sampling(mod, data = data_list, chains = 8, cores = 8, init = "0", iter = n_iter)
  
  # extract log lik arrays
  log_lik <- extract.samples(fit, pars = "log_lik")
  log_lik_null <- extract.samples(fit_null, pars = "log_lik")
  
  # subset by species group
  
  ### Apes #######
  log_lik_ape <- log_lik$log_lik[, data_list$species_id %in% c(1, 2, 5)]
  log_lik_ape_null <- log_lik_null$log_lik[, data_list$species_id %in% c(1, 2, 5)]
  
  loo_ape_null <- loo(log_lik_ape_null, r_eff = loo::relative_eff(exp(log_lik_ape_null), chain_id = rep(1:8, each = n_iter/2)))
  loo_ape <- loo(log_lik_ape, r_eff = loo::relative_eff(exp(log_lik_ape), chain_id = rep(1:8, each = n_iter/2)))
  
  weights_ape <- loo::loo_model_weights(list(loo_ape_null, loo_ape), method = "stacking")
  ################
  
  ### Human adult #######
  log_lik_adult <- log_lik$log_lik[, data_list$species_id == 4]
  log_lik_adult_null <- log_lik_null$log_lik[, data_list$species_id == 4]
  
  loo_adult_null <- loo(log_lik_adult_null, r_eff = loo::relative_eff(exp(log_lik_adult_null), chain_id = rep(1:8, each = n_iter/2)))
  loo_adult <- loo(log_lik_adult, r_eff = loo::relative_eff(exp(log_lik_adult), chain_id = rep(1:8, each = n_iter/2)))
  
  weights_adult <- loo::loo_model_weights(list(loo_adult_null, loo_adult), method = "stacking")
  ################
  
  ### Human infant #######
  log_lik_infant <- log_lik$log_lik[, data_list$species_id == 3]
  log_lik_infant_null <- log_lik_null$log_lik[, data_list$species_id == 3]
  
  loo_infant_null <- loo(log_lik_infant_null, r_eff = loo::relative_eff(exp(log_lik_infant_null), chain_id = rep(1:8, each = n_iter/2)))
  loo_infant <- loo(log_lik_infant, r_eff = loo::relative_eff(exp(log_lik_infant), chain_id = rep(1:8, each = n_iter/2)))
  
  weights_infant <- loo::loo_model_weights(list(loo_infant_null, loo_infant), method = "stacking")
  ################
  
  return(list(
    loo_ape = loo_ape,
    loo_ape_null = loo_ape_null,
    weights_ape = weights_ape,
    
    loo_adult = loo_adult,
    loo_adult_null = loo_adult_null,
    weights_adult = weights_adult,
    
    loo_infant = loo_infant,
    loo_infant_null = loo_infant_null,
    weights_infant = weights_infant
  ))
}



