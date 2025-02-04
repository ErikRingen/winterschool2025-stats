data{
  int K; // number of AOI categories (Agent, Patient, Other)
  int N_obs; // number of observations
  int N_species; // number of species
  int N_fs; // number of footage species
  int N_pid; // number of participants
  int N_stimuli; // number of stimulus videos
  int N_trials; // number of trials
  int N_times; // number of unique time points
  matrix[N_times, N_times] time_dist_mat; // euclidian distance between time points
  array[N_obs] int time_id; // ID of time
  array[N_obs] real time; // value of time
  array[N_obs] int species_id; // ID of participant species
  array[N_obs] int fs_id; // ID of footage species
  array[N_obs] int pid; // ID of participant
  array[N_obs] int stimulus_id; // ID of stimulus
  array[N_obs] int trial_id; // ID of trial
  array[N_obs, 3] real y_prop; // proportion of Agent/Patient/Other fixations
  array[N_obs] int stimulus_cat; // indicator of inanimate (not food), inanimate (food), or social
  array[N_obs] real log_AP_size; // ln(Agent Size / Patient Size)
  array[N_obs] int AP_md; // indicator of agent moves more than patient (1) or not (0)
}

parameters{
  // overall effects, col 1 belongs to agent, col 2 belongs to patient
  // rows correspond to stimulus category 
  // intercepts
  matrix[3, K - 1] b0; 
  // fixed effects
  matrix[3, K - 1] b_AP_size;
  matrix[3, K - 1] b_AP_movement;
  // Varying effects, time-invariant
  matrix<lower=0>[3, (K - 1)] sigmasq_b0; // variance components
  matrix<lower=0>[3, (K - 1)] sigmasq_AOI_size; // variance components
  matrix<lower=0>[3, (K - 1)] sigmasq_AP_md; // variance components
  // Dirichlet decomposition parameters
  array[3, (K - 1)] simplex[6] phi_b0;
  array[3, (K - 1)] simplex[2] phi_AOI_size;
  array[3, (K - 1)] simplex[2] phi_AP_md;
  
  ////////////////////////////////////////
  // Species varying effects /////////////
  array[N_species, 3, (K - 1)] real species_b0_z; // species-level varying effects, unscaled and uncorrelated
  array[N_species, 3, (K - 1)] real species_AP_size_z; // species-level varying effects, unscaled and uncorrelated
  array[N_species, 3, (K - 1)] real species_AP_movement_z; // species-level varying effects, unscaled and uncorrelated

  // Gaussian processes for species
  array[N_species, 3] matrix[N_times, (K - 1)] time_species_z;
  matrix<lower=0>[3, K - 1] alpha_time_species;
  ////////////////////////
  
  // PID varying effects /////////////////
  array[N_pid, 3, (K - 1)] real pid_b0_z; // species-level varying effects, unscaled and uncorrelated
  array[N_pid, 3, (K - 1)] real pid_AP_size_z; // species-level varying effects, unscaled and uncorrelated
  array[N_pid, 3, (K - 1)] real pid_AP_movement_z; // species-level varying effects, unscaled and uncorrelated
  
  // Gaussian processes for PIDs
  array[N_pid, 3] matrix[N_times, (K - 1)] time_pid_z;
  matrix<lower=0>[3, K - 1] alpha_time_pid;
  ////////////////////////
  
  // Stimulus varying effects ////////////
  matrix[N_stimuli, (K - 1)] stimulus_b0_z; // stimulus-level varying effects, unscaled and uncorrelated

  // Stimulus GPs
  array[N_stimuli] matrix[N_times, (K - 1)] time_stimulus_z;
  vector<lower=0>[K - 1] alpha_time_stimulus;
  ////////////////////////
  
  // Footage Species varying effects
  matrix[N_fs, (K - 1)] fs_b0_z;

  // Footage Species x Species varying effects
  array[N_fs, N_species, (K - 1)] real fs_species_b0_z;

  // Trial varying effects
  matrix[N_trials, (K - 1)] trial_b0_z; // trial-level varying effects, unscaled and uncorrelated

  // Trial GP
  array[N_trials] matrix[N_times, 2] time_trial_z;
  vector<lower=0>[(K - 1)] alpha_time_trial;
  
  // Gaussian process parameters /////
  // Smooth functions for time //
  array[3] matrix[N_times, (K - 1)] time_z;
  matrix<lower=0>[3, K - 1] alpha_time;
  
  matrix<lower=0>[3, (K - 1)] sigmasq_time_mu;
  matrix<lower=0>[3, (K - 1)] sigmasq_time;
  array[3, (K - 1)] simplex[4] phi_time;
}

transformed parameters{
  // In this block, we scale correlate the GP time varying effects
  // Gaussian process parameters (smooth time functions)
  array[3] matrix[N_times, (K - 1)] time_v;
  array[N_species, 3] matrix[N_times, (K - 1)] time_species_v;
  array[N_pid, 3] matrix[N_times, (K - 1)] time_pid_v;
  array[N_stimuli] matrix[N_times, (K - 1)] time_stimulus_v;
  array[N_trials] matrix[N_times, 2] time_trial_v;
  
  for (cat in 1:3) {
    for (k in 1:(K - 1)) {
    matrix[N_times, N_times] cov_mat;
    matrix[N_times, N_times] cov_mat_species;
    matrix[N_times, N_times] cov_mat_pid;
    matrix[N_times, N_times] cov_mat_stimulus;
    matrix[N_times, N_times] cov_mat_trial;
    
    for (i in 1:(N_times - 1)) {
      for (j in (i + 1):N_times) {
        cov_mat[i, j] = exp(-alpha_time[cat, k] * (time_dist_mat[i, j]));
        cov_mat[j, i] = cov_mat[i, j]; // symmetrical
        
        cov_mat_species[i, j] = exp(-alpha_time_species[cat, k] * (time_dist_mat[i, j]));
        cov_mat_species[j, i] = cov_mat_species[i, j]; // symmetrical
        
        cov_mat_pid[i, j] = exp(-alpha_time_pid[cat, k] * (time_dist_mat[i, j]));
        cov_mat_pid[j, i] = cov_mat_pid[i, j]; // symmetrical
        
        // stimulus and trial are always just one category
        cov_mat_stimulus[i, j] = exp(-alpha_time_stimulus[k] * (time_dist_mat[i, j]));
        cov_mat_stimulus[j, i] = cov_mat_stimulus[i, j]; // symmetrical
        
        cov_mat_trial[i, j] = exp(-alpha_time_trial[k] * (time_dist_mat[i, j]));
        cov_mat_trial[j, i] = cov_mat_trial[i, j]; // symmetrical
      }
    }
    
      for (d in 1:N_times) {
        cov_mat[d, d] = 1.000001;
        cov_mat_species[d, d] = 1.000001;
        cov_mat_pid[d, d] = 1.000001;
        cov_mat_stimulus[d, d] = 1.000001;
        cov_mat_trial[d, d] = 1.000001;
      }      
      // Now that we've composed the covariance matrices, use them to scale and correlate the time varying effects
      time_v[cat][ , k] = cholesky_decompose(cov_mat) * time_z[cat][ , k];
      
      for (i in 1:N_species) time_species_v[i, cat][, k] = cholesky_decompose(cov_mat_species) * time_species_z[i, cat][, k];
      for (i in 1:N_pid) time_pid_v[i, cat][, k] = cholesky_decompose(cov_mat_pid) * time_pid_z[i, cat][, k];
      for (i in 1:N_stimuli) time_stimulus_v[i][, k] = cholesky_decompose(cov_mat_stimulus) * time_stimulus_z[i][, k];
      for (i in 1:N_trials) time_trial_v[i][, k] = cholesky_decompose(cov_mat_trial) * time_trial_z[i][, k];
  }
  }
  //////////////////////////////////
}

model{
  // Priors ///////////////////////
  // Fixed effects ////
  to_vector(b0) ~ std_normal();
  to_vector(b_AP_size) ~ std_normal();
  to_vector(b_AP_movement) ~ std_normal();
  
  // Varying effects
  for (j in 1:3) {
    for (k in 1:(K - 1)) {
    to_vector(species_b0_z[,j,k]) ~ std_normal();
    to_vector(species_AP_size_z[,j,k]) ~ std_normal();
    to_vector(species_AP_movement_z[,j,k]) ~ std_normal();
    to_vector(pid_b0_z[,j,k]) ~ std_normal();
    to_vector(pid_AP_size_z[,j,k]) ~ std_normal();
    to_vector(pid_AP_movement_z[,j,k]) ~ std_normal();
  }
  }
  
  to_vector(stimulus_b0_z) ~ std_normal();
  to_vector(fs_b0_z) ~ std_normal();
  to_vector(trial_b0_z) ~ std_normal();
  
  for (k in 1:(K-1)) {
  for (j in 1:N_species) {
    to_vector(fs_species_b0_z[,j,k]) ~ std_normal();
  }
  }
  
  for (j in 1:3) to_vector(time_z[j]) ~ std_normal();
    for (i in 1:N_trials) {
      to_vector(time_trial_z[i]) ~ std_normal();
    }
    
   for (i in 1:N_species) {
    for (j in 1:3) {
      to_vector(time_species_z[i, j]) ~ std_normal();   
    }
  }
  
  for (i in 1:N_pid) {
    for (j in 1:3) {
      to_vector(time_pid_z[i, j]) ~ std_normal();
    }
  }
  
  for (i in 1:N_stimuli) to_vector(time_stimulus_z[i]) ~ std_normal();
  
  // Variance components
  to_vector(sigmasq_b0) ~ exponential(1);
  to_vector(sigmasq_AOI_size) ~ exponential(1);
  to_vector(sigmasq_AP_md) ~ exponential(1);
  to_vector(sigmasq_time_mu) ~ exponential(1);
  to_vector(sigmasq_time) ~ exponential(1);
  
   // simplex priors
  for (j in 1:3)
    for (k in 1:(K - 3)) {
    to_vector(phi_time[j, k]) ~ dirichlet(rep_vector(2, 4)); 
    to_vector(phi_b0[j, k]) ~ dirichlet(rep_vector(2, 6));
    to_vector(phi_AOI_size[j, k]) ~ dirichlet(rep_vector(2, 2));
    to_vector(phi_AP_md[j, k]) ~ dirichlet(rep_vector(2, 2));
  }
  
  // Gaussian process parameters
  // Diffuse priors here allow functions to be either wiggly or rigid
  alpha_time_trial ~ normal(0, 10);
  to_vector(alpha_time_pid) ~ normal(0, 10);
  to_vector(alpha_time) ~ normal(0, 10);
  to_vector(alpha_time_species) ~ normal(0, 10);
  to_vector(alpha_time_stimulus) ~ normal(0, 10);

  //////////////////////////////
  // Likelihood loop //////////
  for (i in 1:N_obs) {
    vector[K] theta; // latent-scale probability of AOI fixation
    vector[K] lp; // log probability of each AOI fixation
    
    for (k in 1:(K - 1)) {
      real intercept_stem;
      real time_stem;
      
      // Time invariant parameters
      // Intercepts (b0)
      intercept_stem = b0[stimulus_cat[i], k] + 
      species_b0_z[species_id[i], stimulus_cat[i], k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][1]) + 
      pid_b0_z[pid[i], stimulus_cat[i], k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][2]) + 
      stimulus_b0_z[stimulus_id[i],k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][3]) + 
      trial_b0_z[trial_id[i], k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][4]) + 
      fs_b0_z[fs_id[i], k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][5]) + 
      fs_species_b0_z[fs_id[i], species_id[i], k] * sqrt(sigmasq_b0[stimulus_cat[i], k] * phi_b0[stimulus_cat[i], k][6]) +
      
      // Covariates (relative AOI size, AP movement diff)
      (b_AP_size[stimulus_cat[i], k] + 
      species_AP_size_z[species_id[i], stimulus_cat[i], k] * sqrt(sigmasq_AOI_size[stimulus_cat[i], k] * phi_AOI_size[stimulus_cat[i], k][1]) + 
      pid_AP_size_z[pid[i], stimulus_cat[i], k] * sqrt(sigmasq_AOI_size[stimulus_cat[i], k] * phi_AOI_size[stimulus_cat[i], k][2])) * log_AP_size[i] +
      
      (b_AP_movement[stimulus_cat[i], k] + 
      species_AP_movement_z[species_id[i], stimulus_cat[i], k] * sqrt(sigmasq_AP_md[stimulus_cat[i], k] * phi_AP_md[stimulus_cat[i], k][1]) + 
      pid_AP_movement_z[pid[i], stimulus_cat[i], k] * sqrt(sigmasq_AP_md[stimulus_cat[i], k] * phi_AP_md[stimulus_cat[i], k][2])) * AP_md[i];
      
      // Time-varying parameters: overall smooth + species smooth + pid smooth + stimulus smooth + trial smooth
      time_stem = time_v[stimulus_cat[i]][time_id[i], k] * sqrt(sigmasq_time_mu[stimulus_cat[i], k]) + 
      time_species_v[species_id[i], stimulus_cat[i]][time_id[i], k] * sqrt(sigmasq_time[stimulus_cat[i], k] * phi_time[stimulus_cat[i], k][1]) + 
      time_pid_v[pid[i], stimulus_cat[i]][time_id[i], k] * sqrt(sigmasq_time[stimulus_cat[i], k] * phi_time[stimulus_cat[i], k][2]) + 
      time_stimulus_v[stimulus_id[i]][time_id[i], k] * sqrt(sigmasq_time[stimulus_cat[i], k] * phi_time[stimulus_cat[i], k][3]) + 
      time_trial_v[trial_id[i]][time_id[i], k] * sqrt(sigmasq_time[stimulus_cat[i], k] * phi_time[stimulus_cat[i], k][4]);
      
      // Put stems together
      theta[k] = intercept_stem + time_stem;
      }
      theta[K] = 0; // ref category (other)

    // mix over latent fixation categories (likelihood * prior, log scale)
    for (k in 1:K) lp[k] = categorical_logit_lpmf(k | theta) + log(y_prop[i,k]);
    target += log_sum_exp(lp);
   }
}
