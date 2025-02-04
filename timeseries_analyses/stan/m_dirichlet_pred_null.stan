functions {
  /* dirichlet-logit log-PDF
  * Args:
    *   y: vector of real response values
  *   mu: vector of category logit probabilities
  *   phi: precision parameter
  * Returns:
    *   a scalar to be added to the log posterior
  */
    real dirichlet_logit_lpdf(vector y, vector mu, real phi) {
      return dirichlet_lpdf(y | softmax(mu) * phi);
    }
}

data{
  int K; // num AOI
  int N; // num obs
  int N_species;
  int N_pid;
  int N_stimuli;
  int N_fs;
  array[N] vector[K] y; // proportion AOI fixation
  vector[N] log_AP_size;
  vector[N] AP_md;
  array[N] int stimulus_cat;
  array[N] int species_id;
  array[N] int pid;
  array[N] int stimulus_id;
  array[N] int fs_id; // footage species
  array[N] int fs_species_id; // footage_species:species
}

parameters{
  real<lower=0> phi; // precision parameter
  matrix[5, K - 2] b; // fixed effects
  
  // random effects
  matrix[7 * (K - 2), N_species] species_z; 
  matrix[6 * (K - 2), N_pid] pid_z;
  
  matrix[(K - 2), N_stimuli] stimulus_z;
  matrix[(K - 2), N_fs] fs_z;
  matrix[(K - 2), N_species * N_fs] fs_species_z;
  
  // variance components
  vector<lower=0>[7 * (K - 2)] sigma_species;
  vector<lower=0>[6 * (K - 2)] sigma_pid;
  vector<lower=0>[K-2] sigma_stimulus;
  vector<lower=0>[K-2] sigma_fs;
  vector<lower=0>[K-2] sigma_fs_species;
}

transformed parameters{
  array[N] vector[K] mu;
  
  // random effects, scaled
  matrix[N_species, 7 * (K - 2)] species_v = (diag_matrix(sigma_species) * species_z)'; 
  matrix[N_pid, 6  * (K - 2)] pid_v = (diag_matrix(sigma_pid) * pid_z)';
  
  matrix[N_stimuli, (K - 2)] stimulus_v = (diag_matrix(sigma_stimulus) * stimulus_z)'; 
  matrix[N_fs, (K - 2)] fs_v = (diag_matrix(sigma_fs) * fs_z)'; 
  matrix[N_species*N_fs, (K - 2)] fs_species_v = (diag_matrix(sigma_fs_species) * fs_species_z)'; 
  
  for (i in 1:N) {
      mu[i, 1] = b[stimulus_cat[i], 1] + species_v[species_id[i], 1] + species_v[species_id[i], 1 + stimulus_cat[i]] + pid_v[pid[i], 1] + stimulus_v[stimulus_id[i], 1] + fs_v[fs_id[i], 1] + fs_species_v[fs_species_id[i], 1] + (b[4, 1] + species_v[species_id[i], 5] + pid_v[pid[i], 5])*log_AP_size[i] + (b[5, 1] + species_v[species_id[i], 6] + pid_v[pid[i], 6])*AP_md[i];

      mu[i, 2] = mu[i, 1]; // null model where agent and patient have the same linear model
      
      mu[i, 3] = 0; // ref category 
    }
}

model{
  phi ~ gamma(0.1, 0.1);
  to_vector(b) ~ std_normal();
  
  to_vector(species_z) ~ std_normal();
  to_vector(pid_z) ~ std_normal();
  to_vector(fs_z) ~ std_normal();
  to_vector(stimulus_z) ~ std_normal();
  to_vector(fs_species_z) ~ std_normal();
  
  sigma_species ~ exponential(1);
  sigma_pid ~ exponential(1);
  sigma_fs ~ exponential(1);
  sigma_stimulus ~ exponential(1);
  sigma_fs_species ~ exponential(1);
  
  for (i in 1:N) target += dirichlet_logit_lpdf(y[i] | mu[i], phi * exp(species_v[species_id[i], 7]));
}

generated quantities{
  vector[N] log_lik;
  
  for (i in 1:N) log_lik[i] = dirichlet_logit_lpdf(y[i] | mu[i], phi * exp(species_v[species_id[i], 7]));
}
