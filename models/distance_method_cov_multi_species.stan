functions {
  real partial_sum_lpmf(array[, ] int slice_abund_per_band, // modifications to avoid Stan warnings at compile
                        int start, int end,
                        int max_intervals,
                        array[] int bands_per_sample,
                        array[, ] real max_dist,
                        array[ , , , ] real log_tau,
                        array[] int method,
                        array[] int roadside,
                        array[] int forest,
                        array[] int species)
  {
    real lp = 0;
    int Pi_size = end - start + 1;
    int Pi_index = 1;
    matrix[Pi_size, max_intervals] Pi = rep_matrix(0, Pi_size, max_intervals);
    
    for (i in start:end)
    {
      for (k in 1:(bands_per_sample[i]-1)) // 
        {
          if(k > 1){
            Pi[Pi_index,k] = ((1 - exp(-(max_dist[i,k]^2 / exp(log_tau[method[i],forest[i],roadside[i],species[i]])^2))) - 
                                (1 - exp(-(max_dist[i,k - 1]^2 / exp(log_tau[method[i],forest[i],roadside[i],species[i]])^2)))) / 
              (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[method[i],forest[i],roadside[i],species[i]])^2)));
          }else{
            Pi[Pi_index,k] = (1 - exp(-(max_dist[i,k]^2 / exp(log_tau[method[i],forest[i],roadside[i],species[i]])^2))) /
              (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[method[i],forest[i],roadside[i],species[i]])^2)));
          }
        } 
      Pi[Pi_index,bands_per_sample[i]] = 1 - sum(Pi[Pi_index,]); // what if the final band was used as the constraint?
        
      lp = lp + multinomial_lupmf(slice_abund_per_band[Pi_index, ] | to_vector(Pi[Pi_index, ]));
      Pi_index = Pi_index + 1;
      
    }
    
    return lp;
  }
  
}


data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 1> n_methods;           // total number of methods
  int<lower = 1> n_species;           // total number of specise  
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> grainsize;           // grainsize for reduce_sum() function
  int<lower = 1> n_forests;           // number of unique values of proportion forest
  array[n_samples] int method;             // method being considered for each sample
  array[n_samples] int<lower = 1, upper = n_forests> forest; // open == 1 forested == n_forest, with n_forest-1 steps between
  array[n_samples] int<lower = 1, upper = 2> roadside;  // roadside == 2 offroad == 1 each sample
  
  array[n_samples] int<lower = 1, upper = n_species> species;  // roadside == 2 offroad == 1 each sample
  array[n_samples, max_intervals] int abund_per_band;// abundance in distance band k for sample i
  array[n_samples] int bands_per_sample; // number of distance bands for sample i
  array[n_samples, max_intervals] real max_dist; // max distance for distance band k

  //array[n_species] int mig_strat;        //migration strategy for each species
  //array[n_species] int habitat;        //habitat preference for each species
  array[n_species] real mass;  //log mass of species
  array[n_species] real pitch; //song pitch of species
}

transformed data{
  // transforming the survey-specific covariates
  array[n_samples] int roadside_z;
  array[n_samples] real forest_z;
  for(j in 1:n_samples){
  roadside_z[j] = roadside[j]-1; //transform integer values of roadside into offroad == 0 and roadside == 1
  forest_z[j] = (forest[j]-1)/10.0; //transform integer values of forest (1 : n_forests(11)) into 0 - 1 proportion forest
  }
  //print(roadside_z);
}


parameters {
  vector[n_methods] log_tau_method_raw;
  vector[n_species] log_tau_species_raw;
  real log_TAU;
    real BETA_forest;
      real BETA_roadside;
        real BETA_interaction;
  vector[n_species] beta_forest_raw;
  vector[n_species] beta_roadside_raw;
  vector[n_species] beta_interaction_raw;
  real<lower=0> sd_forest;
    real<lower=0> sd_roadside;
    real<lower=0> sd_interaction;
  real<lower=0> sd_log_tau_method;
  real<lower=0> sd_log_tau_species;
  
  //real beta_mig_strat;
  //real beta_habitat;
  real beta_mass;
  real beta_pitch;
}

transformed parameters{
  array[n_methods,n_forests,2,n_species] real log_tau;
  vector[n_species] beta_roadside;
  vector[n_species] beta_forest;
  vector[n_species] beta_interaction;
  
  beta_roadside = sd_roadside*beta_roadside_raw;
  beta_forest = sd_forest*beta_forest_raw;
  beta_interaction = sd_interaction*beta_interaction_raw;
  
  
  for(j in 1:n_samples){
  log_tau[method[j],forest[j],roadside[j],species[j]] = sd_log_tau_method*log_tau_method_raw[method[j]] + 
                sd_log_tau_species*log_tau_species_raw[species[j]] + 
                beta_forest[species[j]] * forest_z[j] +
                beta_roadside[species[j]] * roadside_z[j] +
                beta_interaction[species[j]] * roadside_z[j] * forest_z[j] +
                //beta_mig_strat * mig_strat[species[j]] +
                //beta_habitat * habitat[species[j]] + 
                beta_mass * mass[species[j]] +
                beta_pitch * pitch[species[j]] +
                log_TAU;
  }
}

model {
  log_tau_method_raw ~ std_normal();
  //sum(log_tau_method_raw) ~ normal(0,0.001*n_methods); // constraint may not be necessary
  log_tau_species_raw ~ std_normal();
  sum(log_tau_species_raw) ~ normal(0,0.001*n_species); // constraint may not be necessary

  log_TAU ~ normal(0,0.5); // weakly informative prior implying that mean EDRs are likely
  // between 45 and 230 m
  sd_log_tau_method ~ normal(0,0.25); // weakly informative prior implying that th average
  // among method variation is likely between -33% and +50% of the mean
  sd_log_tau_species ~ normal(0,0.1); // weakly informative prior implying that th average

  beta_forest_raw ~ std_normal();
  beta_roadside_raw ~ std_normal();
  beta_interaction_raw ~ std_normal();
  
  BETA_forest ~ normal(0,0.25);
  BETA_roadside ~ normal(0,0.25);
  BETA_interaction ~ normal(0,0.25);
  
  sd_forest ~ normal(0,0.1);
  sd_roadside ~ normal(0,0.1);
  sd_interaction ~ normal(0,0.1);
  

  //beta_mig_strat ~ normal(0,0.25);
  //beta_habitat ~ normal(0,0.25);
  beta_mass ~ normal(0,0.25);
  beta_pitch ~ normal(0,0.25);
  
  
  //print(log_tau);
  
  target += reduce_sum(partial_sum_lupmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_dist,
                       log_tau,
                       method,
                       roadside,
                       forest,
                       species);
}

generated quantities {
  vector[n_species] log_tau_open_offroad;
  vector[n_species] log_tau_forest_offroad;
  vector[n_species] log_tau_open_onroad;
  vector[n_species] log_tau_forest_onroad;
  
  for(s in 1:n_species){
  log_tau_open_offroad[s] = sd_log_tau_species*log_tau_species_raw[species[s]] + 
                //beta_mig_strat * mig_strat[species[s]] +
                //beta_habitat * habitat[species[s]] + 
                beta_mass * mass[species[s]] +
                beta_pitch * pitch[species[s]] +
                log_TAU;
                
  log_tau_forest_offroad[s] = log_tau_open_offroad[s] + beta_forest[species[s]];
  log_tau_open_onroad[s] = log_tau_open_offroad[s] + beta_roadside[species[s]];
  log_tau_forest_onroad[s] = log_tau_open_offroad[s] + beta_forest[species[s]] + beta_roadside[species[s]] + beta_interaction[species[s]];
  
  }
  
  
}
