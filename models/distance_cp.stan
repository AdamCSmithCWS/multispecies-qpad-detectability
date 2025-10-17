functions {
  real partial_sum_lpmf(array[, ] int slice_abund_per_band, // modifications to avoid Stan warnings at compile
                        int start, int end,
                        int max_intervals,
                        array[] int bands_per_sample,
                        array[, ] real max_dist,
                        array[] int species,
                        vector log_tau)
  {
    real lp = 0;
    int Pi_size = end - start + 1;
    int Pi_index = 1;
    matrix[Pi_size, max_intervals] Pi = rep_matrix(0, Pi_size, max_intervals);
    
    for (i in start:end)
    {
      for (k in 1:(bands_per_sample[i]-1)) // what if the final band was usedas the constraint? more effecient?
      {
        if(k > 1){
        Pi[Pi_index,k] = ((1 - exp(-(max_dist[i,k]^2 / exp(log_tau[species[i]])^2))) - 
                      (1 - exp(-(max_dist[i,k - 1]^2 / exp(log_tau[species[i]])^2)))) / 
        (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[species[i]])^2)));
        }else{
        Pi[Pi_index,k] = (1 - exp(-(max_dist[i,k]^2 / exp(log_tau[species[i]])^2))) /
        (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[species[i]])^2)));
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
  int<lower = 1> n_species;           // total number of specise
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> grainsize;           // grainsize for reduce_sum() function
  array[n_samples] int species;             // species being considered for each sample
  array[n_samples, max_intervals] int abund_per_band;// abundance in distance band k for sample i
  array[n_samples] int bands_per_sample; // number of distance bands for sample i
  array[n_samples, max_intervals] real max_dist; // max distance for distance band k
  int<lower = 1> n_mig_strat;        //total number of migration strategies
  array[n_species] int mig_strat;        //migration strategy for each species
  int<lower = 1> n_habitat;        //total number of habitat preferences
  array[n_species] int habitat;        //habitat preference for each species
  array[n_species] real mass;  //log mass of species
  array[n_species] real pitch; //song pitch of species
}

parameters {
  real intercept;
  row_vector[n_mig_strat] mu_mig_strat;
  row_vector[n_habitat] mu_habitat;
  real beta_mass;
  real beta_pitch;
  vector[n_species] log_tau;
  real<lower = 0> sigma;
}

model {
  
  for (sp in 1:n_species)
  {
      log_tau[sp] ~ normal((intercept + mu_mig_strat[mig_strat[sp]] +
                       mu_habitat[habitat[sp]] + //should this be the habitat at the survey?
                       beta_mass * mass[sp] +
                       beta_pitch * pitch[sp]), sigma);
  }

  intercept ~ normal(0.05, 0.1);
  mu_mig_strat ~ normal(0,0.1);
  mu_habitat ~ normal(0,0.1);
  
  beta_mass ~ normal(0,0.25);
  beta_pitch ~ normal(0,0.25);
  
  // beta_mass ~ normal(0.01,0.005);
  // beta_pitch ~ normal(-0.01,0.005);
  
  sigma ~ exponential(5);

  target += reduce_sum(partial_sum_lupmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_dist,
                       species,
                       log_tau);
}
