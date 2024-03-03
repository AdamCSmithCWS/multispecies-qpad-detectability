functions {
  real partial_sum_lpmf(array[, ] int slice_abund_per_band, // modifications to avoid Stan warnings at compile
                        int start, int end,
                        int max_intervals,
                        array[] int bands_per_sample,
                        array[, ] real max_dist,
                        array[] int project,
                        vector log_tau)
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
            Pi[Pi_index,k] = ((1 - exp(-(max_dist[i,k]^2 / exp(log_tau[project[i]])^2))) - 
                                (1 - exp(-(max_dist[i,k - 1]^2 / exp(log_tau[project[i]])^2)))) / 
              (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[project[i]])^2)));
          }else{
            Pi[Pi_index,k] = (1 - exp(-(max_dist[i,k]^2 / exp(log_tau[project[i]])^2))) /
              (1 - exp(-(max_dist[i,bands_per_sample[i]]^2 / exp(log_tau[project[i]])^2)));
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
  int<lower = 1> n_projects;           // total number of specise
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int<lower = 1> grainsize;           // grainsize for reduce_sum() function
  array[n_samples] int project;             // project being considered for each sample
  array[n_samples, max_intervals] int abund_per_band;// abundance in distance band k for sample i
  array[n_samples] int bands_per_sample; // number of distance bands for sample i
  array[n_samples, max_intervals] real max_dist; // max distance for distance band k
}

parameters {
  vector[n_projects] log_tau_raw;
  real log_TAU;
  real<lower=0> sd_log_tau;
}

transformed parameters{
  vector[n_projects] log_tau;
  
  log_tau = sd_log_tau*log_tau_raw;
  
}

model {
  log_tau_raw ~ std_normal();
  log_TAU ~ std_normal();
  sd_log_tau ~ std_normal();
  
  target += reduce_sum(partial_sum_lupmf,
                       abund_per_band,
                       grainsize,
                       max_intervals,
                       bands_per_sample,
                       max_dist,
                       project,
                       log_tau);
}
