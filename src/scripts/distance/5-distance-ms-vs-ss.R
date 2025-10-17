####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 5-distance-ms-vs-ss.R
# Created October 2023
# Last Updated April 2024

####### Import Libraries and External Files #######

library(cmdstanr)

source("src/functions/generate-distance-inits.R")

####### Load Data #################################

load("data/generated/distance_stan_data_cv.rda")

####### Set Constants #############################

# Stan settings
n_iter <- 2000
n_warmup <- 1000
n_chains <- 4
refresh <- 400
threads_per_chain <- 3

####### Data Wrangling ############################

ms_model <- cmdstan_model(stan_file = "models/distance_cp.stan",
                          cpp_options = list(stan_threads = TRUE))

ss_model <- cmdstan_model(stan_file = "models/distance_single.stan",
                          cpp_options = list(stan_threads = TRUE))

####### Run Full Models ###########################

inits <- generate_distance_inits(n_chains = n_chains,
                                 sp_list = distance_stan_data_cv$sp_all,
                                 napops_skip = NULL,
                                 param = "cp")

distance_stan_data_cv$sp_list <- NULL
distance_stan_data_cv$sp_all <- NULL
distance_stan_data_cv$grainsize <- 1

# Scale the maximum distances for computational ease
distance_stan_data_cv$max_dist <- distance_stan_data_cv$max_dist / 100

# Strip away center and scale attributes from mass and pitch for analysis
distance_stan_data_cv$mass <- distance_stan_data_cv$mass[,1]
distance_stan_data_cv$pitch <- distance_stan_data_cv$pitch[,1]

# Multi species model
stan_run_ms <- ms_model$sample(
  data = distance_stan_data_cv,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  init = inits
)
stan_run_ms$save_object(file = paste0("output/model_runs/distance_ms.RDS"))

summ <- stan_run_ms$summary()


# Single species model
distance_stan_data_cv$n_mig_strat <- NULL
distance_stan_data_cv$mig_strat <- NULL
distance_stan_data_cv$n_habitat <- NULL
distance_stan_data_cv$habitat <- NULL
distance_stan_data_cv$mass <- NULL
distance_stan_data_cv$pitch <- NULL

stan_run_ss <- ss_model$sample(
  data = distance_stan_data_cv,
  iter_warmup = n_warmup,
  iter_sampling = n_iter,
  chains = n_chains,
  parallel_chains = n_chains,
  refresh = refresh,
  threads_per_chain = threads_per_chain,
  init = inits
)
stan_run_ss$save_object(file = paste0("output/model_runs/distance_ss.RDS"))
