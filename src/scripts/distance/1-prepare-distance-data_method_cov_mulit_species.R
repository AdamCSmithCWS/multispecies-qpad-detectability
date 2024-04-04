####### Script Information ########################
# Adam C. Smith
# Multi-species QPAD Detectability
# 1-prepare-distance-data_method_cov_multi_species.R
# Last Updated March 2024
# modification of the multi-species data prep script to prepare data for
# model that estimates full survey-specific covariate effects (model 5) while also fitting
# a varying intercept for each field method and 
# covariate effects for species-level variation (pitch, mass, migration, habitat association)

setwd("C:/Users/SmithAC/Documents/GitHub/multispecies-qpad-detectability")
####### Import Libraries and External Files #######

library(ape)
library(Matrix)
# library(plyr)
# library(magrittr)
library(tidyverse)
library(napops)
library(cmdstanr)

#source("src/functions/generate-phylo-corr.R")


# Compile distance observation data ---------------------------------------

  
load("data/raw/dist_count_matrix.rda")
load("data/raw/dist_design.rda")
# load(file = "data/generated/corr_matrix_predict.rda")
# binomial <- read.csv("data/generated/binomial_names.csv")

####### Wrangle Data for Modelling ################

# Change infinite values to some other very large number to avoid Stan issues
dist_design <- do.call(data.frame,lapply(dist_design, function(x) replace(x, is.infinite(x),450)))

# Most of this code modified from original (my apologies)
# Drop method I - weird method without consistently increasing distance bands
dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
dist_design <- dist_design[-which(dist_design$Method == "I"), ]

max_bands <- ncol(dist_design) - 2

cls <- paste0("X",1:max_bands)
for(j in 1:nrow(dist_design)){
  dist_design[j,"n_bands"] <- length(which(!is.na(  dist_design[j,cls])))
}

## replace the NAs with 0s to support Stan
dist_design <- dist_design %>% 
  mutate(X1 = ifelse(is.na(X1),0,X1),
         X2 = ifelse(is.na(X2),0,X2),
         X3 = ifelse(is.na(X3),0,X3),
         X4 = ifelse(is.na(X4),0,X4),
         X5 = ifelse(is.na(X5),0,X5),
         X6 = ifelse(is.na(X6),0,X6),
         X7 = ifelse(is.na(X7),0,X7),
         X8 = ifelse(is.na(X8),0,X8),
         X9 = ifelse(is.na(X9),0,X9),
         X10 = ifelse(is.na(X10),0,X10),
         X11 = ifelse(is.na(X11),0,X11),
         X12 = ifelse(is.na(X12),0,X12),
         X13 = ifelse(is.na(X13),0,X13)) %>% 
  rename_with(., ~gsub("X","max_dist_",.x),
              .cols = starts_with("X"))
dist_count_matrix <- dist_count_matrix %>% 
  rename_with(., ~ paste0("abund_per_band_",.x),
              matches("^[[:digit:]]")) %>%
  relocate(Sample_ID,Distance_Method)

full_df <- dist_count_matrix %>% 
  inner_join(.,dist_design,
             by = c("Distance_Method" = "Method")) %>% 
  filter(n_bands > 1)


highest <- full_df %>% 
  mutate(row_sums = rowSums(select(.,starts_with("abund_per")))) %>% 
  group_by(Species) %>% 
  summarise(max_total_count_sp = max(row_sums),
            q95_sp = quantile(row_sums,0.95))

full_df <- full_df %>% 
  mutate(row_sums = rowSums(select(.,starts_with("abund_per")))) %>% 
  left_join(.,highest,
            by = "Species") %>% 
  filter(.,
         (max_total_count_sp >= 20 & row_sums < q95_sp)|
           (max_total_count_sp < 20),
         Species != "RUHU")


# oddones <- full_dftt %>% 
#   filter(row_sums > 20)




  






# Compile survey covariates -----------------------------------------------



landcover <- readRDS("data/landcover_covariates.rds") %>% 
  select(Sample_ID,roaddist,ForestOnly_5x5) %>% 
  mutate(roadside = as.integer(ifelse(roaddist < 100,2,1)),
         #forest = as.integer(ifelse(ForestOnly_5x5 > 12,2,1)),
         forest = (round(10*ForestOnly_5x5/25))+1) %>% 
  select(-c(roaddist,ForestOnly_5x5))

full_df <- full_df %>% 
  inner_join(., landcover,
             by = "Sample_ID")






# compile traits information ----------------------------------------------

sp_w_data <- unique(full_df$Species)

traits <- read.csv("data/raw/traits.csv") %>% 
  rename(Species = Code) %>% 
  filter(Species %in% sp_w_data) %>% 
  mutate(species = as.integer(factor(Species)),
         mass = as.numeric(scale(log(Mass))),
         pitch = as.numeric(scale(Pitch)),
         migrant = Migrant,
         habitat = Habitat)

sp_w_traits <- traits %>% 
  select(Species,species,mass,pitch,migrant,habitat)

full_df <- full_df %>% 
  inner_join(.,sp_w_traits,
             by = "Species")

n_species <- max(full_df$species)
n_sp_2 <- length(unique(full_df$Species))
n_sp_3 <- max(traits$species)
if(!all(n_species == c(n_sp_2,n_sp_3))){
  stop("Something went wrong with species list join")
}

traits <- traits %>% 
  arrange(species)





# data summaries ----------------------------------------------------------



full_df[,"method"] <- as.integer(factor(full_df$Distance_Method))


method_list <- full_df %>% 
  group_by(Distance_Method,method) %>% 
  summarise(n_obs = n()) %>% 
  #distinct() %>% 
  left_join(.,dist_design,
            by = c("Distance_Method" = "Method")) %>% 
  select(Distance_Method,method,n_bands,n_obs,starts_with("max_dist")) 



if(nrow(method_list) < 2){next}
for(j in 1:nrow(method_list)){
  bb <- as.integer(method_list[j,"n_bands"])
  method_list[j,"method_string"] <- paste(as.character(method_list[j,paste0("max_dist_",1:bb)]),
                                          collapse = ":")
} 

method_list <- method_list %>% 
  select(-starts_with("max_dist_"))








# create Stan data --------------------------------------------------------


stan_data <- list(
  n_methods = max(full_df$method),
  n_samples = nrow(full_df),
  n_species = n_species,
  max_intervals = max_bands,
  grainsize = 1,
  n_forests = max(full_df$forest),
  
  method = full_df$method,
  forest = full_df$forest,
  roadside = full_df$roadside,
  bands_per_sample = full_df$n_bands,
  species = full_df$species,
  abund_per_band = as.matrix(full_df[,paste0("abund_per_band_",1:max_bands)]),
  max_dist = as.matrix(full_df[,paste0("max_dist_",1:max_bands)])/100,
  
  #mig_strat = traits$migrant,
  #habitat = traits$habitat,
  mass = traits$mass,
  pitch = traits$pitch
  )


####### Output ####################################
save(stan_data, file = paste0("data/generated/distance_stan_data_method_cov_multi_species.rda"))



#source("src/functions/generate-distance-inits.R")

# inits <- generate_distance_inits(n_chains = n_chains,
#                                  sp_list = setdiff(as.vector(dis_data$sp_all), pred_drops),
#                                  napops_skip = NULL,
#                                  param = "cp")
# 


mean_dist <- vector("numeric",stan_data$n_methods)

for(p in 1:length(mean_dist)){
  wsel <- which(stan_data$method == p)
  # nb <- matrix(nrow = stan_data$bands_per_sample[wsel],
  #              ncol = max(stan_data$bands_per_sample[wsel]))
  dists <- stan_data$max_dist[wsel,]
  abunds <- stan_data$abund_per_band[wsel,]
  mean_dist[p] <- sum(colSums(dists*abunds),na.rm = TRUE)/sum(abunds,na.rm = TRUE)
  
  max_dist <- max(stan_data$bands_per_sample[wsel])
  
  tmp <- colSums(abunds)[1:max_dist]
  
  method_list[which(method_list$method == p),"total_counts"] <- sum(tmp)
}


mean_dist_sp <- vector("numeric",stan_data$n_species)

for(p in 1:stan_data$n_species){
  j <- which(traits$species == p)
  eng <- traits[j,"Species"]
  
  tmp <- napops::edr(species = eng, model = 1,
                     road = FALSE,
                     forest = 0.5)$EDR_est/100
  
  if(tmp < 0.25){
  mean_dist_sp[p] <- NA
  }else{
    mean_dist_sp[p] <- tmp
  }
  rm(tmp)
}

miss <- which(is.na(mean_dist_sp))
mean_dist_sp[miss] <- mean(mean_dist_sp,na.rm = TRUE)


inits <- vector("list",4)
for(i in 1:length(inits)){
  inits[[i]] <- list(log_tau_method_raw = rnorm(stan_data$n_methods,
                                         mean = log(mean_dist/2),
                                         sd = 0.01),
                     log_TAU = rnorm(1,0,sd = 0.01),
                     sd_log_tau_method = abs(rnorm(1,1,0.01)),
                     sd_log_tau_species = abs(rnorm(1,1,0.01)),
                     beta_forest = 0,
                     beta_roadside = 0,
                     beta_interaction = 0,
                     #beta_mig_strat = 0,
                     #beta_habitat = 0,
                     beta_pitch = -0.02,
                     beta_mass = 0.02,
                     log_tau_species_raw = rnorm(stan_data$n_species,
                                                 mean = log(mean_dist_sp)/5,
                                                 sd = 0.01))
}


mod_file <- "models/distance_method_cov_multi_species.stan"

model <- cmdstan_model(mod_file,#,stanc_options = list("Oexperimental"),
                       cpp_options = list(stan_threads = TRUE))

stanfit <- model$sample(
  data=stan_data,
  refresh=100,
  # iter_sampling=3000,
  # iter_warmup=1000,
  #init = 0.01,
  init = inits,
  threads_per_chain = 8,
  parallel_chains = 4,
  output_dir = "output",
  output_basename = "cov_second_temp")#,
  #pars = parms,
  #adapt_delta = 0.8,
  #max_treedepth = 12,
  #seed = 123)


  csv_files <- paste0("output/cov_third_temp-",c(1:4),".csv")
  stanfit <- cmdstanr::as_cmdstan_fit(files = csv_files)

summ <- stanfit$summary()
saveRDS(stanfit,paste0("output/stanfit_multi_cov_method.rds"))

stanfit <- readRDS(paste0("output/stanfit_multi_cov_method.rds"))
stanarray <- cmdstanr::as_mcmc.list(stanfit)
shinystan::launch_shinystan(as.shinystan(stanarray))

# 
# source("c:/github/handy_functions/extract_stan_dimensions.R")
# edrs_only <- summ %>% 
#   filter(grepl("log_tau[",variable,fixed = TRUE)) %>% 
#   mutate(edr = exp(mean)*100,
#          edr_lci = exp(q5)*100,
#          edr_uci = exp(q95)*100) %>% 
#   drop_na() %>% 
#   mutate(method = dim_ext(dim = 1,var = "log_tau",variable),
#          forest = dim_ext(dim = 2,var = "log_tau",variable),
#          roadside = dim_ext(dim = 3,var = "log_tau",variable)) %>% 
#   inner_join(.,method_list,
#              by = "method")
# 
# EDR <- summ %>% 
#   filter(variable == "log_TAU")%>% 
#   mutate(edr = exp(mean)*100,
#          edr_lci = exp(q5)*100,
#          edr_uci = exp(q95)*100,
#          Distance_Method = "Hyperparameter",
#          forest = 1,
#          roadside = 1)
# 
# EDR_covs <- summ %>% 
#   filter(variable %in% c("log_tau_open_offroad",
#                          "log_tau_open_onroad",
#                          "log_tau_forest_offroad",
#                          "log_tau_forest_onroad"))%>% 
#   mutate(edr = exp(mean)*100,
#          edr_lci = exp(q5)*100,
#          edr_uci = exp(q95)*100,
#          Distance_Method = "Hyperparameter",
#          forest = ifelse(grepl("forest",variable),stan_data$n_forests,1),
#          roadside = ifelse(grepl("onroad",variable),2,1))
# 
# 
# edrs <- edrs_only %>% 
#   filter(forest %in% c(1,stan_data$n_forests)) %>% 
#   bind_rows(.,EDR) %>% 
#   bind_rows(.,EDR_covs) %>% 
#   mutate(method = as.character(method),
#          forest = ifelse(forest == stan_data$n_forests,"forest","open"),
#          roadside = ifelse(roadside == 2,"roadside","offroad"),
#          obs_location = paste(roadside,forest,sep = "-"),
#          Distance_Method2 = as.character(ifelse(Distance_Method == "Hyperparameter",
#                                                 Distance_Method,method_string)),
#          total_counts = ifelse(is.na(total_counts),
#                                exp(mean(log(total_counts),na.rm = TRUE)),
#                                total_counts))
# 
# edr_plot <- ggplot(data = edrs,
#                    aes(x = Distance_Method2,y = edr,
#                        colour = obs_location))+
#   geom_errorbar(data = EDR,
#                 aes(x = Distance_Method,y = edr,
#                     ymin = edr_lci,
#                     ymax = edr_uci),
#                 width = 0,
#                 colour = "darkorange",
#                 inherit.aes = FALSE)+
#   geom_point(data = EDR,
#              colour = "darkorange",
#              aes(x = Distance_Method,y = edr),
#              inherit.aes = FALSE)+
#   geom_errorbar(aes(ymin = edr_lci,
#                     ymax = edr_uci),
#                 width = 0,
#                 alpha = 0.5,
#                 position = position_dodge(width = 1))+
#   geom_point(position = position_dodge(width = 1),
#              aes(size = total_counts))+
#   scale_colour_viridis_d(guide = guide_legend(reverse = TRUE))+
#   scale_y_continuous(limits = c(0,NA))+
#   scale_size_continuous(trans = "log10",
#                         range = c(0.2,4))+
#   theme_bw()+
#   theme(text = element_text(size = 6))+
#   ylab("EDR (m)")+
#   xlab("Distance Method")+
#   coord_flip()
# 
# 
# edr_plot
# 
# pdf(paste0(species_sel,"_method_cov_level_estimates.pdf"),height = 11,
#     width = 8.5)
# print(edr_plot)
# dev.off()
# 
# 
# }
