####### Script Information ########################
# Adam C. Smith
# Multi-species QPAD Detectability
# 1-prepare-distance-data_method_cov_1species.R
# Last Updated March 2024
# modification of the multi-species data prep script to prepare data for
# model that estimates full covariate effects (model 5) while also fitting
# a varying intercept for each field method
# output is a Stan data list for a single species with method information

####### Import Libraries and External Files #######

library(ape)
library(Matrix)
# library(plyr)
# library(magrittr)
library(tidyverse)

#source("src/functions/generate-phylo-corr.R")

# c("AMRO","BOBO","REVI",
#   "CORA","STJA","OSFL")

####### Read Data #################################
for(species_sel in c("CHRA","SEWR","BAIS",
                     "STGR","GRPC","RNEP")){
  
  
  
dist_count_matrix <- readRDS("data/raw/dist_count_matrix_project.rds")
load("data/raw/dist_design.rda")
# load(file = "data/generated/corr_matrix_predict.rda")
# binomial <- read.csv("data/generated/binomial_names.csv")
# traits <- read.csv("data/raw/traits.csv")

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


# Removing all but one species --------------------------------------------


dist_count_matrix <- dist_count_matrix %>% 
  rename_with(., ~ paste0("abund_per_band_",.x),
              matches("^[[:digit:]]")) %>%
  filter(Species == species_sel) %>% 
  select(-c(Species,proj1,proj2,proj)) %>% 
  #mutate(proj = Distance_Method) %>% 
  relocate(Sample_ID,Distance_Method)
  


full_df <- dist_count_matrix %>% 
  inner_join(.,dist_design,
             by = c("Distance_Method" = "Method")) %>% 
  filter(n_bands > 1)


n_c_Distance_Method <- full_df %>% 
  group_by(Distance_Method) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 10)

full_df <- full_df %>% 
  filter(Distance_Method %in% n_c_Distance_Method$Distance_Method)


method_list1 <- full_df %>% 
  select(Distance_Method) %>%
  distinct()

n_bands_full <- method_list1 %>% 
  mutate( n_band = 0,
                      flag_no_obs_near = NA,
                      flag_no_obs_far = NA,
          flag_no_obs_mid = NA,
          SumInt1 = NA,
          SumInt2 = NA,
          SumInt3 = NA,
          SumInt4 = NA,
          SumInt5 = NA,
          SumInt6 = NA,
          SumInt7 = NA,
          SumInt8 = NA,
          SumInt9 = NA,
          SumInt10 = NA,
          SumInt11 = NA,
          SumInt12 = NA,
          SumInt13 = NA,
          total_obs = 0) %>% 
  as.data.frame()

for(i in 1:nrow(n_bands_full)){
  
  #p <- as.character(n_bands_full[i,"proj"])
  m <- as.character(n_bands_full[i,"Distance_Method"])
  tmp1 <- full_df %>% 
    filter(Distance_Method == m)
  
  intrvls <- tmp1 %>% 
    select(starts_with("max_dist_")) %>% 
    distinct() %>% 
    unlist(use.names = FALSE) 
  
  n_intrvls <- length(which(intrvls > 0))
  
  abunds <- as.integer(colSums(tmp1[,paste0("abund_per_band_",1:n_intrvls)]))
  sum_abunds <- sum(abunds,na.rm = TRUE)
  p_dist <- abunds/sum_abunds
  if(length(p_dist) > 2){
    p_dist <- p_dist[-c(1,length(p_dist))]
  }

  n_bands_full[i,paste0("SumInt",1:n_intrvls)] <- abunds
  n_bands_full[i,"n_band"] <- n_intrvls
  n_bands_full[i,"flag_no_obs_near"] <- ifelse(abunds[1] == 0,TRUE,FALSE)
  n_bands_full[i,"flag_no_obs_far"] <- ifelse(abunds[n_intrvls] == 0,TRUE,FALSE)
  n_bands_full[i,"flag_no_obs_mid"] <- ifelse(sum_abunds > 100 & any(p_dist < 0.001),TRUE,FALSE)
  n_bands_full[i,"total_obs"] <- sum(abunds,na.rm = TRUE)
  
}

n_bands_full <- n_bands_full %>% 
  arrange(-flag_no_obs_far,-total_obs)

method_problem <- n_bands_full %>% 
  filter(flag_no_obs_near | flag_no_obs_far | flag_no_obs_mid) %>% 
  distinct()

write_csv(method_problem,"possible_prob_methods_March5.csv")

full_df <- full_df %>% 
  filter(!Distance_Method %in% method_problem$Distance_Method)

if(nrow(full_df) < 10){
  next
}





landcover <- readRDS("data/landcover_covariates.rds") %>% 
  select(Sample_ID,roaddist,ForestOnly_5x5) %>% 
  mutate(roadside = as.integer(ifelse(roaddist < 100,2,1)),
         #forest = as.integer(ifelse(ForestOnly_5x5 > 12,2,1)),
         forest = (round(10*ForestOnly_5x5/25))+1) %>% 
  select(-c(roaddist,ForestOnly_5x5))

full_df <- full_df %>% 
  inner_join(., landcover,
             by = "Sample_ID")

n_c_method3 <- full_df %>%
  group_by(Distance_Method) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 10)

full_df <- full_df %>% 
  filter(Distance_Method %in% n_c_method3$Distance_Method)

full_df[,"method"] <- as.integer(factor(full_df$Distance_Method))


method_list <- full_df %>% 
  select(Distance_Method,method) %>% 
  distinct() %>% 
  left_join(.,dist_design,
            by = c("Distance_Method" = "Method")) %>% 
  select(Distance_Method,method,n_bands,starts_with("max_dist")) 
 

if(nrow(method_list) < 2){next}
for(j in 1:nrow(method_list)){
  bb <- as.integer(method_list[j,"n_bands"])
  method_list[j,"method_string"] <- paste(as.character(method_list[j,paste0("max_dist_",1:bb)]),
                                          collapse = ":")
} 
  
method_list <- method_list %>% 
  select(-starts_with("max_dist_"))



stan_data <- list(
  n_methods = max(full_df$method),
  n_samples = nrow(full_df),
  max_intervals = max_bands,
  grainsize = 1,
  n_forests = max(full_df$forest),
  
  method = full_df$method,
  forest = full_df$forest,
  roadside = full_df$roadside,
  bands_per_sample = full_df$n_bands,
  abund_per_band = as.matrix(full_df[,paste0("abund_per_band_",1:max_bands)]),
  max_dist = as.matrix(full_df[,paste0("max_dist_",1:max_bands)])/100
  )


####### Output ####################################
save(stan_data, file = paste0("data/generated/distance_stan_data",species_sel,"_method_cov_1species.rda"))



library(cmdstanr)
#source("src/functions/generate-distance-inits.R")

# inits <- generate_distance_inits(n_chains = n_chains,
#                                  sp_list = setdiff(as.vector(dis_data$sp_all), pred_drops),
#                                  napops_skip = NULL,
#                                  param = "cp")
# 


mean_dist <- vector("numeric",stan_data$n_methods)
n_bands <- data.frame(method = 1:stan_data$n_methods,
                      n_band = 0,
                      flag_no_obs_near = NA,
                      flag_no_obs_far = NA,
                      total_counts = NA)

for(p in 1:length(mean_dist)){
  wsel <- which(stan_data$method == p)
  # nb <- matrix(nrow = stan_data$bands_per_sample[wsel],
  #              ncol = max(stan_data$bands_per_sample[wsel]))
  dists <- stan_data$max_dist[wsel,]
  abunds <- stan_data$abund_per_band[wsel,]
  mean_dist[p] <- sum(colSums(dists*abunds),na.rm = TRUE)/sum(abunds,na.rm = TRUE)
  
  max_dist <- max(stan_data$bands_per_sample[wsel])
  
  tmp <- colSums(abunds)[1:max_dist]
  
  n_bands[p,c(5:(max_dist +4))] <- colSums(abunds)[1:max_dist]
  n_bands[p,"n_band"] <- max_dist
  n_bands[p,"flag_no_obs_near"] <- ifelse(colSums(abunds)[1] == 0,TRUE,FALSE)
  n_bands[p,"flag_no_obs_far"] <- ifelse(colSums(abunds)[max_dist] == 0,TRUE,FALSE)
  n_bands[p,"flag_no_obs_far"] <- ifelse(colSums(abunds)[max_dist] == 0,TRUE,FALSE)
  method_list[which(method_list$method == p),"total_counts"] <- sum(tmp)
}

inits <- vector("list",4)
for(i in 1:length(inits)){
  inits[[i]] <- list(log_tau_raw = rnorm(stan_data$n_methods,
                                         mean = log(mean_dist/2),
                                         sd = 0.01),
                     log_TAU = rnorm(1,0,sd = 0.01),
                     sd_log_tau = abs(rnorm(1,1,0.01)),
                     beta_forest = 0,
                     beta_roadside = 0,
                     beta_interaction = 0)
}


mod_file <- "models/distance_method_cov_1species.stan"

model <- cmdstan_model(mod_file,#,stanc_options = list("Oexperimental"),
                       cpp_options = list(stan_threads = TRUE))

stanfit <- model$sample(
  data=stan_data,
  refresh=100,
  # iter_sampling=3000,
  # iter_warmup=1000,
  #init = 0.01,
  init = inits,
  threads_per_chain = 3,
  parallel_chains = 4,
  output_dir = "output",
  output_basename = "cov_temp")#,
  #pars = parms,
  #adapt_delta = 0.8,
  #max_treedepth = 12,
  #seed = 123)


  # csv_files <- paste0("output/temp-",c(1:4),".csv")
  # stanfit <- cmdstanr::as_cmdstan_fit(files = csv_files)


summ <- stanfit$summary()
saveRDS(stanfit,paste0("output/stanfit_",species_sel,"_cov_method.rds"))



source("c:/github/handy_functions/extract_stan_dimensions.R")
edrs_only <- summ %>% 
  filter(grepl("log_tau[",variable,fixed = TRUE)) %>% 
  mutate(edr = exp(mean)*100,
         edr_lci = exp(q5)*100,
         edr_uci = exp(q95)*100) %>% 
  drop_na() %>% 
  mutate(method = dim_ext(dim = 1,var = "log_tau",variable),
         forest = dim_ext(dim = 2,var = "log_tau",variable),
         roadside = dim_ext(dim = 3,var = "log_tau",variable)) %>% 
  inner_join(.,method_list,
             by = "method")

EDR <- summ %>% 
  filter(variable == "log_TAU")%>% 
  mutate(edr = exp(mean)*100,
         edr_lci = exp(q5)*100,
         edr_uci = exp(q95)*100,
         Distance_Method = "Hyperparameter",
         forest = 1,
         roadside = 1)

EDR_covs <- summ %>% 
  filter(variable %in% c("log_tau_open_offroad",
                         "log_tau_open_onroad",
                         "log_tau_forest_offroad",
                         "log_tau_forest_onroad"))%>% 
  mutate(edr = exp(mean)*100,
         edr_lci = exp(q5)*100,
         edr_uci = exp(q95)*100,
         Distance_Method = "Hyperparameter",
         forest = ifelse(grepl("forest",variable),stan_data$n_forests,1),
         roadside = ifelse(grepl("onroad",variable),2,1))


edrs <- edrs_only %>% 
  filter(forest %in% c(1,stan_data$n_forests)) %>% 
  bind_rows(.,EDR) %>% 
  bind_rows(.,EDR_covs) %>% 
  mutate(method = as.character(method),
         forest = ifelse(forest == stan_data$n_forests,"forest","open"),
         roadside = ifelse(roadside == 2,"roadside","offroad"),
         obs_location = paste(roadside,forest,sep = "-"),
         Distance_Method2 = as.character(ifelse(Distance_Method == "Hyperparameter",
                                                Distance_Method,method_string)),
         total_counts = ifelse(is.na(total_counts),
                               exp(mean(log(total_counts),na.rm = TRUE)),
                               total_counts))

edr_plot <- ggplot(data = edrs,
                   aes(x = Distance_Method2,y = edr,
                       colour = obs_location))+
  geom_errorbar(data = EDR,
                aes(x = Distance_Method,y = edr,
                    ymin = edr_lci,
                    ymax = edr_uci),
                width = 0,
                colour = "darkorange",
                inherit.aes = FALSE)+
  geom_point(data = EDR,
             colour = "darkorange",
             aes(x = Distance_Method,y = edr),
             inherit.aes = FALSE)+
  geom_errorbar(aes(ymin = edr_lci,
                    ymax = edr_uci),
                width = 0,
                alpha = 0.5,
                position = position_dodge(width = 1))+
  geom_point(position = position_dodge(width = 1),
             aes(size = total_counts))+
  scale_colour_viridis_d(guide = guide_legend(reverse = TRUE))+
  scale_y_continuous(limits = c(0,NA))+
  scale_size_continuous(trans = "log10",
                        range = c(0.2,4))+
  theme_bw()+
  theme(text = element_text(size = 6))+
  ylab("EDR (m)")+
  xlab("Distance Method")+
  coord_flip()


edr_plot

pdf(paste0(species_sel,"_method_cov_level_estimates.pdf"),height = 11,
    width = 8.5)
print(edr_plot)
dev.off()


}
