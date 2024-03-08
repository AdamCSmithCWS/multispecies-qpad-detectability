####### Script Information ########################
# Adam C. Smith
# Multi-species QPAD Detectability
# 1-prepare-distance-data_project_1species.R
# Last Updated March 2024
# modification of the multi-species data prep script
# output is a Stan data list for a single species with project information

####### Import Libraries and External Files #######

library(ape)
library(Matrix)
# library(plyr)
# library(magrittr)
library(tidyverse)

#source("src/functions/generate-phylo-corr.R")

####### Read Data #################################

dist_count_matrix <- readRDS("data/raw/dist_count_matrix_project.rds")
load("data/raw/dist_design.rda")
# load(file = "data/generated/corr_matrix_predict.rda")
# binomial <- read.csv("data/generated/binomial_names.csv")
# traits <- read.csv("data/raw/traits.csv")

####### Wrangle Data for Modelling ################

# Change infinite values to some other very large number to avoid Stan issues
dist_design <- do.call(data.frame,lapply(dist_design, function(x) replace(x, is.infinite(x),450)))

# Most of this code adopted from Edwards et al. 2022

# Drop method I
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

species_sel <- "AMRO"

dist_count_matrix <- dist_count_matrix %>% 
  rename_with(., ~ paste0("abund_per_band_",.x),
              matches("^[[:digit:]]")) %>%
  filter(Species == species_sel) %>% 
  select(-c(Species,proj1,proj2)) %>% 
  #mutate(proj = Distance_Method) %>% 
  relocate(Sample_ID,proj,Distance_Method)
  


full_df <- dist_count_matrix %>% 
  inner_join(.,dist_design,
             by = c("Distance_Method" = "Method")) %>% 
  filter(n_bands > 1)


n_c_proj <- full_df %>% 
  group_by(proj) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 100)

full_df <- full_df %>% 
  filter(proj %in% n_c_proj$proj)


pr_list1 <- full_df %>% 
  select(proj,Distance_Method) %>%
  distinct()

n_bands_full <- pr_list1 %>% 
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
  
  p <- as.character(n_bands_full[i,"proj"])
  m <- as.character(n_bands_full[i,"Distance_Method"])
  tmp1 <- full_df %>% 
    filter(proj == p,
           Distance_Method == m)
  
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
  n_bands_full[i,"flag_no_obs_mid"] <- ifelse(sum_abunds > 100 & any(p_dist < 0.01),TRUE,FALSE)
  n_bands_full[i,"total_obs"] <- sum(abunds,na.rm = TRUE)
  
}

n_bands_full <- n_bands_full %>% 
  arrange(-flag_no_obs_far,-total_obs)

proj_problem <- n_bands_full %>% 
  filter(flag_no_obs_near | flag_no_obs_far | flag_no_obs_mid) %>% 
  distinct()

write_csv(proj_problem,"possible_prob_projs_March5.csv")

full_df <- full_df %>% 
  filter(!proj %in% proj_problem$proj)






n_c_proj3 <- full_df %>% 
  group_by(proj) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 100)

full_df <- full_df %>% 
  filter(proj %in% n_c_proj3$proj)

landcover <- readRDS("data/landcover_covariates.rds") %>% 
  select(Sample_ID,roaddist,ForestOnly_5x5) %>% 
  mutate(roadside = ifelse(roaddist < 30,1,0),
         forest = ifelse(ForestOnly_5x5 > 12,1,0)) %>% 
  select(-c(roaddist,ForestOnly_5x5))

full_df <- full_df %>% 
  inner_join(., landcover,
             by = "Sample_ID")

full_df[,"project"] <- as.integer(factor(full_df$proj))

stan_data <- list(
  n_projects = max(full_df$project),
  n_samples = nrow(full_df),
  max_intervals = max_bands,
  grainsize = 1,
  
  project = full_df$project,
  forest = full_df$forest,
  roadside = full_df$roadside,
  bands_per_sample = full_df$n_bands,
  abund_per_band = as.matrix(full_df[,paste0("abund_per_band_",1:max_bands)]),
  max_dist = as.matrix(full_df[,paste0("max_dist_",1:max_bands)])
  )


####### Output ####################################
save(stan_data, file = "data/generated/distance_stan_data_project_cov_1species.rda")



library(cmdstanr)
source("src/functions/generate-distance-inits.R")

# inits <- generate_distance_inits(n_chains = n_chains,
#                                  sp_list = setdiff(as.vector(dis_data$sp_all), pred_drops),
#                                  napops_skip = NULL,
#                                  param = "cp")
# 


mean_dist <- vector("numeric",stan_data$n_projects)
n_bands <- data.frame(proj = 1:stan_data$n_projects,
                      n_band = 0,
                      flag_no_obs_near = NA,
                      flag_no_obs_far = NA)

for(p in 1:length(mean_dist)){
  wsel <- which(stan_data$project == p)
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
  
}

inits <- vector("list",4)
for(i in 1:length(inits)){
  inits[[i]] <- list(log_tau_raw = rnorm(stan_data$n_projects,
                                         mean = log(mean_dist/2),
                                         sd = 0.01),
                     log_TAU = rnorm(1,0,sd = 0.01),
                     sd_log_tau = abs(rnorm(1,1,0.01)),
                     beta_forest = 0,
                     beta_roadside = 0,
                     beta_interaction = 0)
}


mod_file <- "models/distance_project_cov_1species.stan"

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


  csv_files <- paste0("output/temp-",c(1:4),".csv")
  stanfit <- cmdstanr::as_cmdstan_fit(files = csv_files)


summ <- stanfit$summary()
saveRDS(stanfit,"output/stanfit_AMRO_method.rds")


edrs <- summ %>% 
  filter(grepl("log_tau[",variable,fixed = TRUE)) %>% 
  mutate(edr = exp(mean)*100,
         edr_lci = exp(q5)*100,
         edr_uci = exp(q95)*100,
         proj = row_number()) %>% 
  inner_join(.,pr_list,
             by = "proj")

EDR <- summ %>% 
  filter(variable == "log_TAU")%>% 
  mutate(edr = exp(mean)*100,
         edr_lci = exp(q5)*100,
         edr_uci = exp(q95)*100,
         project = "Hyperparameter")

edrs <- edrs %>% 
  bind_rows(.,EDR) %>% 
  mutate(proj = as.character(proj))

edr_plot <- ggplot(data = edrs,
                   aes(x = project,y = edr))+
  geom_errorbar(aes(ymin = edr_lci,
                    ymax = edr_uci),
                width = 0,
                alpha = 0.5)+
  geom_point(aes(colour = n_band))+
  geom_errorbar(data = EDR,
                aes(x = project,y = edr,
                    ymin = edr_lci,
                    ymax = edr_uci),
                width = 0,
                colour = "darkorange",
                inherit.aes = FALSE)+
  geom_point(data = EDR,
             colour = "darkorange",
             aes(x = project,y = edr),
             inherit.aes = FALSE)+
  scale_colour_viridis_c()+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  #theme(text = element_text(size = 6))+
  ylab("EDR (m)")+
  xlab("Distance Method")+
  coord_flip()


edr_plot

pdf(paste0(species_sel,"_method_level_estimates.pdf"),height = 8,
    width = 8.5)
print(edr_plot)
dev.off()

