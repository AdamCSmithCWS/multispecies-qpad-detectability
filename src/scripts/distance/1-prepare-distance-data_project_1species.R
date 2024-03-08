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


# Removing all but one species --------------------------------------------

species_sel <- "AMRO"

dist_count_matrix <- dist_count_matrix %>% 
  filter(Species == species_sel) %>% 
  select(-c(Species,proj1,proj2)) %>% 
  #mutate(proj = Distance_Method) %>% 
  relocate(Sample_ID,proj,Distance_Method)


n_c_proj <- dist_count_matrix %>% 
  group_by(proj) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 100)

dist_count_matrix <- dist_count_matrix %>% 
  filter(proj %in% n_c_proj$proj)


n_c_proj2 <- dist_count_matrix %>% 
  group_by(proj) %>% 
  summarise(n_counts = n()) 

count_names <- c("Sample_ID", "proj", "Distance_Method",
                 paste0(rep("Int", times = max_bands), 1:max_bands))
names(dist_count_matrix) <- count_names

design_names <- c("Distance_Method", "Max_Distance",
                  paste0(rep("Interval", times = max_bands), 1:max_bands))
names(dist_design) <- design_names

# Join data
count_design <- plyr::join(dist_count_matrix, dist_design,
                           by = "Distance_Method", type = "left")

# Filter out methods with only 1 interval
to_remove <- which(is.na(count_design$Interval2))
count_design <- count_design[-c(to_remove), ]


# exploring for projects with odd count distributions ---------------------

projs <- count_design %>% 
  select(proj,Distance_Method) %>% 
  distinct()

n_bands_full <- projs %>% 
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
          total_obs = 0)

for(i in 1:nrow(n_bands_full)){
  
  p <- n_bands_full[i,"proj"]
  m <- n_bands_full[i,"Distance_Method"]
  tmp1 <- count_design %>% 
    filter(proj == p,
           Distance_Method == m)
  
  intrvls <- tmp1 %>% 
    select(starts_with("Interval")) %>% 
    distinct() %>% 
    unlist(use.names = FALSE) 
  
  n_intrvls <- length(which(!is.na(intrvls)))
  
  abunds <- as.integer(colSums(tmp1[,paste0("Int",1:n_intrvls)]))
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

count_design <- count_design %>% 
  filter(!proj %in% proj_problem$proj)


n_c_proj3 <- count_design %>% 
  group_by(proj) %>% 
  summarise(n_counts = n()) %>% 
  filter(n_counts > 100)

count_design <- count_design %>% 
  filter(proj %in% n_c_proj3$proj)

landcover <- readRDS("data/landcover_covariates.rds") %>% 
  select(Sample_ID,roaddist,ForestOnly_5x5) %>% 
  mutate(roadside = ifelse(roaddist < 30,1,0),
         forest = ifelse(ForestOnly_5x5 > 12,1,0)) %>% 
  select(-c(roaddist,ForestOnly_5x5))

count_design <- count_design %>% 
  inner_join(., landcover,
             by = "Sample_ID")
# Create separate data frames
counts <- count_design[,c("Sample_ID", "proj", "Distance_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "proj", "Distance_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

covs <- count_design[,c("Sample_ID", "proj", "Distance_Method",
                        "roadside","forest")]
# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Get proj codes and order of proj codes for the Prediction Matrix
proj_pred_code <- unique(count_design$proj)

# # Create subset of traits dataset for the prediction matrix
# traits_pred <- traits[which(traits$Code %in% species_pred_code), ]
# traits_pred <- traits_pred[match(species_pred_code, traits_pred$Code), ]

# List of input counts
Y_pred <- vector(mode = "list", length = length(proj_pred_code))
names(Y_pred) <- proj_pred_code

# List of input design
D_pred <- vector(mode = "list", length = length(proj_pred_code))
names(D_pred) <- proj_pred_code

# project indicator list
pr_list_pred <- vector(mode = "list", length = length(proj_pred_code))
names(pr_list_pred) <- proj_pred_code

for (s in proj_pred_code)
{
  n_counts <- nrow(counts[counts$proj == s, ])
  
  Y_pred[[s]] <- as.matrix(counts[counts$proj==s, col_names])
  D_pred[[s]] <- as.matrix(design[design$proj==s, col_names])
  pr_list_pred[[s]] <- data.frame(proj = rep(s, n_counts))
}

Y_pred <- Y_pred[lengths(Y_pred) != 0]; Y_pred <- do.call(rbind, Y_pred)

D_pred <- D_pred[lengths(D_pred) != 0]; D_pred <- do.call(rbind, D_pred)

pr_list_pred <- do.call(rbind, pr_list_pred)

#' Corresponds with "bands_per_sample" in distance.stan
dist_bands_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(!is.na(x))))

#' Total proj abundance per sampling event.
#' I.e., this is the sum of Y_sik over k
#' Corresponds with "abund_per_sample" in distance.stan
total_abund_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(x, na.rm = TRUE)))

#' Factored list of proj
#' Corresponds with "proj" in distance_project_1speices.stan
pr_pred_numeric <- data.frame(proj = proj_pred_code,
                              num = seq(1, length(proj_pred_code)))
pr_pred_numeric <- inner_join(pr_list_pred, pr_pred_numeric, by= "proj")



max_distances_proj <- dist_design %>% 
  rowwise() %>% 
  select(Distance_Method,Max_Distance) %>% 
  distinct() %>% 
  inner_join(.,n_bands_full,
            by = "Distance_Method") %>% 
  select(Distance_Method, n_band, Max_Distance) 

#' Create vector of indices corresponding to species modelled by centred and non-
#' centred parameterizations
count_per_pr <- data.frame(table(pr_pred_numeric$num))
proj_ncp <- as.numeric(count_per_pr[which(count_per_pr$Freq <= 1000), "Var1"])
proj_cp <- as.numeric(setdiff(count_per_pr$Var1, proj_ncp))
# # Account for any proj with zero counts, add them into non-centred proj
# species_ncp <- sort(c(species_ncp,
#                  setdiff(as.numeric(count_per_sp$Var1), c(species_ncp, species_cp))))

pr_list <- pr_pred_numeric %>% 
  distinct() %>% 
  bind_cols(count_per_pr) %>% 
  rename(project = proj,
         proj = num,
         n_counts = Freq) %>% 
  select(-Var1) %>% 
  left_join(.,max_distances_proj,
            by = c("project" = "Distance_Method"))
#' Corresponds with "abund_per_band" in distance.stan
abundance_per_band_pred <- Y_pred
abundance_per_band_pred[is.na(abundance_per_band_pred)] <- 0

#' Corresponds with "max_dist" in distance.stan
max_dist_pred <- D_pred
max_dist_pred[is.na(max_dist_pred)] <- 0

n_samples_pred <- nrow(Y_pred)

n_proj_pred <- max(pr_pred_numeric$num)

max_intervals_pred <- ncol(Y_pred)
# 
# # a 1 corresponds with resident, a 2 corresponds with a migrant
# mig_strat_pred <- traits_pred$Migrant + 1
# 
# # a 1 corresponds with open habitat, a 2 corresponds with closed habitat
# habitat_pred <- traits_pred$Habitat + 1
# 
# # centre and scale mass
# mass_pred <- scale(log(traits_pred$Mass))
# 
# # centre and scale pitch
# pitch_pred <- scale(traits_pred$Pitch)

distance_stan_data <- list(n_samples = n_samples_pred,
                           n_projects = n_proj_pred,
                           # species_cp = species_cp,
                           # n_species_cp = length(species_cp),
                           # species_ncp = species_ncp,
                           # n_species_ncp = length(species_ncp),
                           max_intervals = max_intervals_pred,
                           project = pr_pred_numeric$num,
                           forest = covs$forest,
                           roadside = covs$roadside,
                           abund_per_band = abundance_per_band_pred,
                           bands_per_sample = dist_bands_per_sample_pred,
                           max_dist = max_dist_pred/100,
                           #pr_list = pr_list_pred$proj,
                           # n_mig_strat = max(mig_strat_pred),
                           # mig_strat = mig_strat_pred,
                           # n_habitat = max(habitat_pred),
                           # habitat = habitat_pred,
                           # mass = mass_pred,
                           # pitch = pitch_pred,
                           grainsize = 1)

####### Output ####################################
save(distance_stan_data, file = "data/generated/distance_stan_data_project_cov_1species.rda")



library(cmdstanr)
source("src/functions/generate-distance-inits.R")

# inits <- generate_distance_inits(n_chains = n_chains,
#                                  sp_list = setdiff(as.vector(dis_data$sp_all), pred_drops),
#                                  napops_skip = NULL,
#                                  param = "cp")
# 


mean_dist <- vector("numeric",distance_stan_data$n_projects)
n_bands <- data.frame(proj = 1:distance_stan_data$n_projects,
                      n_band = 0,
                      flag_no_obs_near = NA,
                      flag_no_obs_far = NA)

for(p in 1:length(mean_dist)){
  wsel <- which(distance_stan_data$project == p)
  # nb <- matrix(nrow = distance_stan_data$bands_per_sample[wsel],
  #              ncol = max(distance_stan_data$bands_per_sample[wsel]))
  dists <- distance_stan_data$max_dist[wsel,]
  abunds <- distance_stan_data$abund_per_band[wsel,]
  mean_dist[p] <- sum(colSums(dists*abunds),na.rm = TRUE)/sum(abunds,na.rm = TRUE)
  
  max_dist <- max(distance_stan_data$bands_per_sample[wsel])
  
  tmp <- colSums(abunds)[1:max_dist]
  
  n_bands[p,c(5:(max_dist +4))] <- colSums(abunds)[1:max_dist]
  n_bands[p,"n_band"] <- max_dist
  n_bands[p,"flag_no_obs_near"] <- ifelse(colSums(abunds)[1] == 0,TRUE,FALSE)
  n_bands[p,"flag_no_obs_far"] <- ifelse(colSums(abunds)[max_dist] == 0,TRUE,FALSE)
  n_bands[p,"flag_no_obs_far"] <- ifelse(colSums(abunds)[max_dist] == 0,TRUE,FALSE)
  
}

inits <- vector("list",4)
for(i in 1:length(inits)){
  inits[[i]] <- list(log_tau_raw = rnorm(distance_stan_data$n_projects,
                                         mean = log(mean_dist/2),
                                         sd = 0.01),
                     log_TAU = rnorm(1,0,sd = 0.01),
                     sd_log_tau = abs(rnorm(1,1,0.01)))
}


mod_file <- "models/distance_project_1species.stan"

model <- cmdstan_model(mod_file,#,stanc_options = list("Oexperimental"),
                       cpp_options = list(stan_threads = TRUE))

stanfit <- model$sample(
  data=distance_stan_data,
  refresh=100,
  # iter_sampling=3000,
  # iter_warmup=1000,
  #init = 0.01,
  init = inits,
  threads_per_chain = 3,
  parallel_chains = 4,
  output_dir = "output",
  output_basename = "temp")#,
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

