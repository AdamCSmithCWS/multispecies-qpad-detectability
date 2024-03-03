####### Script Information ########################
# Brandon P.M. Edwards
# Multi-species QPAD Detectability
# 1-prepare-distance-data.R
# Created August 2022
# Last Updated October 2023

####### Import Libraries and External Files #######

library(ape)
library(Matrix)
library(plyr)
library(magrittr)

source("src/functions/generate-phylo-corr.R")

####### Read Data #################################

dist_count_matrix <- readRDS("data/raw/dist_count_matrix_project.rds")
load("data/raw/dist_design.rda")
load(file = "data/generated/corr_matrix_predict.rda")
binomial <- read.csv("data/generated/binomial_names.csv")
traits <- read.csv("data/raw/traits.csv")

####### Wrangle Data for Modelling ################

# Change infinite values to some other very large number to avoid Stan issues
dist_design <- do.call(data.frame,lapply(dist_design, function(x) replace(x, is.infinite(x),450)))

# Most of this code adopted from Edwards et al. 2022

# Drop method I
dist_count_matrix <- dist_count_matrix[-which(dist_count_matrix$Distance_Method == "I"), ]
dist_design <- dist_design[-which(dist_design$Method == "I"), ]

max_bands <- ncol(dist_design) - 2
count_names <- c("Sample_ID", "Species", "Distance_Method",
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

# Create separate data frames
counts <- count_design[,c("Sample_ID", "Species", "Distance_Method", count_names[4:length(count_names)])]
design <- count_design[,c("Sample_ID", "Species", "Distance_Method", design_names[3:length(design_names)])]
names(design) <- names(counts)
col_names <- names(counts)[4:length(names(counts))]

# Change 0s to NA in counts table where appropriate based on design table
for (i in col_names)
{
  indices <- which(counts[, i] == 0)
  
  counts[indices, i] <- ifelse(is.na(design[indices, i]), NA, 0)
}

# Get species codes and order of species codes for the Prediction Matrix
species_pred <- gsub(pattern = "_",
                     replacement= " ",
                     x = rownames(corr_matrix_predict))
binomial_pred <- binomial[which(binomial$Scientific_BT %in% species_pred), ]
binomial_pred <- binomial_pred[match(species_pred, binomial_pred$Scientific_BT), ]
#' There is a non-zero chance that Rufous hummingbird is weirdly a problematic species
#' so we are going to get rid of it.
binomial_pred <- binomial_pred[-which(binomial_pred$Code == "RUHU"), ]
species_pred_code <- binomial_pred$Code

# Create subset of traits dataset for the prediction matrix
traits_pred <- traits[which(traits$Code %in% species_pred_code), ]
traits_pred <- traits_pred[match(species_pred_code, traits_pred$Code), ]

# List of input counts
Y_pred <- vector(mode = "list", length = length(species_pred_code))
names(Y_pred) <- species_pred_code

# List of input design
D_pred <- vector(mode = "list", length = length(species_pred_code))
names(D_pred) <- species_pred_code

# Species indicator list
sp_list_pred <- vector(mode = "list", length = length(species_pred_code))
names(sp_list_pred) <- species_pred_code

for (s in species_pred_code)
{
  n_counts <- nrow(counts[counts$Species == s, ])
  
  Y_pred[[s]] <- as.matrix(counts[counts$Species==s, col_names])
  D_pred[[s]] <- as.matrix(design[design$Species==s, col_names])
  sp_list_pred[[s]] <- data.frame(Species = rep(s, n_counts))
}

Y_pred <- Y_pred[lengths(Y_pred) != 0]; Y_pred <- do.call(rbind, Y_pred)

D_pred <- D_pred[lengths(D_pred) != 0]; D_pred <- do.call(rbind, D_pred)

sp_list_pred <- do.call(rbind, sp_list_pred)

#' Corresponds with "bands_per_sample" in distance.stan
dist_bands_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(!is.na(x))))

#' Total species abundance per sampling event.
#' I.e., this is the sum of Y_sik over k
#' Corresponds with "abund_per_sample" in distance.stan
total_abund_per_sample_pred <- unname(apply(Y_pred, 1, function(x) sum(x, na.rm = TRUE)))

#' Factored list of species
#' Corresponds with "species" in removal.stan
sp_pred_numeric <- data.frame(Species = species_pred_code,
                              num = seq(1, length(species_pred_code)))
sp_pred_numeric <- join(sp_list_pred, sp_pred_numeric, by= "Species")

#' Create vector of indices corresponding to species modelled by centred and non-
#' centred parameterizations
count_per_sp <- data.frame(table(sp_pred_numeric$num))
species_ncp <- as.numeric(count_per_sp[which(count_per_sp$Freq <= 1000), "Var1"])
species_cp <- as.numeric(setdiff(count_per_sp$Var1, species_ncp))
# Account for any species with zero counts, add them into non-centred species
species_ncp <- sort(c(species_ncp,
                 setdiff(as.numeric(count_per_sp$Var1), c(species_ncp, species_cp))))

#' Corresponds with "abund_per_band" in distance.stan
abundance_per_band_pred <- Y_pred
abundance_per_band_pred[is.na(abundance_per_band_pred)] <- 0

#' Corresponds with "max_dist" in distance.stan
max_dist_pred <- D_pred
max_dist_pred[is.na(max_dist_pred)] <- 0

n_samples_pred <- nrow(Y_pred)

n_species_pred <- nrow(binomial_pred)

max_intervals_pred <- ncol(Y_pred)

# a 1 corresponds with resident, a 2 corresponds with a migrant
mig_strat_pred <- traits_pred$Migrant + 1

# a 1 corresponds with open habitat, a 2 corresponds with closed habitat
habitat_pred <- traits_pred$Habitat + 1

# centre and scale mass
mass_pred <- scale(log(traits_pred$Mass))

# centre and scale pitch
pitch_pred <- scale(traits_pred$Pitch)

distance_stan_data <- list(n_samples = n_samples_pred,
                           n_species = n_species_pred,
                           species_cp = species_cp,
                           n_species_cp = length(species_cp),
                           species_ncp = species_ncp,
                           n_species_ncp = length(species_ncp),
                           max_intervals = max_intervals_pred,
                           species = sp_pred_numeric$num,
                           abund_per_band = abundance_per_band_pred,
                           bands_per_sample = dist_bands_per_sample_pred,
                           max_dist = max_dist_pred,
                           sp_list = sp_list_pred$Species,
                           n_mig_strat = max(mig_strat_pred),
                           mig_strat = mig_strat_pred,
                           n_habitat = max(habitat_pred),
                           habitat = habitat_pred,
                           mass = mass_pred,
                           pitch = pitch_pred,
                           sp_all = species_pred_code)

####### Output ####################################
save(distance_stan_data, file = "data/generated/distance_stan_data.rda")
