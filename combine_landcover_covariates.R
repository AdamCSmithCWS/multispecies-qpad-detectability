## compile landcover covariates
library(tidyverse)

cov_dir <- "C:/github/na-pops-covariates/landcover"

file_list <- dir(cov_dir)

covs_out <- NULL
for(ff in file_list){
  tmp <- read_csv(paste0(cov_dir,"/",ff),
                  col_types = c("cddddddddddddd"))
  
  covs_out <- bind_rows(covs_out,tmp)
}

covs_out <- covs_out %>% 
  filter(!is.na(Needleleaf_3x3))

saveRDS(covs_out,"data/landcover_covariates.rds")
