library(tidyverse)

load("data/raw/dist_count_matrix.rda")

dist_count_matrix1 <- dist_count_matrix %>% 
  slice_sample(.,n = 10000)

tmp <- as.character(dist_count_matrix1[1,1])

str_extract_all(tmp,".+(?=:)")

string_split_mine <- function(x,nn = 1){
out <- str_split_fixed(x,":", n = nn+1)[nn]
return(out)
}

string_split_mine(tmp,1)

dist_count_matrix1 <- dist_count_matrix %>% 
  rowwise() %>% 
  mutate(proj1 = string_split_mine(Sample_ID,1),
         proj2 = string_split_mine(Sample_ID,2),
         proj = paste(proj1,proj2,sep = ":"))

saveRDS(dist_count_matrix1,
        "data/raw/dist_count_matrix_project.rds")

nProj <- dist_count_matrix1 %>% 
  group_by(proj) %>% 
  summarise(nrows = n())
