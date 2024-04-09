library(tidyverse)
library(janitor)
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/09_USI_generation/scripts")
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/09_USI_generation/out/of_interest_strict/"

# old USIs generated need modifying
files <- list.files(data_dir, full.names = TRUE)

# function to check if a string contains both letters and numbers - 
# as it would in the case of a PTM code
contains_letters_and_numbers <- function(s) {
  grepl("[A-Za-z]", s) && grepl("[0-9]", s)
}

PTM_codes <- list.files(data_dir, full.names = FALSE) %>% 
  gsub("added_USIs_to", "", .)  %>% 
  gsub(".csv", "", .) %>%
  strsplit(., "_") %>%
  unlist() %>%
  unique() %>%
  Filter(contains_letters_and_numbers, .) %>%
  gsub("(\\d+)", "[\\1]", .) # surround numbers with [ ] 

# "N[115]" - 18 for H2O + 1 for deamidation
# "Y[243]" - 18 for H2O + 80 for phosphorylation
# "T[181]" - 18 for H2O + 80 for phosphorylation
# "M[147]" - 18 for H2O + 16 for oxidation
# "S[167]" - 18 for H2O + 80 for phosphorylation
# "n[230]"  - assume [TMT6plex] - add at start and end of peptidoform
# "n[145]" - assume [TMT6plex] - add at start and end of peptidoform



# read in a file touse as exampe;  skip first column as it was added because i didn't sue row.names= FALSE
example <- read.csv(files[[1]])[,2:30] %>% clean_names()


colnames_to_use <- colnames(example)

# create e,mpty df to store all data and set names to the same as the example
all_data <- data.frame(matrix(ncol = length(colnames_to_use), nrow = 0))
colnames(all_data) <- colnames_to_use
# test <- rbind(empty_df, example, example)

# aggregate all data we have of interret together ( bind all rows)

for (file in files) {
  
  temp_dat <- read.csv(file)[,2:30] %>% clean_names()
  all_data <- rbind(all_data, temp_dat)
}

# sanity check - do we have data for all 27 peptidoforms? 
unique(all_data$peptidoform_id)
# yes

