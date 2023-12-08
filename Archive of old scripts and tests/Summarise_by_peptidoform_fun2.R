######### load libraries ##########
library(tidyverse)
library(stringr)
library(data.table)
library(stringi)
library(microbenchmark)
library(hash)
# this function requires a data frame with labelled sample identifier rownames 
# the data frame must contain peptide modification  and peptide sequence data 
# by default the expected columns are named 'mod_peptide' and 'peptide' respectively
# the function assembles the data for all possible peptidofroms (not strict -e.g. position of PTMs not taken into account)

Summarise_By_Peptidoform <- function(df, 
                         columns_to_summarise = c("calibrated_error") ,# must have cols for histogram plotting
                         pepcols = c("peptide", "mod_peptide"),
                         outname = NA){ # cols that contain peptide data
  
  
  if (is.na(outname)) {
    df_name <- deparse(substitute(df)) # store the df name to use for file namings later
    # for some reason calling this at a later stage results in an error when trying to write files
  } else {
    df_name <- outname
  }
  
  
  print("Getting list of unique PTM modifications present in the dataset...")
  modifications_withamino <- str_extract_all(df$mod_peptide, "[A-Za-z]\\[[0-9]+\\]") %>%
    unlist() %>% # flatten list to vector
    unique() %>% # get only unique ones
    sort() # sort them
  
  # define a function to use for PTM counting for each type of PTM present in the data
  count_PTMs <- function(peptide, PTMs_vector) {
    # for a vector that contains list of PTMs detected in the data
   sapply(PTMs_vector, function(ptm) str_count(peptide, fixed(ptm))) #count the number of times each PTM appears
  }
  
  print("Counting PTMs for each peptide... This may take a little while.")
  # # getting the matrix of counts for each PTM
  # PTM_counts <- apply(df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
  # # this part of the function is a bit slow! lets test it with parallel computing?
  # Setting up parallel processing
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  
  # Load required libraries on each worker node
  clusterEvalQ(cl, library(stringr)) # or library(stringi)
  
  # Export necessary objects and functions to each worker
  clusterExport(cl, varlist = c("count_PTMs", "modifications_withamino", "df"))
  
  # Perform the counting in parallel
  PTM_counts <- parApply(cl, df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
  
  # Stop the cluster
  stopCluster(cl)
  
  # transposing the matrix so rownames match rownames of data 
  PTM_counts_matrix <- t(PTM_counts)
  
    # now we can generate a unique identified of the peptidoform by collapsing these counts 
  # into a single column named modificaitn1_modification2_....modificationN which gives the 
  # counts as e.g. 1_0_0_0_3 
  
  # collapsing the matrix to a vector of strings
  PTM_counts_string <- apply(PTM_counts_matrix, 1, paste, collapse="_")
  # getting the scipher for which position correcponds to which PTM in case we need to decode in the future
  PTM_scipher <- paste(colnames(PTM_counts_matrix), collapse="_")
  
  
  # Add a column to the dataframe with PTM counts code
  df$PTM_counts <- PTM_counts_string
  # to keep the schipher, we populate an entire column with it
  n <- nrow(df)
  df$PTM_counts_scipher <- rep(PTM_scipher, n)
  
  # next we need to assemble peptidoforms by peptide sequence and PTM count
  print("Generating peptidoform identifiers...")
  
  peptidoform_ids <- paste(df$peptide, df$PTM_counts, sep = "_")
  # check if it works correctly
  # head(peptidoform_id)
  df$peptidoform_ids <- peptidoform_ids
  # write this dataframe as csv in case you need it later; name it after the df passed
  # onto the function
  # df_name <- deparse(substitute(df)) # this bugs out here, moved to beginnign of function
  filename <- paste0("out/", df_name, "_peptidoforms.csv")
  # write.csv(df, file = filename, row.names = FALSE)
  print(paste("Data with added peptidoform IDs written to", filename))
  
  
  ######### string approach for aggreagating within dataset ########
  # result is a named vector of strings
  # Convert shifts to a single string for each id
  # dict_vector <- tapply(df$calibrated_error, df$peptidoform_id, FUN = function(x) paste(x, collapse = ","))
  # return(dict_vector)
  
  
  dict2 <- split(df$calibrated_error, df$peptidoform_ids)
  return(dict2)
  
  
}

# old testing code


# print("summarising data by peptidoform...")
# 
# str(dict2)
# head(dict2)
# str(dict2[1])
# 
# dict_vector[1:6]
# str(dict_vector[1])
# # Aggregate shifts into a single string for each id
# dict3 <- aggregate(calibrated_error ~ peptidoform_id, data = df, FUN = function(x) paste(x, collapse = ","))
# head(dict3)
# str(dict3[1,1])
# 
# 
# # Convert the result to a named vector or list if needed
# dict_vector <- setNames(dict$shift, dict$id)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# filename <- paste0("out/", df_name, "_peptidoform_groups.csv")
# print("Writing peptidoform Group Data")
# write.csv(as.data.frame(peptidoform_groups), file = filename, row.names = FALSE)
# return(peptidoform_groups)
















# # run function to test if it works.
# out <- AssemblePSMs(df = sub_filtered)

############## code testing leftovers ############



# count_PTMs <- function(peptide, PTMs_vector) {
#   # for a vector that contains list of PTMs detected in the data
#   # for testing only: 
#   # print(peptide)
#   # print(PTMs_vector)
#   sapply(PTMs_vector, function(ptm) str_count(peptide, fixed(ptm))) #count the number of times each PTM appears
#   
# }

# test the function
# output <- count_PTMs (peptide = df$mod_peptide[1],
#             PTMs_vector = modifications_withamino)
# the output is a named vector of integers, with positions named by peptide 
# modification, and integers showing number of times each PTM is present

# using count_PTMs on each row, we get a matrix where each column corresponds to
# a row in the df, and each row is a PTM. We want to transpose to be able to 
# merge this with the actual data lets test it on a small subset:

# we want this to work on a filtered dataframe that only includes peptides with 
# potential true phospho sites. (e.g. the output of the FilterPTMs function)
# df <- sub_filtered
# test_df <- df[1:10,]
# 
# # getting the matrix of counts for each PTM
# PTM_counts <- apply(test_df, 1, function(row) count_PTMs(row["mod_peptide"], modifications_withamino))
# # transposing the matrix so rownames match ronames of data 
# PTM_counts_matrix <- t(PTM_counts)











# Group by peptidoform and summarize calibrated mass shifts
# try on one column
# peptidoform_groups <- df %>%
#   group_by(peptidoform_id) %>%
#   summarise(calibrated_error = list(calibrated_error))