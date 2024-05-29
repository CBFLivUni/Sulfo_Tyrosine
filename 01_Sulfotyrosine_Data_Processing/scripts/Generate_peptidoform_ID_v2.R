######### load libraries ##########
library(tidyverse)
library(stringr)

# this function requires a data frame with labelled sample identifier rownames 
# the data frame must contain peptide modification  and peptide sequence data 
# by default the expected columns are named 'mod_peptide' and 'peptide' respectively
# the assigns ID for all possible peptidofroms (not strict -e.g. position of PTMs not taken into account)

df <- filtered_data[[1]]

# helper function to generate a single peptidoform IF
generate_peptidoform_id <- function(peptide, modifications) {
  # Normalise the peptide sequence by removing the PTMs to treat them separately
  normalized_peptide <- gsub("\\[[0-9]+\\]", "", peptide)
  ptms <- str_extract_all(peptide, "[A-Za-z]\\[[0-9]+\\]") %>% unlist()
  
  # Count and sort the PTMs
  ptm_counts <- table(ptms)
  sorted_ptms <- sort(names(ptm_counts))
  
  # Initialise the peptidoform ID with the normalised peptide sequence
  peptidoform_id <- normalized_peptide
  
  # Loop through each sorted modification and append its count to the peptidoform ID
  for(mod in sorted_ptms) {
    count <- ptm_counts[mod]
    peptidoform_id <- paste(peptidoform_id, mod, count, sep = "_")
  }
  
  return(peptidoform_id)
}


Generate_Peptidoform_ID <- function(df, outname = NA) {
  if (is.na(outname)) {
    df_name <- deparse(substitute(df))
  } else {
    df_name <- outname
  }
  
  print("Getting list of unique PTM modifications present in the dataset...")
  modifications_withamino <- str_extract_all(df$mod_peptide, "[A-Za-z]\\[[0-9]+\\]") %>%
    unlist() %>%
    unique() %>%
    sort()
  
 
  print("Generating peptidoform identifiers...")
  # generate all IDs in the df
  df$peptidoform_id <- sapply(df$mod_peptide, function(peptide) generate_peptidoform_id(peptide, modifications_withamino))
  # save all peptidoform IDs to csv
  df$peptidoform_id[1]
  filename <- paste0("out/", df_name, "_peptidoforms.csv")
  write.csv(df, file = filename, row.names = FALSE)
  print(paste("Data with added peptidoform IDs written to", filename))
  
  return(df)
}


