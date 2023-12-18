######### load libraries ##########
library(tidyverse)
library(stringr)

# this function requires a data frame with labelled sample identifier rownames 
# the data frame must contain peptide modification  and peptide sequence data 
# by default the expected columns are named 'mod_peptide' and 'peptide' respectively
# the default PTMs of interest include phosphorylation of serine, threonine, and tyrosine


FilterPTMs <- function(df, 
                       pepcols = NULL,
                       PTMs_of_interest = c("S[167]", "T[181]", "Y[243]")){ # cols to plot histograms for
  print("PTM filtering initiated.")
 
  # subset data if only certain columns are of interest
  if (!is.null(pepcols))
    {
      print("Subsetting data...")
      sub <- df[, pepcols]
  } else {
      sub <- df
    }
  
  print("Getting list of unique PTM modifications present in the dataset...")
  modifications_withamino <- str_extract_all(sub$mod_peptide, "[A-Za-z]\\[[0-9]+\\]") %>% 
    unlist() %>% # flatten list to vector 
    unique() %>% # get only unique ones
    sort() # sort them
  
  print("The current dataset includes the following PTMs:")
  print(modifications_withamino)
  
  # issue a warning if not all PTMs_of_interest are present in the dataset
  if (!all(PTMs_of_interest %in% modifications_withamino)) {
    warning("Not all PTMs_of_interest are present in the dataset. Check the PTMs_of_interest list.")
  }
 
  # issue a warning if there are PTMs present in the dataset not listed in PTMs_of_interest
  extra_modifications <- setdiff(modifications_withamino, PTMs_of_interest)
  if (length(extra_modifications) > 0) {
    warning("Potential loss of data! There are PTMs present in the dataset that are not listed in PTMs_of_interest. \n These will be excluded from the subset. \n Review the PTMs_of_interest to ensure all relevant data is included.")
    print("Currently selected PTMs of interest include:")
    print(PTMs_of_interest)
  }
  
  # Filter the subset to only include PTMs of interest
  # first, create a regular expression pattern that matches any PTM in PTMs_of_interest
  # the pattern should look like "S[167]|T[181]|Y[243]"
  # ptm_pattern <- paste(PTMs_of_interest, collapse = "|")
  # this actually does not work - we get no matches, we need to escape the brackets! 
  
  # modify the PTMs_of_interest to escape the square brackets when the string is parsed by grepl
  ## effectively we want to add two \ before each type of square bracket:
  escaped_PTMs_of_interest <- gsub("\\[", "\\\\[", PTMs_of_interest) 
  escaped_PTMs_of_interest <- gsub("\\]", "\\\\]", escaped_PTMs_of_interest)
  
  
  # now generate the regular expression pattern as before
  ptm_pattern <- paste(escaped_PTMs_of_interest, collapse = "|")
  # test that it works grepl(ptm_pattern, sub$mod_peptide[2])
  
  # Filter the subset data
  sub_filtered <- sub[grepl(ptm_pattern, sub$mod_peptide), ]
  return(sub_filtered)
}

# sub_filtered <- FilterPTMs(test_data)
