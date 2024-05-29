# to get foreground of IDs in Swissrot
processProteinIDs <- function(data_frame, 
                              Swiss_Prot_Human , 
                              neXtProt_list , 
                              RefSeq_list) {
  
  # data_frame = data_for_ORA[[1]]
  # Split protein IDs and clean them
  foreground_id_list <- lapply(strsplit(data_frame$protein, ", "), function(x) unlist(strsplit(x, ", ")))
  clean_foreground_ID_list <- lapply(foreground_id_list, ProteinIDCleaningFunction)
  foreground_cleaned_IDs <- lapply(clean_foreground_ID_list, unique)
  
  # Map IDs to SwissProt entries if available
  foreground_cleaned_IDs <- lapply(foreground_cleaned_IDs, function(sublist) {
    process_sublist(sublist, Swiss_Prot_Human, neXtProt_list, RefSeq_list)
  })
  
  # Flatten the IDs and add them to the dataframe
  flattened_foreground_cleaned_IDs <- sapply(foreground_cleaned_IDs, function(x) paste(x, collapse = ","))
  data_frame$cleaned_protein_IDs <- flattened_foreground_cleaned_IDs
  
  # Filter rows to only include ones with SwissProt Entry hits
  data_frame_inSwissProt <- data_frame[data_frame$cleaned_protein_IDs %in% Swiss_Prot_Human$Entry,]
  
  return(data_frame_inSwissProt)
}

###### Function to get clean filenames, and round and extract bin boundaries ######
clean_round_extract <- function(name, index) {
  
  #separate the name into pars splitting at _
  parts <- unlist(strsplit(name, "_"))
  
  # Rounding the numbers to the fourth decimal place; lower bin boundary is the 2nd part, upper is the 3rd part
  lower <- round(as.numeric(parts[2]), 4)
  upper <- round(as.numeric(parts[3]), 4)
  
  # Storing the boundaries
  # index = 1
  lower_boundaries[length(lower_boundaries) + 1] <<- lower
  upper_boundaries[length(upper_boundaries) + 1] <<- upper
  
  # Reconstructing the cleaned name (and returning it)
  return(paste0("mz_", lower, "_", upper))
}


clean_data <- function(data) {
  for (i in 1:ncol(data)) {
    cleaned_col <- lapply(data[, i], function(x) {
      clean_string <- gsub("\\[|\\]|'", "", x)
      if (i == 5) {
        clean_string <- gsub(";", ",", clean_string)
      }
      elements <- unlist(strsplit(clean_string, ',\\s*'))
      return(elements)
    })
    
    if (i == 2) {
      data[, i] <- sapply(cleaned_col, function(x) paste(x, collapse = ", "))
    } else {
      data[, i] <- sapply(cleaned_col, function(x) paste(unique(x), collapse = ", "))
    }
  }
  return(data)
}

# function to clean the protein IDs in the Input data
ProteinIDCleaningFunction <- function(IDs) {
  
  # step one, deal with NX_ etries
  # Remove "NX_" prefix and any the suffix following a dash
  cleaned_IDs <- gsub("(^NX_[^-]+)-.*", "\\1", IDs)
  
  
  
  # step 2, deal with sp/ID/ cases
  # Keep everything after "sp|" and before the second "|"
  cleaned_IDs <- gsub(".*sp\\|([^|]+)\\|.*", "\\1", cleaned_IDs)
  
  
  # step 3 deal with gi|315259111|ref|NP_001186752.1| cases 
  
  cleaned_IDs <- gsub(".*gi\\|[^|]+\\|ref\\|(NP_[^|]+).*", "\\1", cleaned_IDs)
  # looks for strings starting with gi and keeps 
  # all text after the first occurance of NP_ up until but excluding the first subsequent |
  
  
  # # Step 4: Handle cases starting with "CONTRIB" and keep everything after the last underscore "_"
  # cleaned_IDs <- gsub("^CONTRIB_([^_]+)_.*", "\\1", cleaned_IDs)
  # ### NB: These are gene names and not UniProt Entry IDs, so we actually want to do that step later and convert from gene name to UniProt ID?
  # after cleaning,keep only the unique IDs at this stage
  cleaned_IDs <- unique(cleaned_IDs)
  
  return(cleaned_IDs)
}

## for getting next prot and Ref seq IDs from the SwissProt library
split_ids <- function(string) {
  # First, try splitting with "; "
  split_result <- strsplit(string, "; ")
  # If the length of any split_result is 1, it means there was no "; " to split on
  # In this case, try splitting with just ";"
  if (any(sapply(split_result, length) == 1)) {
    split_result <- strsplit(string, ";")
  }
  return(split_result)
}

## to process a sublist of a list of lists where each sublist is a row cell of 
# the protein IDs.

# for each sublist of IDs:
#
# 1) check if any of the IDs have a hit in Swiss_Prot_Human$Entry
# if yes replace that sublist with only the IDs with a hit (unlikely to be multiple)
# 2) if not check if any of the IDs have a hit in Swiss_Prot_Human$neXtProt
# if yes retrieve the Swiss_Prot_Human$Entry for that row of Siwss_Prot_Human, then retain that Entry ID and replace the sublist with it
# 3) if not, check if any of the IDs have a hit in Swiss_Prot_Human$RefSeq
# if yes retrieve the Swiss_Prot_Human$Entry for that row of Siwss_Prot_Human, then retain that Entry ID and replace the sublist with it
# if not keep the sublist as it is


process_sublist <- function(sublist, Swiss_Prot_Human, neXtProt_list, RefSeq_list) {
  # Check for hits in Swiss_Prot_Human$Entry
  entry_hits <- sublist[sublist %in% Swiss_Prot_Human$Entry]
  if (length(entry_hits) > 0) {
    return(entry_hits)
  }
  
  
  # Check for hits in Swiss_Prot_Human$neXtProt
  for (id in sublist) {
    if (any(sapply(neXtProt_list, function(x) id %in% x))) {
      matching_row <- which(sapply(neXtProt_list, function(x) id %in% x))
      return(Swiss_Prot_Human$Entry[matching_row])
    }
  }
  
  # Check for hits in Swiss_Prot_Human$RefSeq
  for (id in sublist) {
    if (any(sapply(RefSeq_list, function(x) id %in% x))) {
      matching_row <- which(sapply(RefSeq_list, function(x) id %in% x))
      return(Swiss_Prot_Human$Entry[matching_row])
    }
  }
  
  # If no matches, return the original sublist
  return(sublist)
}


#### subsetting swiss prot gene ids amnd returning uniprot ones for enrichment
extractEnsemblGeneIDs <- function(dataSet, uniprotMapping) {
  # subset to only include IDs present in UniProt mapping
  subsetData <- dataSet[dataSet$cleaned_protein_IDs %in% uniprotMapping$uniprotswissprot,]
  
  # match SwissProt IDs in the UniProt data frame to get the position index
  matchedIndices <- match(subsetData$cleaned_protein_IDs, uniprotMapping$uniprotswissprot)
  
  # rxtract corresponding ENSEMBL Gene IDs
  ensemblGeneIDs <- uniprotMapping$ensembl_gene_id[matchedIndices]
  
  # remove NAs if any (in case some SwissProt IDs didn't have a match)
  ensemblGeneIDs <- na.omit(ensemblGeneIDs)
  
  return(ensemblGeneIDs)
}