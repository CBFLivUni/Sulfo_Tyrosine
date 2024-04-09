library(purrr)
library(dplyr)
library(stringr)



modify_USI <- function(usi) {
  # split the USI string to work with the last segment
  parts <- str_split(usi, ":")[[1]]
  last_segment <- parts[6]
  modified <- FALSE
  
  # replace a random occurance of "Y[Phospho]" with "Y{Sulfo}" if Y[Phospho] exists
  if (str_detect(last_segment, "Y\\[Phospho\\]")) {
   
    
     phospho_positions <- gregexpr("Y\\[Phospho\\]", last_segment)[[1]]
    if (length(phospho_positions) > 0) {
      # choose a random position to replace
      chosen_position <- sample(phospho_positions, 1)
      last_segment <- sub("Y\\[Phospho\\]", "Y{Sulfo}", last_segment, perl = TRUE)
      # reconstruct usi and return
      parts[6] <- last_segment
      return(paste(parts, collapse = ":"))
    }
  }
  
  if (!modified && str_detect(last_segment, "Y(?!\\{)") && str_detect(last_segment, "\\[Phospho\\]")) {
    # Find all "Y" not followed by "{"
    y_positions <- gregexpr("Y(?!\\{)", last_segment, perl = TRUE)
    y_positions <- y_positions[[1]][1:length(y_positions[[1]])]
    
    if (length(y_positions) == 0) {
     return("") 
    } else if (length(y_positions) == 1) {
      
      chosen_y <- y_positions[1]
    } else {
      # choose random y to replace
      chosen_y <- sample(y_positions, 1)
    }
     
      
      # split the sting into a vector fo letters
      last_segment <- strsplit(last_segment, split = "") %>% unlist()
      
      sulfotyrosine <- "Y{Sulfo}"
      
      # replace the chosen tyrosine with a sulfotyrosine
      last_segment[chosen_y] <-  sulfotyrosine 
      
      # re-paste together
      new_segment <- paste(last_segment, collapse = "")
      
    
    
    
    # now that we have added a sulfo, we need to remove a phospho
    
    # Find all "[Phospho]" in this new segment
    phospho_positions <- gregexpr("\\[Phospho\\]", new_segment)[[1]][1:length(gregexpr("\\[Phospho\\]", new_segment)[[1]])]
    
    if (length(phospho_positions) == 0) {
      return("") 
    } else if (length(phospho_positions) == 1) {
      
      chosen_phospho <- phospho_positions[1]
    } else {
      # choose random phospho to replace
      chosen_phospho <- sample(phospho_positions, 1)
    }
    
    
    # ger the new segment as vector again
    temp_segment <- strsplit(new_segment, split = "") %>% unlist()
    # replace the chosen phospho with nothing
    temp_segment[chosen_phospho:(chosen_phospho+8)] <- ""
   
    # reconstruct
    new_segment <- paste(temp_segment, collapse = "")
    
      
    modified <- TRUE
    
  }
  
  
  # if no modifications were made, return an empty string
  if (!modified) {
    return("")
  }
  
  # otherwise, reassemble the USI with the modified last segment
  parts[6] <- new_segment
  return(paste(parts, collapse = ":"))
}


   usi
   usi1 <- modify_USI (usi)
   usi1
   usi2 <- modify_USI (usi1)
   usi2
   usi3 <- modify_USI (usi2)
   usi3
# usi<- "mzspec:PXD999957:08CPTAC_GBM_P_PNNL_20190306_B2S4_f02:scan:47380:[TMT6plex]-ATWLSLFSS[Phospho]EES[Phospho]NLGANNYDDYR[TMT6plex]/3"
# modify_USI <- function(usi) {
#   # check if "Y[Phospho]" is present
#   if (str_detect(usi, "Y\\[Phospho\\]")) {
#     # replace "Y[Phospho]" with "Y{Sulfo}"
#     return(str_replace(usi, "Y\\[Phospho\\]", "Y{Sulfo}"))
#   } else {
#     # if "Y[Phospho]" is not present
#     # randomly replace one "Y" with "Y{Sulfo}" 
#     all_Y_positions <- unlist(gregexpr("Y", usi))
#     selected_Y_position <- sample(all_Y_positions, 1)
#     
#     substr_replace(usi, "Y{Sulfo}", selected_Y_position, selected_Y_position)
#     
#     
#     new_usi <- substr_replace(string = usi, replacement = "Y{Sulfo}", 
#                               start = selected_Y_position, end = selected_Y_position)
#     
#     # randomly remove one "[Phospho]" if present
#     all_phospho_positions <- gregexpr("\\[Phospho\\]", usi)[[1]]
#     selected_phospho_position <- sample(all_phospho_positions, 1)
#     usi <- substr_replace(usi, "", selected_phospho_position, selected_phospho_position + nchar("[Phospho]") - 1)
#     
#     
#     return(usi)
#   }
# }




generate_USI <- function(df) {
  # Pre-allocate the USI column to speed up the process
  df$USI <- NA
  
  for (i in 1:nrow(df)) {
    # Extract necessary components for each row
    current_dataset <- sub("([^-]*)(-.*)?", "\\1", df$dataset_ID[i])
    msRun <- gsub("^(.*?)\\..*$", "\\1", df$spectrum[i])
    scan_number <- gsub("^[^.]*\\.([^.]+)\\..*$", "\\1", df$spectrum[i])
    peptidoform <- df$mod_peptide[i] %>%
      # substitute "N[115]" with "{deamidation}"
      gsub("N\\[115\\]", "N[Deamidated]", .) %>%
      # substitute "M[147]" with "{oxidation}"
      gsub("M\\[147\\]", "M[Oxidation]", .) %>%
      # substitute "Y[243]", "T[181]", and "S[167]" with "[Phospho]"
      gsub("Y\\[243\\]", "Y[Phospho]", .) %>%
      gsub("T\\[181\\]", "T[Phospho]", .) %>%
      gsub("S\\[167\\]", "S[Phospho]", .) %>%
      # substitute "n[230]" and "n[145]" with "[TMT6plex]-" at the start and also add "[TMT6plex]" at the end of the string
      gsub("^n\\[230\\]", "[TMT6plex]-", .) %>%
      gsub("^n\\[145\\]", "[TMT6plex]-", .)
    
    #  if peptidoform starts with [TMT6plex]-, add [TMT6plex] to the end
    if(grepl("^\\[TMT6plex\\]-", peptidoform)) {
      peptidoform <- paste0(peptidoform, "[TMT6plex]")
    }
    # complete spectrum interpretation segment
    spectrum_interpretation <- paste(peptidoform, "/", df$z[i], sep = "")
    
    # Generate USI
    df$USI[i] <- paste("mzspec", current_dataset, msRun, "scan", scan_number, spectrum_interpretation, sep = ":")
  }
  
  # Create a new column with the modified USI values
  df$USI_workaround <- gsub("CPTAC_S0", "PXD9999", df$USI)
  
  # apply the modify_USI function to each USI in the dataframe
  df$USI_with_one_sulfo <- sapply(df$USI_workaround, modify_USI)
  df$USI_with_two_sulfo <- sapply(df$USI_with_one_sulfo, modify_USI)
  
  filename <- unique(df$peptidoform_id)
  write.csv(file = paste0("../out/added_USIs_to_", filename, ".csv"), x = df)
  return(df)
}
