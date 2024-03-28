library(purrr)
library(dplyr)


generate_USI <- function(df) {
  # Pre-allocate the USI column to speed up the process
  df$USI <- NA
  
  for (i in 1:nrow(df)) {
    # Extract necessary components for each row
    current_dataset <- sub("([^-]*)(-.*)?", "\\1", df$dataset_ID[i])
    msRun <- gsub("^(.*?)\\..*$", "\\1", df$spectrum[i])
    scan_number <- gsub("^[^.]*\\.([^.]+)\\..*$", "\\1", df$spectrum[i])
    spectrum_interpretation <- paste(df$peptide[i], "/", df$z[i], sep = "")
    
    # Generate USI
    df$USI[i] <- paste("mzspec", current_dataset, msRun, "scan", scan_number, spectrum_interpretation, sep = ":")
  }
  
  # Create a new column with the modified USI values
  df$USI_workaround <- gsub("CPTAC_S0", "PXD9999", df$USI)
  filename <- unique(df$peptidoform_id)
  write.csv(file = paste0("../out/added_USIs_to", filename, ".csv"), x = df)
  return(df)
}
