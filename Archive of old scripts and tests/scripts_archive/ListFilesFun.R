getfiles <- function(folder, extension, get_prefix = TRUE) {
  # List all files recursively in the folder
  
  folder = wd
  extension = ".pep_calibrated.tsv"
  
  # List all files recursively in the folder with the specified extension
  all_files <- list.files(wd, pattern = extension, full.names = TRUE, recursive = TRUE)
  relative_paths <- 
  
    # Get relative paths
  relative_paths <- sapply(all_files, function(file) {
    normalizePath(file, winslash = "/", mustWork = FALSE)
  })
  
  # Initialize prefixes vector
  prefixes <- character()
  
  if (get_prefix) {
    # Extract prefixes
    prefixes <- sapply(relative_paths, function(path) {
      # Replace backslashes with underscores and remove the extension
      prefix <- gsub("\\\\", "_", path)
      prefix <- sub("_interact-ipro-ptm\\.pep\\.xml$", "", prefix)
      prefix <- sub("_interact-prob$", "", prefix)
      dirname(prefix)
    })
  }
  
  if (get_prefix) {
    return(list(relative_paths, prefixes))
  } else {
    return(relative_paths)
  }
}

# Example usage
# files_and_prefixes <- getfiles("path/to/folder", "xml")
# files_only <- getfiles("path/to/folder", "xml", FALSE)
