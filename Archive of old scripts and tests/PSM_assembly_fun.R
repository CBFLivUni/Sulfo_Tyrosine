library(dplyr)
library(purrr)

Assemble_PSMs <- function(df_list) {
  #########
  # this function requires as an input a list of data frames generated using 
  # summarise_by_peptidoform()
  #########
  
  # Merge all dataframes in the list
  merged_full <- reduce(df_list, full_join, by = "peptidoform_id")
  
  # Function to combine and flatten the list of calibrated mass shifts
  flatten_lists <- function(x, y) {
    # Handle NULL values in either x or y
    if (is.null(x)) {
      x <- c()
    }
    if (is.null(y)) {
      y <- c()
    }
    c(unlist(x), unlist(y))
  }
  
  # Aggregate the results
  merged_full_agg <- merged_full %>%
    group_by(peptidoform_id) %>%
    summarise(calibrated_mass_shifts = list(flatten_lists(calibrated_error_all_values.x,
                                                          calibrated_error_all_values.y)))
  
  return(merged_full_agg)
}