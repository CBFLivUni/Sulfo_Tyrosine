library(tidyverse)
library(viridis)
library(ggplot2)


# set wd to script folder
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/10_confident_subset_BOI/scripts/")
# data by peptidoform is in the data dir
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/out/csv_tables_by_peptidoform/"


####### read in data ##############
# peptidoforms of interest metadata - subset of proteins  that could be sulfated based on CC
# peptidoform_data <- read.csv("../in/potential_sY_peptidoforms_BOI.csv")
# stricter subset where histograms were convincing
peptidoform_data <- read.csv("../in/strict_potential_sY_peptidoforms_BOI.csv")


# below will be needed  to add instrument to datasets
datasets_metadata <- read.csv("../in/human_phosphobuild_metadata.csv") %>% distinct()


# convert the issue with multiple datasets split by comma into a long format by adding more rows
# head(datasets_metadata)
datasets_metadata <- datasets_metadata %>% 
  separate_rows(Dataset, sep = ",") %>% 
  distinct()
# rename columns
names(datasets_metadata) <- c("dataset_ID", "experiment_tag", "instrument", "sample_category")

# replace spaces in instrument names by _ to enable proper plotting later - some error in the gfacet wrap otherwise
datasets_metadata$instrument <- gsub(" ", "_", datasets_metadata$instrument)



# 
# # these are for when we have NA 
# manual_instrument_annotations <- read.csv("../in/manually_annotated_exp_tags.csv")

######### get list of files ##############
extension = ".csv"
files_by_peptidoform <- list.files(data_dir, pattern = extension, full.names = FALSE, recursive = TRUE)

# retain peptidofrom id from the filenames
ids_in_file_names <- gsub("^merged_data_|\\.csv$", "", files_by_peptidoform)

# keep only relevant files and then generate absolute paths for input
filtered_files <- files_by_peptidoform[ids_in_file_names %in% peptidoform_data$peptidoform_id]
input_filenames <- paste0(data_dir, filtered_files)

# decide which columns to colour plots by
# columns_of_interest <- c("dataset_ID", "experiment_tag",  "instrument")
columns_of_interest <- c("instrument")

# add source here so we can generate the colour palette within the function

source("plot_colour_coded_histogram_fun.R")

for (column in columns_of_interest) {
  
  # open a pdf file#
  pdf(file = paste0("../out/faceted_plots_by_", column, ".pdf"), 
      width = 8.27, height = 11.69) #A4 format
  
  for (file in input_filenames) {
    
    # add plot to pdf
    plot_colour_coded_histogram(absolute_file_path = file, 
                                peptidoform_metadata = peptidoform_data,
                                column_to_colour_by = column,
                                mz_binwidth = 0.001)
    
    
  }
  
  dev.off()
  
  
  # open a pdf file#
  pdf(file = paste0("../out/stacked_plots_by_", column, ".pdf"), 
      width = 8.27, height = 11.69) #A4 format
  
  for (file in input_filenames) {
    
    # add plot to pdf
    plot_colour_coded_histogram(absolute_file_path = file, 
                                peptidoform_metadata = peptidoform_data,
                                mz_binwidth = 0.001,
                                column_to_colour_by = column,
                                layout = "stacked_bar")
    
    
  }
  
  dev.off()
  
}









