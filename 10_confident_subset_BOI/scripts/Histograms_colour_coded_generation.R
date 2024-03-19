library(tidyverse)
library(ggplot2)
source("plot_colour_coded_histogram_fun.R")

#set wd to script folder
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/10_confident_subset_BOI/scripts/")
# data by peptidoform is in the data dir
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/out/csv_tables_by_peptidoform/"

# peptidoforms of interest metadata - subset of proteins  that could be sulfated based on CC
peptidoform_data <- read.csv("../in/potential_sY_peptidoforms_BOI.csv")


extension = ".csv"
files_by_peptidoform <- list.files(data_dir, pattern = extension, full.names = FALSE, recursive = TRUE)

# retain peptidofrom id from the filenames
ids_in_file_names <- gsub("^merged_data_|\\.csv$", "", files_by_peptidoform)

# keep only relevant files and then generate absolute paths for input
filtered_files <- files_by_peptidoform[ids_in_file_names %in% peptidoform_data$peptidoform_id]
input_filenames <- paste0(data_dir, filtered_files)


plot_colour_coded_histogram(absolute_file_path = input_filenames[37], 
                            peptidoform_metadata = peptidoform_data,
                            mz_binwidth = 0.002)
