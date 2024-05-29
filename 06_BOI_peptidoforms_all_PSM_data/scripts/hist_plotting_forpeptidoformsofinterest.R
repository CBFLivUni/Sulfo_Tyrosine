library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)


##### directories ####
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/scripts/")


project_dir = "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/"

# we want to use the peptidoform ID aggregated data from individual experiments
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/01_Sulfotyrosine_Data_Processing/out/new_peptidoform_IDs/"
gc()


##### clean file names #####
# specify the extension unique to your calibrated files
extension = "_peptidoforms.csv"

# get all files of interest in absolut path
data_files <- list.files(data_dir, pattern = extension, full.names = TRUE, recursive = TRUE)

# and get file names
input_filenames <- basename(data_files)

clean_names <- gsub(extension, "", input_filenames)


#### read in peptidoforms #####
peptidoforms_potentially_sulfated <- read.csv(paste0(project_dir, "in/potential_sulfo_peptidoforms_.csv"))

# we have added a few new ones, lets update the list; we only need to re-run 
# the script for these that were not in the original list
updated_potentially_sulfated <- read.csv(paste0(project_dir, "in/table_S1_final_v2.csv"))
all_Ycont_BOIs <- read.csv(paste0(project_dir, "in/BOIs_all.csv"))
all_Ycont_BOIs <- all_Ycont_BOIs[all_Ycont_BOIs$peptidoform_id %in% updated_potentially_sulfated$peptidoform_ID,]


peptidoforms_potentially_sulfated <- all_Ycont_BOIs[!all_Ycont_BOIs$peptidoform_id %in% peptidoforms_potentially_sulfated$peptidoform_id,]

# let's get a test subset of just 2 rows with few datasets to run the code on - rows 5 and 21
# peptidoforms_potentially_sulfated <- store[c(5,21),]
# 
# 
# peptidoforms_potentially_sulfated$dataset_ID
#### read in datasets and plot pdfs ####

# Loop through each row in peptidoforms_potentially_sulfated
for (i in 1:nrow(peptidoforms_potentially_sulfated)) {
  
  # Extract peptidoform_id and dataset_ID
  current_peptidoform_id <- peptidoforms_potentially_sulfated$peptidoform_id[i]
  print(current_peptidoform_id)
  
  # The dataset ID column contains all dataset IDs separated by a comma
  dataset_ids <- strsplit(peptidoforms_potentially_sulfated$dataset_ID[i], ", ")[[1]]
  
  # Create a PDF for the current peptidoform id
  pdf(paste0("histograms_", current_peptidoform_id, ".pdf"))
  
  # Make an empty dataframe for storing merged data for the current peptidoform_id
  merged_data_for_id <- data.frame()
  
  # Plot the histograms of each peptidoform ID within each dataset
  for (dataset_id in dataset_ids) {
    print(dataset_id)
    
    # Get file name and read the dataset
    file_name <- paste0(data_dir, dataset_id, "_peptidoforms.csv")
    dataset <- read.csv(file_name)
    
    # Clean peptidoform_id in the current dataset
    dataset$peptidoform_id <- gsub("\\[|\\]", "", dataset$peptidoform_id)
    
    # Filter rows with the current peptidoform_id
    filtered_data <- subset(dataset, peptidoform_id == current_peptidoform_id)
    print(nrow(filtered_data))
    
    # Merge filtered data for the current peptidoform_id
    merged_data_for_id <- rbind(merged_data_for_id, filtered_data)
    
    # Add dataset_id and number of rows in filtered_data on the PDF page
    grid.text(paste("Dataset ID:", dataset_id, "\nNumber of rows in filtered data:", nrow(filtered_data)),
              x = 0.1, y = 0.95, just = "left", gp = gpar(fontsize = 10))
    
    # Print to page if fewer than 6 PSMs are present in the dataset
    if (nrow(filtered_data) < 6) {
      text_grob <- grid::textGrob(paste("Fewer than 6 PSMs are present in the dataset for \n",
                                        current_peptidoform_id, "\nin dataset", dataset_id))
      table_grob <- gridExtra::tableGrob(filtered_data$calibrated_error)
      gridExtra::grid.arrange(text_grob, table_grob, ncol = 1)
    } else {
      # Generate and plot histogram
      mean_calibrated_error <- mean(filtered_data$calibrated_error, na.rm = TRUE)
      distvals <- unlist(density(filtered_data$calibrated_error)[2])
      ylim <- max(distvals) * 5 / 3
      
      p <- ggplot(filtered_data, aes(x = calibrated_error)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "green", color = "darkgreen") +
        geom_density(color = "#0072B2", size = 1) +
        geom_vline(aes(xintercept = mean_calibrated_error), color = "#0072B2", linetype = "dashed", size = 1) +
        geom_text(aes(x = mean_calibrated_error, y = ylim), label = round(mean_calibrated_error, 4),
                  color = "#0072B2", vjust = -0.5, hjust = -0.1, size = 3.5) +
        ggtitle(paste("Histogram for", current_peptidoform_id, "\nin dataset", dataset_id)) +
        xlab("Calibrated Error") +
        ylab("Density")
      
      # Draw the ggplot on the current page
      print(p)
    }
    
    # Add a page break in the PDF for each histogram
    grid::grid.newpage()
  }
  
  # Write the merged data for the current peptidoform_id to a CSV file
  write.csv(merged_data_for_id, paste0(project_dir, "out/merged_data_", current_peptidoform_id, ".csv"), row.names = FALSE)
  
 
  dev.off()
}
 