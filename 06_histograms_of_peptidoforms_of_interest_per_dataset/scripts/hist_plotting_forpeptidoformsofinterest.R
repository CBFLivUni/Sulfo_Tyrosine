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


# let's get a test subset of just 2 rows with few datasets to run the code on - rows 5 and 21
# peptidoforms_potentially_sulfated <- store[c(5,21),]
# 
# 
# peptidoforms_potentially_sulfated$dataset_ID
#### read in datasets and plot pdfs ####

# Loop through each row in peptidoforms_potentially_sulfated
for (i in 1:nrow(peptidoforms_potentially_sulfated)) {
  
  # extract peptidoform_id and dataset_ID
  current_peptidoform_id <- peptidoforms_potentially_sulfated$peptidoform_id[i]
  print(current_peptidoform_id)
  # the dataset ID column contains all dataset IDs separated by ,
  dataset_ids <- strsplit(peptidoforms_potentially_sulfated$dataset_ID[i], ", ")[[1]]
  
  # create a PDF for the current peptidoform id
  # pdf(paste0("histograms_", current_peptidoform_id, ".pdf"))
  
  # make an empty dataframe for storing merged data for the current peptidoform_id
  merged_data_for_id <- data.frame()
  
  # plot the histograMs of each peptidoform ID within each dataset
  for (dataset_id in dataset_ids) {
    #to test loop
    # dataset_id <- "PXD015901-s015_breast_cancer_phospho_itraq4"
    print(dataset_id)
    
    # get file name and read the dataset
    file_name <- paste0(data_dir, dataset_id, "_peptidoforms.csv")
    dataset <- read.csv(file_name)
    
    # clean peptidoform_id in the current dataset
    dataset$peptidoform_id <- gsub("\\[|\\]", "", dataset$peptidoform_id)
    # head( dataset$peptidoform_id) - has worked
    
    # Filter rows with the current peptidoform_id
    filtered_data <- subset(dataset, peptidoform_id == current_peptidoform_id)
    print(nrow(filtered_data))
    # Merge filtered data for the current peptidoform_id
    # there should always be > 0 because we know this peptidoform is in the dataset; removed if
    # if(nrow(filtered_data) > 0) {
    merged_data_for_id <- rbind(merged_data_for_id, filtered_data)
    # }
    
   #  # add dataset_id and number of rows in filtered_data on the PDF page too
   #  grid.text(paste("Dataset ID:", dataset_id, "\nNumber of rows in filtered data:", nrow(filtered_data)),
   #            x = 0.1, y = 0.95, just = "left", gp = gpar(fontsize = 10))
   #  
   # #  print to page that fewer than 5 PSMs are present in the dataset and
   # # print the filtered data to the page
   #  if (nrow(filtered_data) < 6) {
   #    text_grob <- grid::textGrob(paste("Fewer than 6 PSMs are present in the dataset for \n",
   #                                      current_peptidoform_id, "\nin dataset", dataset_id))
   #    table_grob <- gridExtra::tableGrob(filtered_data$calibrated_error)
   #    gridExtra::grid.arrange(text_grob, table_grob, ncol = 1)
   #  }
    # Generate and plot histogram
    
    # originally thought I'd only plot the histigrams when a certain number o
    # values is present, but it's better to have all the plots; we do need at least 2 points for the density curve though
    
    # if (nrow(filtered_data) > 1) {
    #   mean_calibrated_error <- mean(filtered_data$calibrated_error, na.rm = TRUE)
    #   
    #   # Create the histogram plot
    #   p <- ggplot(filtered_data, aes(x = calibrated_error)) +
    #     geom_histogram(bins = 30, fill = "green", color = "darkgreen", aes(y = ..count..)) + # Histogram bars with counts
    #     geom_density(aes(y = ..count..), color = "#0072B2", size = 1, adjust = 1/3) + # Overlay density line scaled to counts
    #     geom_vline(aes(xintercept = mean_calibrated_error),
    #                color = "#0072B2", linetype = "dashed", size = 1) +
    #     geom_text(aes(x = mean_calibrated_error, y = Inf), 
    #               label = paste("Mean:", round(mean_calibrated_error, 4)), 
    #               color = "#0072B2", vjust = -1.5, hjust = 1.1, size = 3.5) +
    #     ggtitle(paste("Histogram for", current_peptidoform_id, "\nin dataset", dataset_id)) +
    #     xlab("Calibrated Error") +
    #     ylab("Count")
    #   
    #   # Draw the ggplot on the current page
    #   print(p)  # this is needed to actually draw the plot in the PDF
    # }
    
    
    # 
    # if (nrow(filtered_data) > 1) {
    #   mean_calibrated_error <- mean(filtered_data$calibrated_error, na.rm = TRUE)
    #   distvals <- unlist(density(filtered_data$calibrated_error)[2])
    #   ylim <- max(distvals)*5/3
    #   p <- ggplot(filtered_data, aes(x = calibrated_error)) +
    #     geom_histogram(aes(y = ..density..), bins = 30, fill = "green", color = "darkgreen") + # Histogram bars
    #     geom_density(color = "#0072B2", size = 1) +  # density line
    #     geom_vline(aes(xintercept = mean_calibrated_error),
    #                color = "#0072B2", linetype = "dashed", size = 1) +  
    #     geom_text(aes(x = mean_calibrated_error, y = ylim), 
    #                   label = round(mean_calibrated_error, 4), 
    #               color = "#0072B2", 
    #               vjust = -0.5, 
    #               hjust = -0.1, 
    #               size = 3.5) +
    #     ggtitle(paste("Histogram for", current_peptidoform_id, "\nin dataset", dataset_id)) +
    #     xlab("Calibrated Error") +
    #     ylab("Density")
    #   
    #   # draw the ggplot on the current page
    #   print(p)  # this is needed to actually draw the plot in the PDF
    # }
      # add a page break in the PDF for each histogram
      # grid::grid.newpage()
      
    # } else {
    #   
    #   # print to page that fewer than 5 PSMs are present in the dataset and 
    #   # print the filtered data to the page
    #   text_grob <- grid::textGrob(paste("Fewer than 5 PSMs are present in the dataset for \n", 
    #                                     current_peptidoform_id, "\nin dataset", dataset_id))
    #   table_grob <- gridExtra::tableGrob(filtered_data$calibrated_error)
    #   gridExtra::grid.arrange(text_grob, table_grob, ncol = 1)
    #   # add a page break after the message and table
    #   grid::grid.newpage()
    # }
      
      
   
  }
  
  # # Close the PDF device
  # dev.off()
  # Write the merged data for the current peptidoform_id to a CSV file
  write.csv(merged_data_for_id, paste0(project_dir, "out/merged_data_", current_peptidoform_id, ".csv"), row.names = FALSE)
}

