library(ggplot2)
library(viridis)
# now we want to apply a function where: 

plot_colour_coded_den <- function(absolute_file_path, 
                                        peptidoform_metadata,
                                        column_to_colour_by,
                                        mz_binwidth = 0.02) {
  
  # first, we read in the .csv file for a peptidoform
  temp_data <- read.csv(absolute_file_path)
  
  
  # calculate total density for all data points across datasets
  total_density_data <- temp_data %>% 
    ggplot(aes(x = calibrated_error)) + 
    geom_density(fill = "grey60", alpha = 0.5)
  
  density_plots <- lapply(unique(temp_data$dataset_ID), function(dataset_id) {
    dataset_data <- temp_data[temp_data$dataset_ID == dataset_id, ]
    ggplot(dataset_data, aes(x = calibrated_error)) + 
      geom_density(aes(fill = dataset_ID), alpha = 0.5) + 
      # Apply the same axis limits and breaks as in total_density_data for consistency
      labs(title = dataset_id, x = "Calibrated Error (m/z)", y = "Density") +
      theme_minimal() +
      theme(legend.position = "none", # Hide individual legends for each plot
            plot.title = element_text(hjust = 0.5)) # Center the plot titles
  })
  
  
  # Combine the total density and individual densities into one plot
  combined_plot <- ggplot() + 
    geom_density(data = temp_data, aes(x = calibrated_error), fill = "grey80", alpha = 0.5) + 
    facet_wrap(~ dataset_ID, scales = "free", nrow = 2) +
    theme_minimal() +
    labs(x = "Calibrated Error (m/z)", y = "Density")
  
  for (plot in density_plots) {
    combined_plot <- combined_plot + layer(data = plot$data, mapping = plot$mapping, stat = plot$stat, position = plot$position, geom = "density", fill = plot$fill, alpha = plot$alpha)
  }
  
  # Print the combined plot
  print(combined_plot)
  
}
 
