library(ggplot2)
library(viridis)
# now we want to apply a function where: 

plot_colour_coded_histogram <- function(absolute_file_path, 
                                        peptidoform_metadata,
                                        column_to_colour_by,
                                        mz_binwidth = 0.02) {
  
  # first, we read in the .csv file for a peptidoform
  temp_data <- read.csv(absolute_file_path)
  
  # based on temp_data$dataset_ID we need to bin temp_data$calibrated_error into 0.02 m/z wide bins. 
  # find the min calibrated error, round down to nearest 0.1; and max calibrated error, round up top enarest 0.1
  # this will be our x axis lims
  # split the callibrated errors into bins of mz_binwidth m/z - default for mz_binwidth is 0.02
  # somehow keep track of the number of values for rach dataset ID that fall in a bin, then create a "hiostogram" in the form
  # of a stacked bar chart - for each m/z bin the bar should be split and coloure dby datasrt ID.
  
  full_IDs <- strsplit(temp_data$dataset_ID, split = "-")
  
  # clean dataset ID and keep separate from experiment tag
  temp_data$dataset_ID <- sapply(full_IDs, `[`, 1) # extracts the first element
  temp_data$experiment_tag <- sapply(full_IDs, `[`, 2) # extracts the second element
  
  # precalculate min and max calibrated error
  min_error <- round(min(temp_data$calibrated_error) - 2*mz_binwidth, digits = 2)
  max_error <- round(max(temp_data$calibrated_error) + 2*mz_binwidth, digits = 2)
  

  
  # create bins to plot data by
  bins <- seq(min_error, max_error, by = mz_binwidth)
  
  
  # cut the calibrated errors into bins
  temp_data$bin <- cut(temp_data$calibrated_error, bins, include.lowest = TRUE, labels = FALSE)
  
  ## TODO add functionality to specify column to colour by in the metadata. 
  ## especially interested in colouring by dataset ID, experiment tag, instrument, 
  ## and instrument sensitivity.
  # ask Andy to rate instrument sensitivity as low, medium, high from what he knows
  # ask Andy for input on USI generation
  
  # aggregate data for plotting
  plot_data <- aggregate(cbind(count = calibrated_error) ~ bin + dataset_ID, data = temp_data, FUN = length)
  
  # get bin center for plotting instead of using bin number
  plot_data$binstart <- bins[plot_data$bin]
  plot_data$binend <- bins[plot_data$bin+1]
  plot_data$bincenter <- (plot_data$binstart + plot_data$binend) / 2
  
  head(plot_data)
  # create a stacked bar chart type of histogram
  p <- ggplot(plot_data, aes(x = bincenter, y = count, fill = dataset_ID)) +
    geom_col(position = "stack", width = mz_binwidth) +  # geom_col is used here for explicit bar widths
    scale_x_continuous(name = "Calibrated Error (m/z)", breaks = seq(min(plot_data$binstart), max(plot_data$binend), by = mz_binwidth * 5)) +
    labs(y = "PSM Count", title = gsub("^.*/|\\.csv$", "", absolute_file_path)) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") # Use a distinct color palette
  
  # adjust axis limits 
  p <- p + expand_limits(x = c(min(plot_data$binstart), max(plot_data$binend)))
  print(p)
  return(p)
}
 
