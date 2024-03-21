library(ggplot2)
library(viridis)
library(forcats)
library(dplyr)





# create colour palette for all instruments
instruments <- unique(datasets_metadata$instrument)

# number of unique instruments
num_instruments <- length(instruments)

# create a colorblind-friendly palette
color_palette <- viridis::viridis(n = num_instruments, option = "D")
# color_palette <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
#                    "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70",
#                    "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")[1:num_instruments]
# 
# # assign colors to each instrument
instrument_colors <- setNames(color_palette, instruments)


# pick palette by actually visualising - in future may be worth
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, length(instruments)), axes=FALSE)

# add a rectangle for each color, along with the instrument name as text
for(i in 1:length(instruments)) {
  rect(0, i-1, 1, i, col=instrument_colors[i], border="white")
  text(0.5, i-0.5, instruments[i], cex=0.7)
}



plot_colour_coded_histogram <- function(absolute_file_path, 
                                        peptidoform_metadata,
                                        column_to_colour_by,
                                        mz_binwidth = 0.02,
                                        layout = "facet",
                                        metadata_by_dataset = datasets_metadata ) {# set as default as i wont rly customise
  
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
  
  # replace NA values by "" to avoid future instrument assignment complications
  
  temp_data <- temp_data %>%
    mutate(
      dataset_ID = replace_na(dataset_ID, ""),
      experiment_tag = replace_na(experiment_tag, "")
    )
  #  add instrument metadata to the dataset. 
  temp_data$instrument <- ""
  unknowns <- data.frame()
  for(i in 1:nrow(temp_data)) {
    current_dataset_ID <- gsub("^CPTAC_", "", temp_data$dataset_ID[i]) # get rid of CPTAC_ prefix
    
    
    current_experiment_tag <- temp_data$experiment_tag[i]
    
    match_row <- metadata_by_dataset[
      metadata_by_dataset$dataset_ID == current_dataset_ID & 
        metadata_by_dataset$experiment_tag == current_experiment_tag, 
    ]
    
    if(nrow(match_row) == 1) {
      temp_data$instrument[i] <- match_row$instrument
    } else {
      # warning(paste("Dataset ID", current_dataset_ID, "+ experiment tag not matched, attempting instrument assignment by dataset ID only."))
      
      match_row <- metadata_by_dataset[metadata_by_dataset$dataset_ID == current_dataset_ID, ]
      
      if (nrow(match_row) == 1) {
        temp_data$instrument[i] <- match_row$instrument
      } else if (nrow(match_row) > 1) {
        instrument <- unique(match_row$instrument)
        
        if (length(instrument) == 1) {
          temp_data$instrument[i] <- instrument
        } else {
          # print(paste("Current experiment tag:", current_experiment_tag))
          # print(paste("Current dataset ID:", current_dataset_ID))
          # warning("More than one possible instrument, the most likely was picked, please double-check.")
          
          instr_counts <- table(match_row$instrument)
          most_common_instrument <- names(sort(instr_counts, decreasing = TRUE))[1]
          temp_data$instrument[i] <- most_common_instrument
        }
        
      } else {
        temp_data$instrument[i] <- "not_known"
        warning(paste("No instrument found for dataset ID", current_dataset_ID, "."))
      }
    }
  }
  # 
  # for(i in 1:nrow(temp_data)) {
  #   # extract the current row's dataset_ID and experiment_tag
  #   current_dataset_ID <- temp_data$dataset_ID[i]
  #   # print(current_dataset_ID)
  #   current_experiment_tag <- temp_data$experiment_tag[i]
  #   # print(current_experiment_tag)
  #   
  #   # find a matching row in metadata_by_dataset
  #   match_row <- metadata_by_dataset[metadata_by_dataset$dataset_ID == current_dataset_ID & 
  #                                    metadata_by_dataset$experiment_tag == current_experiment_tag, ]
  #   # print(match_row)
  #   # if a match is found, update the instrument in temp_data
  #   if(nrow(match_row) == 1) {
  #     temp_data$instrument[i] <- match_row$instrument
  #   } else { # usually the problem when no match found is the experimental tag
  #     warning("datasetd ID + experimetn tag not matched, attempting instrument assignment by dataset ID only")
  #     # get the row by just matching the dataset ID 
  #     match_row <- metadata_by_dataset[metadata_by_dataset$dataset_ID == current_dataset_ID ,]
  #     
  #     if (nrow(match_row) == 1) {
  #       
  #       temp_data$instrument[i] <- match_row$instrument
  #       
  #     } else if (nrow(match_row) > 1) { # here we could have multiple matching rows! 
  #       # usually same instrument in a dataset so try to narrow it down to just one
  #       instrument <- unique(match_row$instrument)
  #     
  #       if (length(instrument) == 1) {
  #         
  #         temp_data$instrument[i] <- instrument
  #         
  #       } else { # if more than one leave NA and stop function 
  #         print(current_experiment_tag)
  #         print(current_dataset_ID)
  #         warning("more than one possible instrument, the most likely was picked, please double-check")
  #         # count how many of each instrument in match_row
  #         instr_counts <- table(match_row$instrument)
  #         # pick the one that appears most often
  #         most_common_instrument <- names(sort(instr_counts, decreasing = TRUE))[1]
  #         temp_data$instrument[i] <- most_common_instrument
  #       }
  #       
  #     } else {
  #       # No matches found, instrument unknown
  #       temp_data$instrument[i] <- "unknown"
  #       warning("No instrument found for dataset ID ", current_dataset_ID, ".")
  #     
  #   }
  #   
  # }
  

  
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
  
  if (!column_to_colour_by %in% names(temp_data)) {
    stop("The specified column to color by does not exist in the data frame.")
  }
  
  # aggregate data for plotting
  # here we need to aggregate by ~bin + column to colourr by;
  # old version:
  # plot_data <- aggregate(cbind(count = calibrated_error) ~ bin + dataset_ID, data = temp_data, FUN = length)

  # i guess a work around is this:
  formula_text <- as.formula(paste("cbind(count = calibrated_error) ~ bin +", column_to_colour_by))
  plot_data <- aggregate(formula_text, data = temp_data, FUN = length)
  
  
  # get bin center for plotting instead of using bin number
  plot_data$binstart <- bins[plot_data$bin]
  plot_data$binend <- bins[plot_data$bin+1]
  plot_data$bincenter <- (plot_data$binstart + plot_data$binend) / 2
  
  # calculate total PSM counts per dataset
  total_counts <- plot_data %>%
    group_by(!!sym(column_to_colour_by)) %>% # group by column to colour by 
    summarise(total = sum(count))
  
  # convert to factor
  plot_data[[column_to_colour_by]] <- factor(plot_data[[column_to_colour_by]], levels = total_counts[[column_to_colour_by]])
  
  # add to plot df and custom facet label 
  plot_data <- plot_data %>%
    left_join(total_counts, by = column_to_colour_by) %>%
    mutate(facet_label = paste0(get(column_to_colour_by), "_total_PSMs:", total)) # do for column to colour by - here we need to sue get
  
  head(plot_data)
  
  title_text <- gsub("^.*/|\\.csv$|^", "", absolute_file_path) # remove path and .csv
  title_text <- sub("^merged_data_", "", title_text) # remove merged_data_ prefix

  if (layout == "facet") {
    
    # generate and display the plot with facets for each dataset_ID
    # order to have dataset with most at the top
    p <- ggplot(plot_data, aes(x = bincenter, y = count, fill = get(column_to_colour_by))) +
      geom_col(width = mz_binwidth) +
      facet_wrap(~reorder(facet_label, -total), ncol = 3, labeller = label_parsed) +  
      scale_x_continuous(name = "Calibrated Error (m/z)") +
      labs(y = "Count", title = title_text) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.x = element_text(size = 10, hjust = 1)) +
      geom_vline(xintercept = -0.0095, linetype = "dashed", color = "black") 
 
  } 
  
  if (layout == "stacked_bar") {
    
    
    # create a stacked bar chart type of histogram
    p <- ggplot(plot_data, aes(x = bincenter, y = count, fill = get(column_to_colour_by))) +
      geom_col(position = "stack", width = mz_binwidth) +  # geom_col is used here for explicit bar widths
      scale_x_continuous(name = "Calibrated Error (m/z)", breaks = seq(min(plot_data$binstart), max(plot_data$binend), by = mz_binwidth * 5)) +
      labs(y = "PSM Count", title = title_text) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
      geom_vline(xintercept = -0.0095, linetype = "dashed", color = "black")
    
    # adjust axis limits
    p <- p + expand_limits(x = c(min(plot_data$binstart), max(plot_data$binend)))
    
  }
 
  if (column_to_colour_by == "instrument") {
    
    p <- p + scale_fill_manual(values = instrument_colors)
  } else {
    
    p <- p + scale_fill_brewer(palette = "Set3") 
  }
  
  print(p)
  return(p)
}
 
