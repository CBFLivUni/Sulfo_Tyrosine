library(ggplot2)
library(viridis)
library(forcats)
library(dplyr)
library(RColorBrewer)
library(pals)
library(Polychrome)




# create colour palette for all instruments
instruments <- unique(datasets_metadata$instrument) %>% sort()

# number of unique instruments
num_instruments <- length(instruments)


# create your own color palette (16 colors) based on `seedcolors`
color_palette <- createPalette(num_instruments,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(color_palette)



# # create a colorblind-friendly palette
# color_palette <- RColorBrewer::RColorBrewer(n = num_instruments, option = "D")

barplot(rep(1, num_instruments), col = color_palette, border = "NA", 
        main = "Color Palette Visualization", 
        axes = FALSE)
# color_palette <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
#                    "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70",
#                    "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")[1:num_instruments]
# 
# # assign colors to each instrument
instrument_colors <- setNames(color_palette, instruments)


# pick palette by actually visualising 
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
                                        metadata_by_dataset = datasets_metadata,
                                        totals_only = FALSE) {
  # Read in the .csv file for a peptidoform
  temp_data <- read.csv(absolute_file_path)
  
  # Extract dataset ID and experiment tag
  full_IDs <- strsplit(temp_data$dataset_ID, split = "-")
  temp_data$dataset_ID <- sapply(full_IDs, `[`, 1)
  temp_data$experiment_tag <- sapply(full_IDs, `[`, 2)
  
  # Replace NA values with "" to avoid future complications
  temp_data <- temp_data %>%
    mutate(
      dataset_ID = replace_na(dataset_ID, ""),
      experiment_tag = replace_na(experiment_tag, "")
    )
  
  # Add instrument metadata to the dataset
  temp_data$instrument <- ""
  for(i in 1:nrow(temp_data)) {
    current_dataset_ID <- gsub("^CPTAC_", "", temp_data$dataset_ID[i])
    current_experiment_tag <- temp_data$experiment_tag[i]
    
    match_row <- metadata_by_dataset[
      metadata_by_dataset$dataset_ID == current_dataset_ID & 
        metadata_by_dataset$experiment_tag == current_experiment_tag, 
    ]
    
    if(nrow(match_row) == 1) {
      temp_data$instrument[i] <- match_row$instrument
    } else {
      match_row <- metadata_by_dataset[metadata_by_dataset$dataset_ID == current_dataset_ID, ]
      if (nrow(match_row) == 1) {
        temp_data$instrument[i] <- match_row$instrument
      } else if (nrow(match_row) > 1) {
        instrument <- unique(match_row$instrument)
        if (length(instrument) == 1) {
          temp_data$instrument[i] <- instrument
        } else {
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
  
  # Precalculate the number of PSMs by instrument
  total_counts_by_instrument <- as.data.frame(table(temp_data$instrument))
  names(total_counts_by_instrument) <- c("instrument", temp_data$peptidoform_id[1])
  write.csv(total_counts_by_instrument, 
            file = paste0("../out/total_counts_by_instrument_for_",temp_data$peptidoform_id[1], ".csv"),
            row.names =  FALSE)
  
  if (totals_only) {
    warning("csv written, but no plots generated.")
    return()
  }
  
  # Precalculate min and max calibrated error
  min_error <- round(min(temp_data$calibrated_error) - 2*mz_binwidth, digits = 2)
  max_error <- round(max(temp_data$calibrated_error) + 2*mz_binwidth, digits = 2)
  
  # Create bins to plot data by
  bins <- seq(min_error, max_error, by = mz_binwidth)
  
  # Cut the calibrated errors into bins
  temp_data$bin <- cut(temp_data$calibrated_error, bins, include.lowest = TRUE, labels = FALSE)
  
  if (!column_to_colour_by %in% names(temp_data)) {
    stop("The specified column to color by does not exist in the data frame.")
  }
  
  # Aggregate data for plotting
  formula_text <- as.formula(paste("cbind(count = calibrated_error) ~ bin +", column_to_colour_by))
  plot_data <- aggregate(formula_text, data = temp_data, FUN = length)
  
  # Get bin center for plotting
  plot_data$binstart <- bins[plot_data$bin]
  plot_data$binend <- bins[plot_data$bin+1]
  plot_data$bincenter <- (plot_data$binstart + plot_data$binend) / 2
  
  # Calculate total PSM counts per dataset
  total_counts <- plot_data %>%
    group_by(!!sym(column_to_colour_by)) %>%
    summarise(total = sum(count))
  
  # Convert to factor
  plot_data[[column_to_colour_by]] <- factor(plot_data[[column_to_colour_by]], levels = total_counts[[column_to_colour_by]])
  
  # Add to plot df and custom facet label
  plot_data <- plot_data %>%
    left_join(total_counts, by = column_to_colour_by) %>%
    mutate(facet_label = paste0(column_to_colour_by, "_total_PSMs:", total))
  
  head(plot_data)
  
  title_text <- gsub("^.*/|\\.csv$|^", "", absolute_file_path)
  title_text <- sub("^merged_data_", "", title_text)
  
  # Generate and display the plot with facets for each dataset_ID
  p <- ggplot(plot_data, aes(x = bincenter, y = count, fill = !!sym(column_to_colour_by))) +
    geom_col(width = mz_binwidth) +
    facet_wrap(~reorder(facet_label, -total), ncol = 3, labeller = label_parsed) +  
    scale_x_continuous(name = "Calibrated Error (m/z)") +
    labs(y = "Count", title = title_text) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = element_text(size = 10, hjust = 1)) +
    geom_vline(xintercept = -0.0095, linetype = "dashed", color = "black")
  
  if (column_to_colour_by == "instrument") {
    p <- p + scale_fill_manual(values = instrument_colors)
  } else {
    p <- p + scale_fill_brewer(palette = "Set3") 
  }
  
  print(p)
  return(p)
}

