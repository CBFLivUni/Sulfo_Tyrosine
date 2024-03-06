library("vcd")
library("gplots")
library("graphics")
# ideas - depending on the results for different factors - e.g. more sensitive instrument,
# older instruments, etc, we can reorder the data to better show the presence of patterns with these plots. 


png(filename = "Counts_by_Instrument_baloonplot.png", width = 4000, height = 4000, res = 300)

# Generate the balloonplot with rotated column labels
balloonplot(t(dt_reordered), main ="Counts by Instrument", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE,
            colsrt=90) # Rotate column labels to 90 degrees

# Close the device to save the file
dev.off()


# mosaic plot
# this would look better if data were in the cortrect order.

# Step 1: Extract the "D1" row
d1_row <- dt["D1", ]

# Step 2: Determine the order of columns based on descending values of the "D1" row
ordered_columns <- names(d1_row)[order(d1_row, decreasing = TRUE)]

# Step 3: Reorder the data table based on the determined order
dt_reordered <- dt[, ordered_columns]
png(filename = "Counts_by_Instrument_mosaicplot.png", width = 2000, height = 2000, res = 300)

mosaicplot(t(dt_reordered), shade = TRUE, las=2, main = "Counts by Instrument")

# Close the device to save the file
dev.off()


png(filename = "Counts_by_Instrument_mosaicplot_v2.png", width = 2000, height = 2000, res = 300)

mosaic(t(dt_reordered), shade = TRUE, las=2, main = "Counts by Instrument")

# Close the device to save the file
dev.off()


####
# install.packages("vcd")

# plot just a subset of the table
assoc(dt_reordered, shade = TRUE, las=2, main = "Counts by Instrument")






pdf(file = "residuals.pdf", width = 8, height = 6)

for (i in 1:ncol(standardised_residuals)){
  
  residuals_data <- data.frame(bin_ID = rownames(contingency_counts_for_goodness_of_fit), Residuals = standardized_residuals[,i])
  
  column_name <- colnames(standardized_residuals)[i]
  plot_title <- paste0("Standardized Residuals from Chi-squared Test: ", column_name)
  
  # Plotting the standardized residuals
  p <- ggplot(residuals_data, aes(x = bin_ID, y = Residuals)) +
    geom_point(size = 3) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") + # Thresholds for significance
    theme_minimal() +
    labs(title = plot_title,
         x = "Bin (Ordered by Lower m/z Boundary)",
         y = "Standardized Residual") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display the plot
  print(p)
  
}

dev.off()

