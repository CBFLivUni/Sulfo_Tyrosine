columns_to_plot <- names(results)[11:(ncol(results)-4)]
results <- results[, -c(21:26)]
# start PDF device
pdf("../out/BarPlots_byBin_15pc_summed_final_version_for_manuscript_BLANK_LABELS.pdf", width = 11, height = 8.5)


# loop through each column to create and save the plot
for (col in columns_to_plot) {
  # Create a temporary data frame with an additional column for the formatted labels
  temp_results <- results %>%
    mutate(!!paste0(col, "_formatted") := sprintf("%.2f", !!sym(col)))
  
  # Use the formatted column for the label in geom_text
  p <- ggplot(temp_results, aes(x = reorder(bin_ID, lower_boundary), y = get(col))) +
    geom_bar(stat = "identity", aes(fill = highlight)) +
    scale_fill_manual(values = color_blind_friendly_colors) +
    geom_text(aes(label = get(paste0(col, "_formatted")), color = highlight), vjust = 0.5, 
              size = 11, fontface = "plain", angle = 0, nudge_y = 0.03) +
    scale_color_manual(values = color_blind_friendly_colors) +
    theme_minimal() +
    labs(x = "Bin ID", y = col, title = paste(col, "by Bin")) +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      plot.title = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),  
      axis.text.y = element_text(angle = 0, hjust = 1, size = 24),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.5, "cm"),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(title = "Category"), color = FALSE)
  
  print(p) # Print the plot
}
# Close the PDF device
dev.off()


# 
# ########### boxplots ##########
# # Start PDF device
# pdf("../out/BoxPlots_byBin.pdf", width = 11, height = 8.5)
# 
# for (col in columns_to_plot) {
#   # Subset the dataframe to include only the current column and bin_ID
#   data_to_plot <- data.frame(bin_ID = results$bin_ID, value = results[[col]], PeptideType = col)
#   
#   # Generate the boxplot
#   p <- ggplot(data_to_plot, aes(x = bin_ID, y = value)) +
#     geom_boxplot() +
#     theme_minimal() +
#     labs(x = "Bin ID", y = "Proportion/Fraction", title = paste("Distribution of", col, "Across Bins")) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
#   # Print the plot
#   print(p)
# }
# 
# dev.off()






pdf("../out/scatterplots_bymz_2.pdf", width = 11, height = 8.5)


for (col_name in y_columns) {
  # Create the scatter plot for the filtered dataset
  p <- ggplot(results, aes_string(x = "lower_boundary", y = col_name)) +
    geom_point(aes(color = highlight), size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black", aes_string(x = "lower_boundary", y = col_name)) +
    geom_text_repel(aes_string(label = "bin_ID"), size = 2.5, box.padding = 0.35, point.padding = 0.5) +
    scale_color_manual(values = color_blind_friendly_colors) +
    theme_minimal() +
    labs(x = "Lower m/z Boundary of bin", y = col_name, 
         title = paste("Scatter Plot of", col_name, "vs. Lower Boundary")) +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display the plot
  print(p)
}

dev.off()

pdf("../out/scatterplots_bymz_first3rowsremoved.pdf", width = 11, height = 8.5)

# Sort results by 'lower_boundary' and exclude the three rows with the lowest 'lower_boundary'
sorted_results <- results[order(results$lower_boundary), ]
filtered_results <- sorted_results[-c(1:3), ] # Remove the first three rows after sorting - they stood out as outliers,
# likely due to low number of peptidoforms in therE? 

for (col_name in y_columns) {
  # Create the scatter plot for the filtered dataset
  p <- ggplot(filtered_results, aes_string(x = "lower_boundary", y = col_name)) +
    geom_point(aes(color = highlight), size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "black", aes_string(x = "lower_boundary", y = col_name)) +
    geom_text_repel(aes_string(label = "bin_ID"), size = 2.5, box.padding = 0.35, point.padding = 0.5) +
    scale_color_manual(values = color_blind_friendly_colors) +
    theme_minimal() +
    labs(x = "Lower m/z Boundary of bin", y = col_name, 
         title = paste("Scatter Plot of", col_name, "vs. Lower Boundary")) +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display the plot
  print(p)
}

dev.off()



pdf("../out/ChiSquared_residuals.pdf", width = 11, height = 8.5)

# Loop over the columns and generate the plots
for (i in 4:ncol(observed_counts_for_chi)) {
  observed_col <- names(observed_counts_for_chi)[i]
  expected_col <- names(expected_counts_for_chi)[i]
  
  plot_chi_squared_residuals(observed_counts_for_chi, expected_counts_for_chi, observed_col, expected_col, pdf_device = NULL)
}

# Close the PDF device
dev.off()



pdf("../out/ChiSquared_residuals_proportions_2.pdf", width = 11, height = 8.5)

# Loop over the columns and generate the plots
for (i in 4:ncol(observed_proportions_for_chi)) {
  observed_col <- names(observed_proportions_for_chi)[i]
  expected_col <- names(expected_proportions_for_chi)[i]
  
  plot_chi_squared_residuals(observed_proportions_for_chi, expected_proportions_for_chi, observed_col, expected_col, pdf_device = NULL)
}

# Close the PDF device
dev.off()



# fishers test


pdf("../out/FishersExactTest_pvalues_2.pdf", width = 11, height = 8.5)

# List of peptide types to test
peptide_types <- c("Y", "T", "S", "pY", "pT", "pS")

# Loop over the peptide types and generate the plots
for (peptide in peptide_types) {
  observed_col <- paste("observed", peptide, "counts", sep = "_")
  expected_col <- paste("expected", peptide, "counts", sep = "_")
  
  plot_fishers_exact_test(observed_counts_for_fisher, expected_counts_for_fisher, peptide, pdf_device = NULL)
}

# Close the PDF device
dev.off()
