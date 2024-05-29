columns_to_plot <- names(results)[11:(ncol(results)-4)]
results <- results[, -c(21:26)]
# start PDF device
pdf("../out/BarPlots_for_manuscript_multipleBOIs.pdf", width = 11, height = 8.5)


# loop through each column to create and save the plot
for (col in columns_to_plot) {
  
  temp_results <- results %>%
    mutate(!!paste0(col, "_formatted") := sprintf("%.2f", !!sym(col)))
  
  # use the formatted column for the label in geom_text
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
  
  print(p)
}

dev.off()

########## chi squared residuals ###########

pdf("../out/ChiSquared_residuals.pdf", width = 11, height = 8.5)

# Loop over the columns and generate the plots
for (i in 4:ncol(observed_counts_for_chi)) {
  observed_col <- names(observed_counts_for_chi)[i]
  expected_col <- names(expected_counts_for_chi)[i]
  
  plot_chi_squared_residuals(observed_counts_for_chi, 
                             expected_counts_for_chi, 
                             observed_col, 
                             expected_col, 
                             pdf_device = NULL)
}

dev.off()



pdf("../out/ChiSquared_residuals_proportions.pdf", width = 11, height = 8.5)

for (i in 4:ncol(observed_proportions_for_chi)) {
  observed_col <- names(observed_proportions_for_chi)[i]
  expected_col <- names(expected_proportions_for_chi)[i]
  
  plot_chi_squared_residuals(observed_proportions_for_chi, expected_proportions_for_chi, observed_col, expected_col, pdf_device = NULL)
}

dev.off()


