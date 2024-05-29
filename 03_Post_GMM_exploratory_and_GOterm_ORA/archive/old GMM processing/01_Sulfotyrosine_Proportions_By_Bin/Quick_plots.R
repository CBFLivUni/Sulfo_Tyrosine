w <- read.csv("../out/01_Sulfotyrosine_Proportions_By_Bin/wide_bins.csv", header=TRUE)
n <- read.csv("../out/01_Sulfotyrosine_Proportions_By_Bin/narrow_bins.csv", header=TRUE)
m <- read.csv("../out/01_Sulfotyrosine_Proportions_By_Bin/medium_bins.csv", header=TRUE)

data<-m
data <- data[order(data$lower_boundary, decreasing = FALSE), ]
data$bin_ID <- factor(data$bin_ID, levels = data$bin_ID)

names(data)[17:19] <- c("pY_containing_as_fraction_of_Y_containing", 
                        "pT_containing_as_fraction_of_T_containing", 
                        "pS_containing_as_fraction_of_S_containing")

# Open a PDF device to save all plots to a single file
pdf("All_Bar_Plots_medium_bins.pdf", width = 10, height = 8)

# Loop through each column starting from the fourth column
for (i in 4:ncol(data)) {
  # Create a bar plot for the tyrosine_containing_peptides_count against bin_ID
  p <- ggplot(data, aes(x = bin_ID, y = data[,i])) +
    geom_bar(stat = "identity", fill = "grey") +
    labs(x = "Bin ID", y = names(data)[i]) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Bar Plot of", names(data)[i], "vs Bin ID"))
  
  # Print the plot to the PDF device
  print(p)
}

# Close the PDF device
dev.off()
