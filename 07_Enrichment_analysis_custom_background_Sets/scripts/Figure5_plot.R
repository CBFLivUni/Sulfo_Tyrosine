library(ggplot2)

setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/07_Enrichment_analysis_custom_background_Sets/scripts")

data_to_plot <- read.csv("../out/ORA/ORA_master_results_toplot.csv")



# convert Bin_Tag to a factor and set levels in the desired order
data_to_plot$Bin_Tag <- factor(data_to_plot$Bin_Tag, levels = c("DM1", "BOI", "D1", "D2", "D3"))

# split the 'Max_Genes' and 'Detected_Genes' into numeric values
data_to_plot$Background_detected <- as.integer(sub("/.*", "", data_to_plot$Max_Genes))
data_to_plot$Detected_in_bin <- as.integer(sub("/.*", "", data_to_plot$Detected_Genes))
data_to_plot$Background_total <- as.integer(sub(".*?/", "", data_to_plot$Max_Genes))
data_to_plot$Bin_total <- as.integer(sub(".*?/", "", data_to_plot$Detected_Genes))


# caculate actual ratios for plotting
data_to_plot$Detected_Ratio <- with(data_to_plot, Detected_in_bin / Bin_total)
data_to_plot$Max_Ratio <- with(data_to_plot, Background_detected / Background_total)

# assign significance levels based on qvalue to add as labels
data_to_plot$Significance <- ifelse(data_to_plot$qvalue <= 0.001, "***", 
                                    ifelse(data_to_plot$qvalue <= 0.01, "**", 
                                           ifelse(data_to_plot$qvalue <= 0.05, "*", "")))

# maintain custom color palette from before
custom_colors <- c("DM1" = "grey", "BOI" = "#E69F00", "D1" = "#56B4E9", "D2" = "#56B4E9", "D3" ="#56B4E9")



# generate all plots one by one for each term - to best visualise ratios; facet wrap did not work so well.


########## GOLGI ###############
golgi_data <- data_to_plot[data_to_plot$Category == "Golgi",]

ggplot(golgi_data, aes(x = Bin_Tag, y = Detected_Ratio, fill = Bin_Tag)) +
  geom_bar(stat = "identity", color = "black") + # This creates the colored portion
  geom_hline(aes(yintercept = Max_Ratio, group = Bin_Tag), size = 1, linetype = "dashed") + # add dashed line that represents ration in background
  scale_fill_manual(values = custom_colors) + # use the custom color palette
  geom_text(aes(label = paste0("q=",sprintf("%.3f", qvalue)), y = Detected_Ratio + 0.03), vjust = 1.5, color = "black", size = 5) + # add actual q values
  geom_text(aes(label = Significance, y = Detected_Ratio + 0.03), vjust = -0.5, color = "black", size = 5) + # add significance stars
  ggtitle("Custom Term: Golgi") +
  labs(x = "Bin", y = "Gene Ratio", fill = "Custom Tag") +
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # rotate and change size of x-axis labels
        axis.text.y = element_text(size = 12), # 
        axis.title.x = element_text(size = 14), # 
        axis.title.y = element_text(size = 14), # 
        plot.title = element_text(size = 16, hjust = 0.5)) # 
 


ggsave("../out/ORA/GOLGI.png", width = 5, height = 4, dpi = 600)


############## known sY ################

sY_data <- data_to_plot[data_to_plot$Category == "sY",]



ggplot(sY_data, aes(x = Bin_Tag, y = Detected_Ratio, fill = Bin_Tag)) +
  geom_bar(stat = "identity", color = "black") + # This creates the colored portion
  geom_hline(aes(yintercept = Max_Ratio, group = Bin_Tag), size = 1, linetype = "dashed") + # This adds the dashed line
  scale_fill_manual(values = custom_colors) + # Use the custom color palette
  geom_text(aes(label = paste0("q=",sprintf("%.3f", qvalue)), y = Detected_Ratio + 0.015), vjust = 1.5, color = "black", size = 5) +
  geom_text(aes(label = Significance, y = Detected_Ratio + 0.01), vjust = -0.5, color = "black", size = 5) + # Add significance stars
  ggtitle("Custom Term: Known sY") +
  labs(x = "Bin", y = "Gene Ratio", fill = "Custom Tag") +
  ylim(c(0,0.135)) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # rotate and change size of x-axis labels
        axis.text.y = element_text(size = 12), # 
        axis.title.x = element_text(size = 14), # 
        axis.title.y = element_text(size = 14), # 
        plot.title = element_text(size = 16, hjust = 0.5)) # 



ggsave("../out/ORA/sY.png", width = 5, height = 4, dpi = 600)


################ SECRETED #####################


secreted_data <- data_to_plot[data_to_plot$Category == "Secreted",]


# Create the plot with facet wrap by Category
ggplot(secreted_data, aes(x = Bin_Tag, y = Detected_Ratio, fill = Bin_Tag)) +
  geom_bar(stat = "identity", color = "black") + # This creates the colored portion
  geom_hline(aes(yintercept = Max_Ratio, group = Bin_Tag), size = 1, linetype = "dashed") + # This adds the dashed line
  scale_fill_manual(values = custom_colors) + # Use the custom color palette
  geom_text(aes(label = paste0("q=",sprintf("%.3f", qvalue)), y = Detected_Ratio + 0.07), vjust = 1.5, color = "black", size = 5) +
  geom_text(aes(label = Significance, y = Detected_Ratio + 0.07), vjust = -0.5, color = "black", size = 5) + # Add significance stars
  ggtitle("Custom Term: Secreted") +
  labs(x = "Bin", y = "Gene Ratio", fill = "Custom Tag") +
  theme_minimal() + # Minimal theme
  ylim(c(0, 0.6))+
  theme(legend.position = "none", # Remove legend
              axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # rotate and change size of x-axis labels
              axis.text.y = element_text(size = 12), # 
              axis.title.x = element_text(size = 14), # 
              axis.title.y = element_text(size = 14), # 
              plot.title = element_text(size = 16, hjust = 0.5)) # 

ggsave("../out/ORA/SECRETED.png", width = 5, height = 4, dpi = 600)


########## TRANSMEMBRANE ############


transmembrane_data <- data_to_plot[data_to_plot$Category == "Transmembrane",]


# Create the plot with facet wrap by Category
ggplot(transmembrane_data, aes(x = Bin_Tag, y = Detected_Ratio, fill = Bin_Tag)) +
  geom_bar(stat = "identity", color = "black") + # This creates the colored portion
  geom_hline(aes(yintercept = Max_Ratio, group = Bin_Tag), size = 1, linetype = "dashed") + # This adds the dashed line
  scale_fill_manual(values = custom_colors) + # Use the custom color palette
  geom_text(aes(label = paste0("q=",sprintf("%.3f", qvalue)), y = Detected_Ratio + 0.03), vjust = 1.5, color = "black", size = 5) +
  geom_text(aes(label = Significance, y = Detected_Ratio ), vjust = -0.5, color = "black", size = 5) + # Add significance stars
  ggtitle("Custom Term: Transmembrane") +
  labs(x = "Bin", y = "Gene Ratio", fill = "Custom Tag") +
  theme_minimal() + # Minimal theme
  ylim(c(0, 0.3))+
  theme(legend.position = "none", # Remove legend
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # rotate and change size of x-axis labels
        axis.text.y = element_text(size = 12), # 
        axis.title.x = element_text(size = 14), # 
        axis.title.y = element_text(size = 14), # 
        plot.title = element_text(size = 16, hjust = 0.5)) # 

ggsave("../out/ORA/TRANSMEMBRANE.png", width = 5, height = 4, dpi = 600)


############# unlikely sY ###########

unlikelysY_data <- data_to_plot[data_to_plot$Category == "Unlikely",]


# Create the plot with facet wrap by Category
ggplot(unlikelysY_data, aes(x = Bin_Tag, y = Detected_Ratio, fill = Bin_Tag)) +
  geom_bar(stat = "identity", color = "black") + # This creates the colored portion
  geom_hline(aes(yintercept = Max_Ratio, group = Bin_Tag), size = 1, linetype = "dashed") + # This adds the dashed line
  scale_fill_manual(values = custom_colors) + # Use the custom color palette
  geom_text(aes(label = paste0("q=",sprintf("%.3f", qvalue)), y = Detected_Ratio + 0.13), vjust = 1.5, color = "black", size = 5) +
  geom_text(aes(label = Significance, y = Detected_Ratio + 0.07), vjust = -0.5, color = "black", size = 5) + # Add significance stars
  ggtitle("Custom Term: Unlikely_sY") +
  labs(x = "Bin", y = "Gene Ratio", fill = "Custom Tag") +
  theme_minimal() + # Minimal theme
  ylim(c(0, 1.1))+
  theme(legend.position = "none", # Remove legend
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # rotate and change size of x-axis labels
        axis.text.y = element_text(size = 12), # 
        axis.title.x = element_text(size = 14), # 
        axis.title.y = element_text(size = 14), # 
        plot.title = element_text(size = 16, hjust = 0.5)) # 

ggsave("../out/ORA/Unlikely_sY.png", width = 5, height = 4, dpi = 600)
