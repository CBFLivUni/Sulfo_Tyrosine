library(patchwork)

# create blank plot for when no hits w/ q<= 0.1
blank_plot <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "No hits with q<=0.1", size = 5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "black")) # white background


############  BOI ################
# run code in qmd, then generate plot
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_MF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_MF$GeneRatio <- gene_ratio

BOI_MF <- ggplot(GO_term_dataframe_MF, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  labs(
    title = "BOI_GO_MF",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(hjust = 1, size = 20, face = "bold"),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16)
  )

# make consistent across
BOI_MF <- BOI_MF + scale_size(range = c(3, 8)) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 1, face = "plain", size = 10),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )
  
# in figure i have each plot at 4.5 * 5.6 cm, but this is too small to fit the plot? 
width_in_inches <- 5.6 
height_in_inches <- 4.5 
BOI_MF
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/BOI_GO_MF.png", 
       plot = BOI_MF, 
       width = width_in_inches, 
       height = height_in_inches, units = "in",
       dpi = 600)


############ Decoy 1 ######################

##### MF ######
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_MF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_MF$GeneRatio <- gene_ratio


GO_term_dataframe_MF$Description <- str_wrap(GO_term_dataframe_MF$Description, width = 50)


D1_MF <- ggplot(GO_term_dataframe_MF, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_MF",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(angle = 0,hjust = 1, size = 10),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D1_MF
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D1_GO_MF.png", 
       plot = D1_MF, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)





####### BP ##########
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_BP$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_BP$GeneRatio <- gene_ratio

GO_term_dataframe_BP$Description <- str_wrap(GO_term_dataframe_BP$Description, width = 50)


D1_BP <- ggplot(GO_term_dataframe_BP, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_BP",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12, face = "plain", lineheight = 1),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D1_BP
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D1_GO_BP.png", 
       plot = D1_BP, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)




####### CC ##########
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_CC$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_CC$GeneRatio <- gene_ratio

GO_term_dataframe_CC$Description <- str_wrap(GO_term_dataframe_CC$Description, width = 50)


D1_CC <- ggplot(GO_term_dataframe_CC, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_CC",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12, face = "plain", lineheight = 1),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D1_CC
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D1_GO_CC.png", 
       plot = D1_CC, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)






############# DECOY 3 #################


##### MF ######
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_MF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_MF$GeneRatio <- gene_ratio


GO_term_dataframe_MF$Description <- str_wrap(GO_term_dataframe_MF$Description, width = 50)


D3_MF <- ggplot(GO_term_dataframe_MF, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_MF",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(angle = 0,hjust = 1, size = 10),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D3_MF
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D3_GO_MF.png", 
       plot = D3_MF, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)




####### CC ##########
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_CC$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_CC$GeneRatio <- gene_ratio

GO_term_dataframe_CC$Description <- str_wrap(GO_term_dataframe_CC$Description, width = 50)


D3_CC <- ggplot(GO_term_dataframe_CC, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_CC",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12, face = "plain", lineheight = 1),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D3_CC
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D3_GO_CC.png", 
       plot = D3_CC, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)



####### BP ##########
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_BP$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_BP$GeneRatio <- gene_ratio

GO_term_dataframe_BP$Description <- str_wrap(GO_term_dataframe_BP$Description, width = 30)


D3_BP <- ggplot(GO_term_dataframe_BP, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_BP",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12, face = "plain", lineheight = 1),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

D3_BP
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/D3_GO_BP.png", 
       plot = D3_BP, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)








############ DECOY MINUS 1 ##################


##### MF ######
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_MF$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_MF$GeneRatio <- gene_ratio


GO_term_dataframe_MF$Description <- str_wrap(GO_term_dataframe_MF$Description, width = 50)


DM1_MF <- ggplot(GO_term_dataframe_MF, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_MF",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(angle = 0,hjust = 1, size = 10),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

DM1_MF
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/DM1_GO_MF.png", 
       plot = DM1_MF, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)



####### CC ##########
gene_ratio <- as.numeric(sapply(strsplit(GO_term_dataframe_CC$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
GO_term_dataframe_CC$GeneRatio <- gene_ratio

GO_term_dataframe_CC$Description <- str_wrap(GO_term_dataframe_CC$Description, width = 50)


DM1_CC <- ggplot(GO_term_dataframe_CC, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() + # Create the dots
  scale_color_gradient(low = "red", high = "blue") + # Color gradient from low (good) to high (bad) p-values
  theme_bw() +
  scale_size(range = c(3, 8)) +
  labs(
    title = "BOI_GO_CC",
    x = "Gene Ratio",
    size = "Count",
    color = "p.adjusted"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(hjust = 1, size = 12, face = "plain", lineheight = 1),
    axis.ticks.length.x = unit(2, "mm"),
    axis.ticks = element_line(linewidth = 1),
    axis.title.x = element_blank(),
    # plot.title = element_text(size = 12),
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10)
  )

DM1_CC
ggsave("../out/15pc_sum_AUC_results/plots_to_show/pathway enrichment/DM1_GO_CC.png", 
       plot = DM1_CC, 
       width = width_in_inches , 
       height = height_in_inches, units = "in",
       dpi = 600)



############## blanks ############

blank_plot

combined_plots <- blank_plot + blank_plot + blank_plot + 
  blank_plot + blank_plot + blank_plot +
  
  plot_layout(nrow = 52, ncol = 3)



combined_plots

a4_width_in <- 210 / 25.4
a4_height_in <- 297 / 25.4

ggsave("combined_plots_a4.png", combined_plots, width = a4_width_in, height = a4_height_in, dpi = 300)



