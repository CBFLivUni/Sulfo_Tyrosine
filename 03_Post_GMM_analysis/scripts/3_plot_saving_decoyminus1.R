# Start PDF device
pdf("../out/EnrichmentPlots_DECOY_minus1_abithacky.pdf", width=8.27, height=11.69)

# "Title PAGE"
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
text(0.5, 0.5, "GO Molecular Function, DECOY_minus1", cex = 2.5) 

goplot(cluster_GO_MF)  
barplot(cluster_GO_MF, showCategory = 30)
treeplot(cluster_GO_tree_MF)

# "Title PAGE"
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
text(0.5, 0.5, "GO Biological Process, DECOY_minus1", cex = 2.5) 
goplot(cluster_GO_BP)  
barplot(cluster_GO_BP, showCategory = 30)
treeplot(cluster_GO_tree_BP)


# "Title PAGE"
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
text(0.5, 0.5, "GO Cellular Component, DECOY_minus1", cex = 2.5) 
goplot(cluster_GO_CC)  
barplot(cluster_GO_CC, showCategory = 30)
treeplot(cluster_GO_tree_CC)


# "Title PAGE"
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
text(0.5, 0.5, "KEGG, DECOY_minus1_mz_-0.0175_-0.0125", cex = 2.5) 
ggplot(KEGG_Significant_Pathways, aes(y = reorder(Description, -pvalue), x = GeneRatio)) +
  geom_point(aes(size = Count, color = pvalue)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(y = "Description", x = "Gene Ratio", color = "pvalue", size = "Count") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("KEGG Pathways with p<0.05 but not p.adj")

dev.off()
