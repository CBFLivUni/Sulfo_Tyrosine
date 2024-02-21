
sY_known <- read.csv("../data/known_sY_SP_2024_02_20.tsv", header = TRUE, sep = "\t")

sY_potential <- read.csv("../data/potential_sY_SP_2024_02_20.tsv", header = TRUE, sep = "\t")

## searching for secreted on swiss prot was apparently not enough to get all proteins with signal peptide. 
has_signal_peptide <- sY_potential[sY_potential$Signal.peptide != "",]
# the above list comprises 3382 proteins base don a subset of trnasmembrane or secreted as part of the search terms. 

# 3613 proteins in swissprot have a signal peptide
has_signal_peptide2 <- Swiss_Prot_Human[Swiss_Prot_Human$Signal.peptide != "",]

# quick look at differences
# setdiff(has_signal_peptide2$Protein.names, has_signal_peptide$Protein.names)
# setdiff(has_signal_peptide2$Subcellular.location..CC., has_signal_peptide$Subcellular.location..CC.)


sY_potential$Subcellular.location..CC.[1:15]
names(sY_potential)

# custom protein sets
gene_sets <- list(
  confirmed_sY = sY_known$Entry,
  potential_sY_general = sY_potential$Entry,
  potential_sY_Golgi = c("geneA", "geneB", "geneC"),
  potential_sY_transmembrane = sY_potential$Transmembrane[sY_potential$Transmembrane != ""],
  potential_sY_secreted = sY_potential$Subcellular.location..CC.["Secreted" %in% sY_potential$Subcellular.location..CC.],
  potential_sY_signal_peptide = sY_potential$Signal.peptide[sY_potential$Signal.peptide != ""]
  unlikely_sY_SwissProt = Swiss_Prot_Human
)

# Create a function to format each gene set
format_gmt_line <- function(set_name, genes) {
  paste(set_name, "NA", paste(genes, collapse = "\t"), sep = "\t")
}

# Apply the function to each gene set and combine into a single character vector
gmt_lines <- unlist(lapply(names(gene_sets), function(x) format_gmt_line(x, gene_sets[[x]])))

# Write the GMT lines to a file
writeLines(gmt_lines, "./custom_gene_sets.gmt")