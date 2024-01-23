library(biomaRt)
library(org.Hs.eg.db)
# load data from 03_postGMManalyses

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("ensembldb", "AnnotationDbi"))


## we get the fasta files for all background cleaned IDs with mathc in SwissPRot
proteins_to_download_FASTA_aaseq_for <- flattened_background_cleaned_IDs[flattened_background_cleaned_IDs %in% Swiss_Prot_Human$Entry]
#make sure we inly keep only unique

proteins_to_download_FASTA_aaseq_for <- unique(proteins_to_download_FASTA_aaseq_for)

protein_ids_df <- data.frame(ID = proteins_to_download_FASTA_aaseq_for)
write.table(protein_ids_df, "../data/protein_ids_for_FASTA_download.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE)




# we run a powerSHell script to collec the associated FASTA sequences from the web



# create mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Get the corresponding ENSEMBL IDs for all UniProt IDs that we can
Uniprot = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot'), 
  mart = ensembl)

## get the genes in our DECOY bin and save each to a file

for (id in foreground_genes) {
  sequence <- getSequence(id = id, type = "ensembl_gene_id", seqType = "peptide", mart = mart)
  if (nrow(sequence) > 0) {
    file_name <- paste0("../out/",id, ".fasta")
    writeLines(paste(">", id, "\n", sequence$sequence, sep = ""), file_name)
  }
}

test <- foreground_genes[230:240]

id <- test[1]
sequence <- getSequence(id = id, type = "ensembl_gene_id", seqType = "peptide", mart = ensembl)
# this returns a df with multiple sequences

longest_sequence <- sequence$peptide[7]


for (i in 1:nrow(sequence)) {
  if (i != 7) {  # Assuming 7 is the index of the longest sequence
    is_contained <- grepl(pattern = sequence$peptide[i], x = longest_sequence)
    print(paste("Sequence", i, "is contained in the longest sequence:", is_contained))
  }
}


if (nrow(sequence) > 0) {
  file_name <- paste0("../out/",id, ".fasta")
  writeLines(paste(">", id, "\n", sequence$sequence, sep = ""), file_name)
}


  
for (id in foreground_genes) {
  sequence <- getSequence(id = id, type = "ensembl_gene_id", seqType = "peptide", mart = mart)
  if (nrow(sequence) > 0) {
    file_name <- paste0("../out/",id, ".fasta")
    writeLines(paste(">", id, "\n", sequence$sequence, sep = ""), file_name)
  }
}
