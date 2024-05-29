library(tidyverse)
library(stringr) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(devtools)
library(ggplot2)

project_dir = "C:/Users/jtzve/Desktop/Sufo_Tyrosine/03_Post_GMM_analysis/"
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/03_Post_GMM_analysis/data/"
setwd(paste0(project_dir, "/scripts"))
source("GO_Enrichment_ORA_helper_functions.R")
########### read in data ###########
# specify the extension unique to your calibrated files
extension = "postGMM"

# get all files of interest
bin_data_files <- list.files(data_dir, pattern = extension, full.names = TRUE, recursive = TRUE)

# and get file path (not reltive!)
input_filenames <- basename(bin_data_files)

# read in all data one by one into an empty list, setting each data frame as a
# new elemnt of the list. 
postGMM_data <- list()
print("reading in data")
for (i in 1:length(bin_data_files)) {
  print(bin_data_files[[i]])
  postGMM_data[[i]] <- read.csv(bin_data_files[[i]])
  print(i/length(input_filenames)*100)
}

# read in all human proteins for data filtering - from Swiss Prot
Swiis_Prot_data_file <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/03_Post_GMM_analysis/data/SwissProtLibrary_2024_01_17.tsv"
Swiss_Prot_Human <- read.csv(file = Swiis_Prot_data_file, sep = "\t")


############ clean names of data #########
# clean the names - remove .csv, round the ranges, leave just the range and put mz before it
# reatin information about bin ID == clean names, lower boundary, upper boundary to then put in a dataframe

# Cleaning the names
cleaned_names <- gsub("postGMM_fitting_binrange", "mz", input_filenames)
cleaned_names <- gsub(".csv", "", cleaned_names, fixed = TRUE)

# Now, cleaned_names should have the format like 'mz_0.300-0.400' 

# Initialize empty vectors to store upper and lower boundaries
lower_boundaries <- numeric()
upper_boundaries <- numeric()


# Applying the function to each name and storing the results
cleaned_names <- sapply(cleaned_names, clean_round_extract, seq_along(cleaned_names))

names(postGMM_data) <- cleaned_names


############## remove pythonic remnants from data ##############
# we need to clean the read in data as it contains python chaacters (e.g. [] for a list)
# rm(cleaned_postGMM_data, test)
gc()
cleaned_postGMM_data <- list()
for (i in 1:length(postGMM_data)) {
  
  cleaned_postGMM_data[[i]] <- clean_data(postGMM_data[[i]])
  print(i/length(postGMM_data))
  
}

names(cleaned_postGMM_data) <- names(postGMM_data)
names(cleaned_postGMM_data)
########### filter data to only include Y-containing ###########
# for our foreground we only really want phosphotyrosine-containing peptides in 
# the bin of interest; reflect that in the background too!
## NB: PTM assignment might not be 100% correct and we want to actually include
# all tyrosine-containing peptidoforms, not just these assigned to have a 
# phosphotyrosine because some phosphogroups assigned to be positionally on a 
# Threonine or Serine might actually be on a Tyrosine. Therefore,
# for all data we are only interested in the tyrosinde containing peptidoforms 
BOI1 <- cleaned_postGMM_data[["mz_-0.0125_-0.0075"]]
BOI1_foreground_df <- BOI1[grepl("Y", BOI1$peptidoform_id), ]


BOI2 <- cleaned_postGMM_data[["mz_-0.0175_-0.0125"]]
BOI2_foreground_df <- BOI2[grepl("Y", BOI2$peptidoform_id), ]

BOI3 <- cleaned_postGMM_data[["mz_-0.0225_-0.0175"]]
BOI3_foreground_df <- BOI3[grepl("Y", BOI3$peptidoform_id), ]

D1 <- cleaned_postGMM_data[["mz_0.0075_0.0125"]]
D1_foreground_df <- D1[grepl("Y", D1$peptidoform_id), ]

D2 <- cleaned_postGMM_data[["mz_0.0125_0.0175"]]
D2_foreground_df <- D2[grepl("Y", D2$peptidoform_id), ]

D3 <- cleaned_postGMM_data[["mz_0.0175_0.0225"]]
D3_foreground_df <- D3[grepl("Y", D3$peptidoform_id), ]

# # for background we only want tyrosines but from across ALL bins, including bin of interest
# first get all bin data together and only keep unique rows (as some peptidoforms that could have bi or trimodal distribution may be captured in multiple bins)
tyrosines_containing_background_df <- bind_rows(cleaned_postGMM_data) %>% distinct()
# filter down to only include phosphotyrosine-containing peptides
tyrosines_containing_background_df <- tyrosines_containing_background_df[grepl("Y", tyrosines_containing_background_df$peptidoform_id), ]


data_for_ORA <- list(
  BOI1 = BOI1_foreground_df,
  BOI2 = BOI2_foreground_df,
  D1 = D1_foreground_df,
  D2 = D2_foreground_df,
  D3 = D3_foreground_df,
  background = tyrosines_containing_background_df
  
)
########## Retain only entires in SwissProt
# Split Swiss_Prot_Human$neXtProt and Swiss_Prot_Human$RefSeq into lists of IDs
# Apply the split_ids function I wrote; we will need these 
neXtProt_list <- lapply(Swiss_Prot_Human$neXtProt, split_ids)
RefSeq_list <- lapply(Swiss_Prot_Human$RefSeq, split_ids)

data_for_ORA_inSwissProt <- list()

for (i in 1:6) {
  
  data_for_ORA_inSwissProt[[i]] <- processProteinIDs(data_for_ORA[[i]],
                                                     Swiss_Prot_Human = Swiss_Prot_Human, 
                                                     neXtProt_list = neXtProt_list, 
                                                     RefSeq_list = RefSeq_list)
  
}

names(data_for_ORA_inSwissProt) <- names(data_for_ORA)

BOI3_inSwissProt <- processProteinIDs(BOI3_foreground_df,
                                      Swiss_Prot_Human = Swiss_Prot_Human, 
                                      neXtProt_list = neXtProt_list, 
                                      RefSeq_list = RefSeq_list)
# only 2 proteins, so not taken forward
data_for_ORA_inSwissProt[[2]]$cleaned_protein_IDs %>% unique() %>% length()
############## run enrichment analyses ###########

## convert the SwissProt IDs to Ensembl gene IDs
# Connect to the Ensembl database
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get the corresponding Entrez IDs for all UniProt IDs that we can
Uniprot = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot'), 
  mart = ensembl)


# convert swiss prot to ensembl gene ids and remove unmatched
background_genes <- extractEnsemblGeneIDs(data_for_ORA_inSwissProt[[6]], Uniprot)

# For foreground datasets do it inside a loop and run the analyses
ORA_results <- list()

for (i in 4:5) {
  
  foreground_genes <- extractEnsemblGeneIDs(data_for_ORA_inSwissProt[[i]], Uniprot)
  bin <- names(data_for_ORA_inSwissProt)[i]
  # run enrichment GO term analyses for MF, CC, and BP
  cluster_GO_MF <- enrichGO(gene         = foreground_genes,
                            universe     = background_genes,
                            OrgDb        = org.Hs.eg.db,  # human annotation
                            keyType      = "ENSEMBL", #just to make sure
                            ont          = "MF",  # 
                            pAdjustMethod = "BH", # 
                            pvalueCutoff = 1,
                            qvalueCutoff = 0.1)
  
 
  cluster_GO_CC <- enrichGO(gene         = foreground_genes,
                            universe     = background_genes,
                            OrgDb        = org.Hs.eg.db,  # 
                            keyType      = "ENSEMBL", #just to make sure
                            ont          = "CC",  #
                            pAdjustMethod = "BH", # 
                            pvalueCutoff = 1,
                            qvalueCutoff = 0.1)
  
  cluster_GO_BP <- enrichGO(gene         = foreground_genes,
                            universe     = background_genes,
                            OrgDb        = org.Hs.eg.db,  # human annotation
                            keyType      = "ENSEMBL", #just to make sure
                            ont          = "BP",  # 
                            pAdjustMethod = "BH", # 
                            pvalueCutoff = 1,
                            qvalueCutoff =  0.1)
  
  ## get the results df for each, assign an identifier column for the analysis type
  MF_results <- as.data.frame(cluster_GO_MF)
  try(MF_results$analysis_type <- "GO_MF", silent = TRUE)
  BP_results <- as.data.frame(cluster_GO_BP)
  try(BP_results$analysis_type <- "GO_BP", silent = TRUE)
  CC_results <- as.data.frame(cluster_GO_CC)
  try(CC_results$analysis_type <- "GO_CC", silent = TRUE)
  
  # join the results together and ass identifier for the analysis type
  results_df <- rbind(MF_results, BP_results, CC_results)
  try(results_df$bin <- bin, silent = TRUE)
  
 try(ORA_results[[bin]] <- results_df, silent = TRUE) 
}

# combine all results for plotting 

names(ORA_results)
df_for_final_plot <- as.data.frame(matrix(ncol = 0, nrow = 0))
for (result in ORA_results) {
  
  df_for_final_plot <- rbind(df_for_final_plot, result)
  
}

data_to_plot <- df_for_final_plot
# split the 'Max_Genes' and 'Detected_Genes' into numeric values
data_to_plot$Background_detected <- as.integer(sub("/.*", "", data_to_plot$BgRatio ))
data_to_plot$Detected_in_bin <- as.integer(sub("/.*", "", data_to_plot$GeneRatio))
data_to_plot$Background_total <- as.integer(sub(".*?/", "", data_to_plot$BgRatio))
data_to_plot$Bin_total <- as.integer(sub(".*?/", "", data_to_plot$GeneRatio))


# caculate actual ratios for plotting
data_to_plot$Detected_Ratio <- with(data_to_plot, Detected_in_bin / Bin_total)
data_to_plot$Max_Ratio <- with(data_to_plot, Background_detected / Background_total)

########## generate and save plot for Fig. 
head(data_to_plot$Detected_Ratio)
# The following columns have the data
# Description will be on the y axis, bin will be on the X axis, circle size will be 
# gene ratio, colour coded will be the p.adjust 

# make descriptionsmaller
data_to_plot$Description <- str_wrap(data_to_plot$Description, width = 60)

# make bin labels nicer
data_to_plot$bin <- gsub("BOI", "BOI_", data_to_plot$bin)
data_to_plot$bin <- gsub("D", "DECOY_", data_to_plot$bin)





###########  MF  ##########
MF_to_plot <- data_to_plot[data_to_plot$analysis_type == "GO_MF", c("Description",
                                                                    "bin", 
                                                                    "Detected_Ratio",
                                                                    "p.adjust",
                                                                    "analysis_type")]

# Arrange the data by bin and then by Detected_Ratio within each bin
MF_to_plot <- MF_to_plot %>%
  arrange(bin, Detected_Ratio) %>%
  mutate(Description = fct_inorder(Description))

# add D2 as NAs 
if(!"DECOY_2" %in% MF_to_plot$bin){
  MF_to_plot <- rbind(MF_to_plot, data.frame(Description="", bin="DECOY_2", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_MF"))
}

# Specify the order of the bins
MF_to_plot$bin <- factor(MF_to_plot$bin, levels = c("BOI_2", 
                                                    "BOI_1", 
                                                    "DECOY_1", 
                                                    "DECOY_2", 
                                                    "DECOY_3"))


p <- ggplot(MF_to_plot, aes(x = bin, y = Description)) + 
  geom_point(aes(size = Detected_Ratio, color = p.adjust)) +
  scale_color_gradient(low = "firebrick", high = "#FFCCCC") + 
  theme_minimal() +
  labs(
    title = "GO Molecular Function Significant Terms by Bin",
    size = "Gene Ratio",
    color = "p Adjusted"
  ) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1)
        )


print(p)
ggsave("../out/MF_GO_terms.png", width = 8, height = 8, units = "in", dpi = 600,
       plot = p)

######### BP ##########
BP_to_plot <- data_to_plot[data_to_plot$analysis_type == "GO_BP", c("Description",
                                                                    "bin", 
                                                                    "Detected_Ratio",
                                                                    "p.adjust",
                                                                    "analysis_type")]

# Arrange the data by bin and then by Detected_Ratio within each bin
BP_to_plot <- BP_to_plot %>%
  arrange(bin, Detected_Ratio) %>%
  mutate(Description = fct_inorder(Description))

# add D2 as NAs 
if(!"DECOY_2" %in% BP_to_plot$bin){
  BP_to_plot <- rbind(BP_to_plot, data.frame(Description="", bin="DECOY_2", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_BP"))
}

if(!"BOI_1" %in% BP_to_plot$bin){
  BP_to_plot <- rbind(BP_to_plot, data.frame(Description="", bin="BOI_1", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_BP"))
}

if(!"BOI_2" %in% BP_to_plot$bin){
  BP_to_plot <- rbind(BP_to_plot, data.frame(Description="", bin="BOI_2", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_BP"))
}
# Specify the order of the bins
BP_to_plot$bin <- factor(BP_to_plot$bin, levels = c("BOI_2", 
                                                    "BOI_1", 
                                                    "DECOY_1", 
                                                    "DECOY_2", 
                                                    "DECOY_3"))


p <- ggplot(BP_to_plot, aes(x = bin, y = Description)) + 
  geom_point(aes(size = Detected_Ratio, color = p.adjust)) +
  scale_color_gradient(low = "firebrick", high = "#FFCCCC") + 
  theme_minimal() +
  labs(
    title = "GO Biological Process Significant Terms by Bin",
    size = "Gene Ratio",
    color = "p Adjusted"
  ) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1)
  )


print(p)
ggsave("../out/BP_GO_terms.png", width = 8, height = 5, units = "in", dpi = 600,
       plot = p)

###########  CC  ##########
CC_to_plot <- data_to_plot[data_to_plot$analysis_type == "GO_CC", c("Description",
                                                                    "bin", 
                                                                    "Detected_Ratio",
                                                                    "p.adjust",
                                                                    "analysis_type")]

# Arrange the data by bin and then by Detected_Ratio within each bin
CC_to_plot <- CC_to_plot %>%
  arrange(bin, Detected_Ratio) %>%
  mutate(Description = fct_inorder(Description))

# add D2 as NAs 
if(!"DECOY_2" %in% CC_to_plot$bin){
  CC_to_plot <- rbind(CC_to_plot, data.frame(Description="", bin="DECOY_2", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_CC"))
}

if(!"BOI_1" %in% CC_to_plot$bin){
  CC_to_plot <- rbind(CC_to_plot, data.frame(Description="", bin="BOI_1", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_CC"))
}
# Specify the order of the bins
CC_to_plot$bin <- factor(CC_to_plot$bin, levels = c("BOI_2", 
                                                    "BOI_1", 
                                                    "DECOY_1", 
                                                    "DECOY_2", 
                                                    "DECOY_3"))


p <- ggplot(CC_to_plot, aes(x = bin, y = Description)) + 
  geom_point(aes(size = Detected_Ratio, color = p.adjust)) +
  scale_color_gradient(low = "firebrick", high = "#FFCCCC") + 
  theme_minimal() +
  labs(
    title = "GO Cellular Component Significant Terms by Bin",
    size = "Gene Ratio",
    color = "p Adjusted"
  ) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1)
  )


print(p)
ggsave("../out/CC_GO_terms.png", width = 8, height = 4, units = "in", dpi = 600,
       plot = p)

##### all terms by rbinding ####
all_terms_to_plot_v2 <- rbind(MF_to_plot, CC_to_plot, BP_to_plot)

p <- ggplot(all_terms_to_plot_v2, aes(x = bin, y = Description)) + 
  geom_point(aes(size = Detected_Ratio, color = p.adjust)) +
  scale_color_gradient(low = "firebrick", high = "#FFCCCC") + 
  theme_minimal() +
  labs(
    title = "GO Cellular Component Significant Terms by Bin",
    size = "Gene Ratio",
    color = "p Adjusted"
  ) +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1)
  )


print(p)
ggsave("../out/ALL_GO_terms_v2.png", width = 8, height = 11, units = "in", dpi = 600,
       plot = p)

# ########## combined plot of all terms ##########
# all_terms_to_plot <- data_to_plot[, c("Description",
#                                                                     "bin", 
#                                                                     "Detected_Ratio",
#                                                                     "p.adjust",
#                                                                     "analysis_type")]
# 
# # Arrange the data by bin and then by Detected_Ratio within each bin
# all_terms_to_plot <- all_terms_to_plot %>%
#   arrange(bin, Detected_Ratio) %>%
#   mutate(Description = fct_inorder(Description))
# 
# # add D2 as NAs 
# if(!"DECOY_2" %in% all_terms_to_plot$bin){
#   all_terms_to_plot <- rbind(all_terms_to_plot, data.frame(Description="", bin="DECOY_2", Detected_Ratio=NA, p.adjust=NA, analysis_type="GO_CC"))
# }
# 
# 
# # Specify the order of the bins
# all_terms_to_plot$bin <- factor(all_terms_to_plot$bin, levels = c("BOI_2", 
#                                                     "BOI_1", 
#                                                     "DECOY_1", 
#                                                     "DECOY_2", 
#                                                     "DECOY_3"))
# 
# all_terms_to_plot$Description <- str_wrap(all_terms_to_plot$Description, width = 50)
# 
# p <- ggplot(all_terms_to_plot, aes(x = bin, y = Description)) + 
#   geom_point(aes(size = Detected_Ratio, color = p.adjust)) +
#   scale_color_gradient(low = "firebrick", high = "#FFCCCC") + 
#   theme_minimal() +
#   labs(
#     title = "GO All Significant Terms by Bin",
#     size = "Gene Ratio",
#     color = "p Adjusted"
#   ) +
#   theme(legend.position = "right",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1)
#   )
# 
# 
# print(p)
# ggsave("../out/ALL_GO_terms.png", width = 8, height = 11, units = "in", dpi = 600,
#        plot = p)
