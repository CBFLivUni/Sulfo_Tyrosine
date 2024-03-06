library(pdftools)


pdf_convert("../in/BarPlots_peptidoform_distributions_byBin.pdf", format = "png", pages = NULL, filenames = NULL, dpi = 600, opw = "", upw = "", verbose = TRUE)



# individual aucs for bin of interest
# pages as follows: 

pages_15pc_individual <- c(43, 106, 122, 29)
# 43 = MSNYSLLSVDYVVDK_M147_1_S167_1_Y243_1  - sulfated only example; Sulfotransferase 2A1
# 106 = nATWLSLFSSEESNLGANNYDDYR_n230_1_S167_3 - doubly and singly sulfated mix; Vitronectin, known sY protein
# 122 = nDSYETSQLDDQSAETHSHK_n230_1_S167_1_T181_1 - mix sulfated and phosphorylated; Osteopontin - secreted
# 29 = KAYSFCGTVEYMAPEVVNR_M147_1_Y243_1 - false positive that slipped through filtering; 


pdf_convert("../in/histograms_post_auc_filtering_15pc_individualAUCbin_-0.0125_-0.0075.pdf", 
            format = "png", 
            pages = pages_15pc_individual, 
            filenames = NULL, 
            dpi = 600, opw = "", upw = "", verbose = TRUE)
