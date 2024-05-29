library(pdftools)


pdf_convert("../in/BarPlots_for_manuscript_multipleBOIs.pdf", format = "png", pages = NULL, filenames = NULL, dpi = 600, opw = "", upw = "", verbose = TRUE)



# individual aucs for bin of interest
# pages as follows: 

pages_15pc_individual <- c(100, # 100 = nASEEEPEYGEEIK_n230_1_S167_1_Y243_1 - Secretogranin-1 - KNOWN SULFATED
                           106, # 106 = nATWLSLFSSEESNLGANNYDDYR_n230_1_S167_3 - Vitronectin, KNOWN SULFATED
                           148, # 148 = nGDVFTMPEDEYTVYDDGEEK_n230_1_Y243_1 - Vitronectin, KNOWN SULFATED
                           
                           122, # 122 = nDSYETSQLDDQSAETHSHK_n230_1_S167_1_T181_1 - mix sulfated and phosphorylated; Osteopontin - secreted
                           43, # 43 = MSNYSLLSVDYVVDK_M147_1_S167_1_Y243_1  - sulfated only example; Sulfotransferase 2A1
                           398, # 398 = nVHNDAQSFDYDHDAFLGAEEAK_n230_1_S167_1_Y243_1 - sulfated; example for Calumenin - Golgi
                           
                           29, # 29 = KAYSFCGTVEYMAPEVVNR_M147_1_Y243_1 - false positive that slipped through filtering; 
                           40, # 40 =  LSRGSIDREDGSLQGPIGNQHIYQPVGKPDPAAPPK_S167_1 - false positive
                           420 # 420 = nYVGFGNTPPPQKK_n145_1_T181_1 - false positive
                          ) 



pdf_convert("../in/histograms_post_auc_filtering_15pc_individualAUCbin_-0.0125_-0.0075.pdf", 
            format = "png", 
            pages = pages_15pc_individual, 
            filenames = NULL, 
            dpi = 600, opw = "", upw = "", verbose = TRUE)
getwd()



# get another example for vitronectin, good for demonstrating the need for 2 gaussians
pdf_convert("../in/histograms_post_auc_filtering_15pc_individualAUCbin_-0.0175_-0.0125.pdf", 
            format = "png", 
            pages = 20, 
            filenames = NULL, 
            dpi = 600, opw = "", upw = "", verbose = TRUE)
