library(tidyverse)
source("generate_USI_fun.R")
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/09_USI_generation/scripts")
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_BOI_peptidoforms_all_PSM_data/out/csv_tables_by_peptidoform/"

### we need to construct a USI for each peptidoform we are interested in
# online example: 
# mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_05_2Feb12_Cougar_11-10-09.mzML:scan:12298:[iTRAQ4plex]-LHFFM[Oxidation]PGFAPLTSR/3


# the format is: 
#   mzspec:<CollectionID>:<msRunComponent>:<IndexFlag>:<ScanNumber>:<SpectrumInterpretation>/<charge>


# get the files we are interested in read in. 
# we only want ones with convincing histogram and likely protein ID based on cellular localisation
# metadata <- read.csv("../in/potential_sulfo_peptidoforms_annotated.csv")

# use the updated table
metadata <- read.csv("../in/Table_S1_final_v2.csv")
peptidoforms_of_interst <- metadata[metadata$Histogram_Evaluation == "Convincing" &
                                      metadata$Biological_context_tags != "No prior knowledge", "peptidoform_ID"]
files_of_interst <- paste0(data_dir, "merged_data_", peptidoforms_of_interst, ".csv")


files_list <- list()

for (i in (1: length(files_of_interst))) {
  
  files_list[[i]] <- read.csv(files_of_interst[i])
  
}

test <- files_list[[1]]


files_with_USIs_list <- list()

for (i in (1: length(files_of_interst))) {
  
  files_with_USIs_list[[i]] <- generate_USI(files_list[[i]])
  
}

all_data <- files_with_USIs_list[[1]]

for (i in (2: length(files_of_interst))) {
  
  all_data <- rbind(all_data,files_with_USIs_list[[i]] )
  
}

# add columns that we will need

all_data$protein_ID <- ""
all_data$protein_name <- ""
all_data$known_sulfated <- ""



all_colnames <- colnames(all_data)
collnames_to_add_to_end <- all_colnames[1:24]

all_data_reorganised <- all_data[, c("peptidoform_id",
                                     "known_sulfated",
                                    
                                     "USI_workaround", 
                                     "USI_with_one_sulfo", 
                                     "USI_with_two_sulfo",
                                     "calibrated_error",
                                     
                                     "protein_ID",
                                     "protein_name",
                                     "dataset_ID",
                                     collnames_to_add_to_end)]



# populate protein ID and protein name based on matching peptidoform ID:
# if all_data_reorganised$peptidoform_id matches  metadata$peptidoform_id,
# 
# populate all_data_reorganised$protein_ID with values from the same row for metadata$cleaned_protein_IDs
# and populate all_data_reorganised$protein_name with values from the same row for metadata$Protein.names
# 


all_data_reorganised <- all_data_reorganised %>%
  left_join(metadata %>% 
              select(peptidoform_id, cleaned_protein_IDs, Protein.names), 
            by = "peptidoform_id") %>%
  mutate(protein_ID = cleaned_protein_IDs,
         protein_name = Protein.names) %>%
  select(-cleaned_protein_IDs, -Protein.names)



write.csv(all_data_reorganised, file = "../out/USIs_to_share_v3.csv", row.names = FALSE)

# use each row as a current data frame to feed into generate USI
new_USIs <- generate_USI(test)

# pick a spectrum to test that has a calibrated error close to the expected  -0.0095
spectra_to_choose_from <- test[-0.01 < test$calibrated_error & test$calibrated_error < -0.009,]

manual_example <- spectra_to_choose_from[7,] # pick random one

# extract dataset ID
current_dataset <- sub("([^-]*)(-.*)?", "\\1", manual_example$dataset_ID)

# i think ms run is in the spectrum column, but only includes the name up to the first '.'
msRun <- gsub("^(.*?)\\..*$", "\\1", manual_example$spectrum)

# the stuff after the first and second  . is the scan number

scan_number <- gsub("^[^.]*\\.([^.]+)\\..*$", "\\1", manual_example$spectrum)

# peptide/charge is the format the spectrum interpretation is needed in.
spectrum_interpretation <- paste0(manual_example$peptide, "/", manual_example$z)


USI <- paste("mzspec", current_dataset, msRun, "scan", scan_number, spectrum_interpretation, sep = ":")
USI
