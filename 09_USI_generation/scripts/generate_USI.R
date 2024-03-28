library(tidyverse)

setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/09_USI_generation/scripts")
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/out/csv_tables_by_peptidoform/"

### we need to construct a USI for each peptidoform we are interested in
# online example: 
# mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_05_2Feb12_Cougar_11-10-09.mzML:scan:12298:[iTRAQ4plex]-LHFFM[Oxidation]PGFAPLTSR/3


# the format is: 
#   mzspec:<CollectionID>:<msRunComponent>:<IndexFlag>:<ScanNumber>:<SpectrumInterpretation>/<charge>


# get the files we are interested in read in. 
# we only want ones with convincing histogram and likely protein ID based on cellular localisation
metadata <- read.csv("../in/potential_sulfo_peptidoforms_annotated.csv")
peptidoforms_of_interst <- metadata[metadata$convincing_histogram_alldatasets == "y" &
                                      metadata$convincing_protein_subcellular_location == "y", "peptidoform_id"]
files_of_interst <- paste0(data_dir, "merged_data_", peptidoforms_of_interst, ".csv")


files_list <- list()

for (i in (1: length(files_of_interst))) {
  
  files_list[[i]] <- read.csv(files_of_interst[i])
  
}


files_with_USIs_list <- list()

for (i in (1: length(files_of_interst))) {
  
  files_with_USIs_list[[i]] <- generate_USI(files_list[[i]])
  
}

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
