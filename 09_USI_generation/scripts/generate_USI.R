library(tidyverse)

setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/09_USI_generation/scripts")
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/06_histograms_of_peptidoforms_of_interest_per_dataset/out/csv_tables_by_peptidoform/"

### we need to construct a USI for each peptidoform we are interested in

# the format is: 
#   mzspec:<CollectionID>:<DataSource>:<ScanType>:<ScanNumber>:<IndexType>=<IndexValue>:<Interpretation>


# online example: 
# mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_05_2Feb12_Cougar_11-10-09.mzML:scan:12298:[iTRAQ4plex]-LHFFM[Oxidation]PGFAPLTSR/3


# manually constructed example by me:
# read in a test dataset
test<- read.csv("../../06_histograms_of_peptidoforms_of_interest_per_dataset/out/csv_tables_by_peptidoform/merged_data_nIYEFPETDDEEENK_n145_1_T181_1.csv")
# pick a spectrum == Abemaciclib_01201_C07_P012160_B00_A00_R1.02646.02646.2
hist(test$calibrated_error, breaks = 50)
manual_example <- test[7,]
current_dataset <- sub("([^-]*)(-.*)?", "\\1", manual_example$dataset_ID)
current_data_source <- manual_example$spectrum # this may be incorrect
scan_number_string <- manual_example$native_id
# Use regmatches and regexpr to extract the scan number
scan_number <- as.numeric(regmatches(scan_number_string, regexpr("(?<=scan=)\\d+", scan_number_string, perl = TRUE)))


paste("mzspec", current_dataset, current_data_source, "scan", scan_number, sep = ":")
