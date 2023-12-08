######### load libraries ##########










# now we need to assemble the data from two datasets! 














## let's do it manually on two dataset first: 
# PXD037549 is one of the smaller ones # DATASET1 #############
test_data <- read.csv("../Sulfotyrosine/PXD037549/PXD037549_interact-prob.pep_calibrated.tsv", sep = "\t")

# lets subset to only include the peptide sequence data for now; check if we can assign row names to spectrum name
length(unique(test_data$spectrum)) == nrow(test_data)
# TRUE => we can set rown
rownames(test_data) <- test_data$spectrum
colstokeep <- c("ppm_error", "da_error", "calibrated_error", "peptide", "mod_peptide")

sub <- test_data[, colstokeep]

# filter to only include peptides that have a phosphorylated amino acid
sub_filtered <- FilterPTMs(sub, 
                           pepcols = colstokeep,
                           PTMs_of_interest = c("S[167]", "T[181]", "Y[243]"))

# get a summary od calibrated Da shift by each unique nonstrict  peptidoform (e.g. PTM position not taken into account) 
sub_caldasum <- Summarise_By_Peptidoform(sub_filtered, 
                                         columns_to_summarise = c("calibrated_error"), # columns to summarise
                                         pepcols = c("peptide", "mod_peptide")) # columns with peptide data

########## DATASET2 - PXD036069, also quite small ##################
test_data2 <- read.csv("../Sulfotyrosine/PXD036069/PXD036069_interact-prob.pep_calibrated.tsv", sep = "\t")

# lets subset to only include the peptide sequence data for now; check if we can assign row names to spectrum name
length(unique(test_data2$spectrum)) == nrow(test_data2)
# TRUE => we can set rown
rownames(test_data2) <- test_data2$spectrum
sub2 <- test_data2[, colstokeep]

# filter to only include peptides that have a phosphorylated amino acid
sub_filtered2 <- FilterPTMs(sub2, 
                           pepcols = colstokeep,
                           PTMs_of_interest = c("S[167]", "T[181]", "Y[243]"))
# much fewer kept here! 

# get a summary od calibrated Da shift by each unique nonstrict  peptidoform (e.g. PTM position not taken into account) 
sub_caldasum2 <- Summarise_By_Peptidoform(sub_filtered2, 
                                         columns_to_summarise = c("calibrated_error"), # columns to summarise
                                         pepcols = c("peptide", "mod_peptide")) # columns with peptide data

## lets see how many of one data set overlap with the other? 
sum(sub_caldasum$peptidoform_id %in% sub_caldasum2$peptidoform_id)
# only 33!!!! 

# lets set rownames to the peptidoform ID and merge the two 
rownames(sub_caldasum) <- sub_caldasum$peptidoform_id
rownames(sub_caldasum2) <- sub_caldasum2$peptidoform_id
# ok, lets bind the data based on rows and investigate where the overlap is - not really of interest.
# Perform inner join to keep only overlapping rows
# merged_inner <- inner_join(sub_caldasum, sub_caldasum2, by = "peptidoform_id")


# Perform full join to combine all data 
merged_full <- full_join(sub_caldasum, sub_caldasum2, by = "peptidoform_id")

# Function to combine and flatten the list of calibrated mass shifts
flatten_lists <- function(x, y) {
  c(unlist(x), unlist(y))
}

# Aggregate the results
# Assuming 'calibrated_error_all_values' is the column with the list of calibrated mass shifts
merged_full_agg <- merged_full %>%
  group_by(peptidoform_id) %>%
  summarise(calibrated_mass_shifts = list(flatten_lists(calibrated_error_all_values.x,
                                                        calibrated_error_all_values.y)))




