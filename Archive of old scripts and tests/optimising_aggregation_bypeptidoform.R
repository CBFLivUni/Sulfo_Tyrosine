source("PTM_filtering_fun.R") # filter peptides by PTMs to only include phosphorylated - my function
source("Summarise_by_peptidoform_fun2.R") # for each dataset summarise the mass shift for each unique peptidoform - my function


##### get datasets calibrated data to work with ##########
# get all calibrated tsvs
wd <- paste0(getwd(), "/")
extension = ".pep_calibrated.tsv"
calibrated_files <- list.files(wd, pattern = extension, full.names = FALSE, recursive = TRUE)
# extract the folder paths for hist plot f-n
folders <- dirname(calibrated_files)
# and the filenames for input in plot histograms f-n by Andy
input_filenames <- basename(calibrated_files)


# read in all data one by one into an empty list, setting each data frame as a
# new elemnt of the list. 
calibrated_data <- list()

for (i in 1:length(calibrated_files)) {
  
  calibrated_data[[i]] <- read.csv(calibrated_files[[i]], sep = "\t")
  
}

# we want to add names to the dataframes so they are used e.g. when writing csvs. 
# first strip everythign after and including _interact .... .tsv
stripped_filenames <- gsub("_interact.*", "", input_filenames)


########## filter all data to only include peptides with phospho PTMs ##########
#determine which columns of the data you would like to continue workinbg with.
colstokeep <- c("ppm_error", "da_error", "calibrated_error", "peptide", 
                "mod_peptide") # must have this col to filter!

# determine what PTMs are of interest
PTMs_of_interest = c("S[167]", "T[181]", "Y[243]") # here phosphorylation of S, T, and Y

#both of the above are default settings, left here for clarity/readability

# for every calibrated dataset
filtered_data <- list()
for (i in 1:length(calibrated_data)) {
  #apply the PTM filtering function to only include peptides with phosphtylation
  filtered_data[[i]] <- FilterPTMs(df = calibrated_data[[i]],
                                   pepcols = colstokeep,
                                   PTMs_of_interest = PTMs_of_interest)
}


####### summarise data by unique peptidoform ##########
# define which parameters you want summarised by peptidoform
columns_of_interest <- c("calibrated_error") # e.g. cols for histogram plotting
pepcols <-  c("peptide", "mod_peptide") # cols with peptide and PTM data
# these are defaults, here just for ease
peptidoform_id_summaries <- list()
names(filtered_data) <- stripped_filenames

for (i in 1:length(filtered_data)) {
  #apply the PTM filtering function to only include peptides with phosphtylation
  output_name <- names(filtered_data[i])
  print(output_name)
  peptidoform_id_summaries[[i]] <- Summarise_By_Peptidoform(df = filtered_data[[i]],
                                                            columns_to_summarise =  columns_of_interest,
                                                            pepcols = pepcols,
                                                            outname = output_name)
  
}


peptidoform_id_fun2 <- peptidoform_id_summaries

peptidoform_id_fun3 <- peptidoform_id_summaries 
# fun3 is ~ 13% more memory efficient for storing the output, haven't tested for 
# memory efficiency while running. 




# Function to safely get values from a vector (handles NA)
safe_get <- function(vec, key) {
  if (is.na(vec[key]) || is.null(vec[key])) {
    return("")
  } else {
    return(vec[key])
  }
}

# Function to merge two dict_vectors
merge_two_dict_vectors <- function(dict_vector_1, dict_vector_2) {
  all_keys <- unique(c(names(dict_vector_1), names(dict_vector_2)))
  merged_dict_vector <- sapply(all_keys, function(key) {
    paste(c(safe_get(dict_vector_1, key), safe_get(dict_vector_2, key)), collapse = ",")
  })
  names(merged_dict_vector) <- all_keys
  return(merged_dict_vector)
}

# Merge all dict_vectors in the list

test_list <- peptidoform_id_fun2[1:3]
merged_dict_vector <- Reduce(merge_two_dict_vectors, test_list)

