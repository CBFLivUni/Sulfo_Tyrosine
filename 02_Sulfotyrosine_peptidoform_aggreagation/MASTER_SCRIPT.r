### source functions and other scripts#####
source("plot_hist_fun.R") # plot all histograms - #plotting with function from Andy 
source("PTM_filtering_fun.R") # filter peptides by PTMs to only include phosphorylated - my function
source("Summarise_by_peptidoform_fun.R") # for each dataset summarise the mass shift for each unique peptidoform - my function
source("PSM_assembly_fun.R") # function to assemble the summarised data across multiple datasets


##### get datasets calibrated data to work with ##########
# get all calibrated tsvs
wd <- paste0(getwd(), "/")
extension = ".pep_calibrated.tsv"
calibrated_files <- list.files(wd, pattern = extension, full.names = FALSE, recursive = TRUE)
# extract the folder paths for hist plot f-n
folders <- dirname(calibrated_files)
# and the filenames for input in plot histograms f-n by Andy
input_filenames <- basename(calibrated_files)

########## plot initial hostograms ############
# for all files, ploit the histogram using andy's modified function 
## added a line or two of code in his f-n to deal with subdirectories
# for (i in 1:length(input_filenames)) {
#   
#   plot_histograms(folder = folders[[i]],
#                   input_file = input_filenames[[i]],
#                   wd = wd)
# }


# read in all data one by one into an empty list, setting each data frame as a
# new elemnt of the list. 
calibrated_data <- list()
print("reading in data")
for (i in 1:length(calibrated_files)) {
  print(calibrated_files[[i]])
  calibrated_data[[i]] <- read.csv(calibrated_files[[i]], sep = "\t")
  
}
# we want to add names to the dataframes so they are used e.g. when writing csvs. 
# first strip everythign after and including _interact .... .tsv
stripped_filenames <- gsub("_interact.*", "", input_filenames)
# assign this to our list
names(calibrated_data) <- stripped_filenames


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
  print(i)
  filtered_data[[i]] <- FilterPTMs(df = calibrated_data[[i]],
                                   pepcols = colstokeep,
                                   PTMs_of_interest = PTMs_of_interest)
}


# add names just in case
names(filtered_data) <- stripped_filenames
####### summarise data by unique peptidoform ##########
# define which parameters you want summarised by peptidoform
columns_of_interest <- c("calibrated_error") # e.g. cols for histogram plotting
pepcols <-  c("peptide", "mod_peptide") # cols with peptide and PTM data
# these are defaults, here just for ease
peptidoform_id_summaries <- list()

for (i in 1:length(filtered_data)) {
  #apply the PTM filtering function to only include peptides with phosphtylation
  output_name <- names(filtered_data[i])
  print(output_name)
  peptidoform_id_summaries[[i]] <- Summarise_By_Peptidoform(df = filtered_data[[i]],
                                   columns_to_summarise =  columns_of_interest,
                                   pepcols = pepcols,
                                   outname = output_name)
  
}

names(peptidoform_id_summaries) <- stripped_filenames

########### Assemble PSM data ############
# function we need to merge the outputs of peptidoform_id_summaries

######### code below needs optimising! 

  # Use reduce with full_join to join all data frames in the peptidoform summary list
  # this results in a tibble with peptidoform id in the first column and one #
  # column for each summary in the input
  merged_full <- reduce(peptidoform_id_summaries, full_join, by = "peptidoform_id") %>% 
                 mutate_all( ~replace(., lengths(.)==0, NA)) #replace NULL values with NAs



# trying data.table to speed it up because it might be faster than dplyr>? 
library(data.table)
setDT(merged_full)


concatenate_values <- function(...) {
  # Combine all arguments into a single vector
  vals <- c(...)
  
  # Initialize an empty vector to store the results
  result <- character(0)
  
  # Iterate over each value
  for (val in vals) {
    # Check if the value is NA or not a character, and skip if it is
    if (is.na(val) || !is.character(val)) {
      next
    }
    
    # Split comma-separated values and unlist them
    val_split <- unlist(strsplit(val, split = ",", fixed = TRUE), use.names = FALSE)
    
    # Remove NA and trim spaces
    val_clean <- trimws(val_split[!is.na(val_split)])
    
    # Append to the result
    result <- c(result, val_clean)
  }
  
  return(result)
}


merged_full[, peptidoform_id := unlist(peptidoform_id)]
cols_to_convert <- grep("^calibrated_error", names(merged_full), value = TRUE)
merged_full[, (cols_to_convert) := lapply(.SD, as.character), .SDcols = cols_to_convert]
# this takes ages.... I am still figuring out if there's a way to speed it up, ChatGPT gave this to me :D 
result <- merged_full[, .(concatenated = concatenate_values(.SD)), by = peptidoform_id, .SDcols = cols_to_convert]
