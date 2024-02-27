### we need to manually assingn some of the instruments! for this reason we take
# the background data and take all relevent dataset_experimentTag combinations from there

data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/08_instrument_enrichment/in/"
all_bins_background <- read.csv(file = paste0(data_dir, "tyrosine_containing_background_inSwissProt.tsv"), 
                                header = TRUE,
                                sep = "\t"
)


metadata <- read.csv(file = paste0(data_dir, "human_phosphobuild_metadata.csv"), 
                     header = TRUE,
)

# dataset combinations - currently data is aggregated by peptiodofrm ID

dataset_combos <- all_bins_background$dataset_ID %>% unique()

# 8015 combos

# split at , and merge into a large list, but keep only unique records
all_datasets_split <- unique(unlist(strsplit(as.character(all_bins_background$dataset_ID), ", ")))

# 246 combos

expanded_data <- data.frame(Experiment_Tag = character(),
                            dataset_id = character(),
                            stringsAsFactors=FALSE)


# for each tag try to match it to the metadata
for(i in seq_along(all_datasets_split)) {
  # get the vector of dataset_experiemtnalTags
  
  datasets <- all_datasets_split[[i]]
  # populate the peptidoform id column 
  # peptidoform_id <- rep(BOI$peptidoform_id[i], length(datasets))
  
  # split each item at the first '-' - the bit on the left is dataset ID, everything else on the right including subsequent - is Experiment Tag 
  dataset_ids <- sapply(datasets, function(x) strsplit(x, "-", fixed=TRUE)[[1]][1])
  experiment_tags <- sapply(datasets, function(x) {
    parts <- strsplit(x, "-", fixed=TRUE)[[1]]
    paste(parts[-1], collapse="-")
  })
  
  temp_data <- data.frame( 
                          Experiment_Tag = experiment_tags,
                          dataset_id = dataset_ids,
                          stringsAsFactors=FALSE)
  
  expanded_data <- rbind(expanded_data, temp_data)
}

# set instrument name to NA
expanded_data$Instrument_Name <- NA

# try match the experiment tag and replace the NA where there is a match
for(i in 1:nrow(expanded_data)) {
  experiment_tag <- expanded_data$Experiment_Tag[i]
  matching_row <- metadata[metadata$Experiment.Tag == experiment_tag,]
  if(nrow(matching_row) > 0) {
    expanded_data$Instrument_Name[i] <- matching_row$Instrument.Name[1]
  }
}

# subset to only include the rows with NAs 

mismatched_datasets <- expanded_data[is.na(expanded_data$Instrument_Name),]
# save the data and manually deal with NAs. 
write.csv(mismatched_datasets, file = "../out/mismatched_datasets_for_manual_annotation.csv", row.names = TRUE)
