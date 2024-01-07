library(stringr)

wd <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/01_Sulfotyrosine_Data_Processing/scripts"
setwd(wd)

# read in the aggregated data
test <- read.csv("../out/_agg")
names(test)
test1<-test
# Convert the string representations of lists to actual lists and then apply 
# unique to see where we're actually going wrong - are we recording all information in the
# aggregated_by_peptidoform_ID.csv files? 

# testing some stuff to write cleaning function in loop:
# clean_string <- gsub("\\[|\\]|'", "", test[3,5]) 
# clean_string <- gsub(";", ",", test[2,5])# Remove brackets and single quotes
# elements <- unlist(strsplit(clean_string, ',\\s*'))  # Split into elements
# unique(elements)
# length(unique(elements))
# i = 3


# for columns that are not the peptidoform ID one
for (i in 2:ncol(test)){
  # Apply a cleaning function to:
  cleaned_col <- lapply(test[,i], function(x) {
    
    # 1) # Remove brackets and single quotes from the string
    clean_string <- gsub("\\[|\\]|'", "", x)  
    
    # 2) replace ; with , for the alt protein column (number 5)
    clean_string <- gsub(";", ",", clean_string)
    
    # 3) # Split into elements at every comma removing any white space to follow it
    elements <- unlist(strsplit(clean_string, ',\\s*'))
    
    # 4) return all elements
    return(elements)
      
  })
  
  
  # 5) Convert the list of vectors into a single character string per row
  # by collapsing the elements
  # for the calibrated mass errors column (number 2), retain all elements.
  if (i==2){
    
    test[,i] <- sapply(cleaned_col, function(x) paste(x, collapse = ", "))
    
  } else {
  # for all other columnsm reatin only the unique values
  test[,i] <- sapply(cleaned_col, function(x) paste(unique(x), collapse = ", "))
    
  }
  
}


# sanity check: how many peptidoforms were found in more than one dataset? 
# easist way is to count how many rows contain a ',' in their 3rd column
has_comma <- grepl(",", test$dataset_ID)
sum(has_comma)

## conclusion - b the looks of it our code that outputs all aggregated_by_peptidoform_ID.csv
# funcitons as intended. The trouble we run into must be an artefact of subsequent
# scripts.