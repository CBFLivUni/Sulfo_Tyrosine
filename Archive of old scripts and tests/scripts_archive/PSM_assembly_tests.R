## first Discard peptides that do not have at least one phospho modification to 
# speed up follow-up computaiton 

## next, we need to assemble all peptide-spectrum matches (PSMs) for a given peptidoform 
# not strict (i.e. position does not matter)
# we would therefore need to keep track of the peptide sequence and number of 
# modifications per modification type
# record all Da mass shift for the given peptide in the calibrated data. plot this 
# as histogram to capture the distribution - if normally distributed, it is likely phospho,
# if shofted to left, it is likely sulfor misidentified as phospho,
# if small bump to left, it is likely sulfo in some cases, phospho in other cases. 


### ideally this will be a function that works on a list of data frames, each of
# which contains the calibrated data from a dataset. 

## let's do it on one dataset first: # PXD037549 is one of the smaller ones
test_data <- read.csv("../Sulfotyrosine/PXD037549/PXD037549_interact-prob.pep_calibrated.tsv", sep = "\t")

head(test_data)

# lets subset to only include the peptide sequence data for now; check if we can assign row names to spectrum name

length(unique(test_data$spectrum)) == nrow(test_data)
# TRUE => we can
 rownames(test_data) <- test_data$spectrum

 sub <- test_data[, c("ppm_error", "da_error", "calibrated_error", "peptide", "mod_peptide")]

 ## I need to get an understanding of what mdoifications are within the data
 library(stringr)
# # Extract all modifications as a vector
# modifications <-str_extract_all(sub$mod_peptide, "\\[[0-9]+\\]") %>% 
#                 unlist() %>% # flatten list to vector 
#                 unique() %>% # get only unique ones
#                 sort() # sort them
# 
# print(modifications) #print to inspect 
# # this is all good, but not very informative - the number represents the amino acid + modification mass


modifications_withamino <-str_extract_all(sub$mod_peptide, "[A-Za-z]\\[[0-9]+\\]") %>% 
  unlist() %>% # flatten list to vector 
  unique() %>% # get only unique ones
  sort() # sort them

print(modifications_withamino) #print to inspect 
# [1] "A[151]" "C[143]" "E[111]" "M[147]" "N[115]" "n[43]"  "Q[111]" "Q[129]" "S[167]" "T[181]" "Y[243]"

# of these, "S[167]" "T[181]" "Y[243]" correspond to phosphorylated serine (S), threonine (T), and tyrosine (Y)

# Alanine 151 is also the correct mass for phosphorylation or sulfation, so keep for now.

