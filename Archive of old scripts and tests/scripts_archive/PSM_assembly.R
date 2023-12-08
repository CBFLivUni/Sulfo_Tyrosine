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

## let's do it on one dataset first. 
test_data <- read.csv("/")