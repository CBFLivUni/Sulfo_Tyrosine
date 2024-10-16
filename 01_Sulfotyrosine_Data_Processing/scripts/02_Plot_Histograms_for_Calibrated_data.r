### source functions and other scripts#####
source("plot_hist_fun.R") # plot all histograms - #plotting with function from Andy 

##### get datasets calibrated data to work with ##########
# get all calibrated tsvs
# wd <- paste0(getwd(), "/")
wd <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/01_Sulfotyrosine_Data_Processing/batch2_9/"
extension = ".pep_thresholded_calibrated.tsv"
calibrated_files <- list.files(wd, pattern = extension, full.names = FALSE, recursive = TRUE)
# extract the folder paths for hist plot f-n
folders <- dirname(calibrated_files)
# and the filenames for input in plot histograms f-n by Andy
input_filenames <- basename(calibrated_files)
########## plot initial hostograms ############
# for all files, ploit the histogram using andy's modified function 
## added a line or two of code in his f-n to deal with subdirectories
for (i in 16:length(input_filenames)) {
  print(input_filenames[[i]])
  plot_histograms(folder = folders[[i]],
                  input_file = input_filenames[[i]],
                  wd = wd)
}
