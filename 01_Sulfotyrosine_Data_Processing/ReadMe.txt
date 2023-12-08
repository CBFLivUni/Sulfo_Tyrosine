This folder contains all scripts for processing .pep.xml data using Andy's TPP_Process_Python scripts.

Prior to processing, the .pep.xml data was downloaded using FileZilla from peptide atlas. Data was provided
by Eric Deutsch as a compressed .zip file. Each zipped file was extracted and the .pep.xml was moved to the 
parent folder, then all folders were moved to this folder. Processing was carried out as follows: 

1) Run 01_Data_Conversion.py to convert all pep.xml files to .tsv files.
2) Run 02_FDR_Thresholding_andRecalibration.py to re-calibrate the data and calculate FDR, then threshold based on FDR. 
3) Plot the histograms of the claibrated data by running 03_Plot_Histograms_for_Calibrated_data.r
