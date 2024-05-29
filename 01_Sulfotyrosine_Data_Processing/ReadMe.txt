This folder contains all scripts for processing .pep.xml data to .tsv using Andy's TPP_Process_Python scripts.

Prior to processing, the .pep.xml data was downloaded using FileZilla from PeptideAtlas. Data was provided
by Eric Deutsch as a compressed .zip file. Each zipped file was extracted and the .pep.xml was moved to the 
parent folder, then all folders were moved to this folder. 

Processing was carried using the scripts in the 'scripts' folder as follows: 

1) Run 01_Data_Processing.py to convert all pep.xml files to .tsv files, perform FDR thresholding, and re-calibrate
the Da error mass shift

2) Plot the histograms of the claibrated data by running 02_Plot_Histograms_for_Calibrated_data.r

3) Filter data to include phosphorylated peptidoforms only and generate non-strict peptidofotm IDs by running
03_Generate_peptidoform_IDs.Rmd

4) Aggregate all data by peptidoform ID within a batch, then aggregate the resulting .csv files across all batches 
by running 04_Aggregate_batch_data_bypeptidoformID.py and specifying the corresponding batch folder in the script.
- per batch .csf files are in the out/new_peptidoformIDs folder
- the final aggregated data is in the out/ folder saved under clean_new_peptidoform_IDs_aggregated_by_peptidoform_ID.csv
- the out/ folder also contains some excluded experiments and reasoning for why these were excluded. 