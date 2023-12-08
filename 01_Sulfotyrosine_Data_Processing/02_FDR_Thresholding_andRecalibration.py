import ListFilesFun as LF # my function to get a list of files based on extension across all folders of the parent directory
import os
from TPP_Process_Python import CalculateFDR_and_threshold as FDR # run Andy's FDR thresholding
from TPP_Process_Python import Calibrate_errors as Cal # calibration script from Andy

folder_path = os.getcwd()  # Get the current directory
extension_to_find = '.pep.tsv'
tsvs = LF.getfiles(folder = folder_path, extension = extension_to_find, get_prefix = False)
print(f"there are {len(tsvs)} tsv files detected, as follows: \n {tsvs}")

print(tsvs[0])
# FDR.calculateFDR_Peptide_Level(results_file = tsvs[0], decoy_string = "DECOY", q_thresh = 0.01, prob_column_header = "pp_prob")
# running the above results in the error when I was testing it on a single file

# run FDR thresholding and calibration on all files
for tsv in tsvs:
    print(f"Calibrating {tsv}")
    FDR.calculateFDR(results_file = tsv, decoy_string = "DECOY", q_thresh = 0.01, prob_column_header = "pp_prob")
    Cal.calibrate_errors(tsv)