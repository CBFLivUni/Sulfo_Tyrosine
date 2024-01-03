import ListFilesFun as LF
import os
import zipfile
from TPP_Process_Python import Convert_pepXML_toCSV as convert
from TPP_Process_Python import CalculateFDR_and_threshold as Fdr  # run Andy's FDR thresholding
from TPP_Process_Python import Calibrate_errors as Cal  # calibration script from Andy

# use function to get a list of all already extracted pep.xml files that need to be converted
chunk_size_bytes = 8 * 1024 * 1024 * 1024  # 8GB in bytes # define chunk size when dealing with large zipped files

# Get the directory where data to process is and resolve absolute path
data_to_process_folder = os.path.abspath(os.path.join(os.getcwd(), '..', 'batch2_10'))
print(f"Processing initiated for all data found in: {data_to_process_folder}")

# Create folder to store new extracted files if needed
extracted_data_folder = os.path.join(data_to_process_folder, 'extracted_files')
if not os.path.exists(extracted_data_folder):
    print(f"Creating {extracted_data_folder}")
    os.makedirs(extracted_data_folder)

# Define file extensions for the data that needs processing
extension_to_find1 = '.pep.xml'
print(f"Looking for {extension_to_find1} files...")
# already extracted .pep.xml files
pep_xml_files, pep_prefixes = LF.getfiles(folder=data_to_process_folder, extension=extension_to_find1, get_prefix=True)
print(f"There are {len(pep_xml_files)} pep.xml files already extracted.")

# not yet extracted data
extension_to_find2 = '.zip'  # most actually originally as zipped.
print(f"Looking for {extension_to_find2} files...")
zip_files, zip_prefixes = LF.getfiles(folder=data_to_process_folder, extension=extension_to_find2, get_prefix=True)
print(f"There are {len(zip_files)} files that need to be extracted.")
#
# print("Processing already extracted files...")
# for pepxml_file, pepxml_prefix in zip(pep_xml_files, pep_prefixes):
#     print(f"Converting {pepxml_file} with prefix {pepxml_prefix}")
#     # make sure absolute paths are used and convert from pepxml to tsv
#     full_file_path = os.path.join(data_to_process_folder, pepxml_file)
#     # print(f"full path: {full_file_path}")
#     convert.convert(input=full_file_path, file_prefix=pepxml_prefix)
#
#     # Create a zip file for the .pep.xml file to save space, then delete the pepxml
#     zip_file_path = full_file_path + '.zip'
#     with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
#         zipf.write(full_file_path, os.path.basename(pepxml_file))
#     # Delete the original .pep.xml file after zipping to save space
#     os.remove(full_file_path)
# #
#
# print("Processing pep.xmls from within zipped folders...")
# for zip_file, zip_prefix in zip(zip_files, zip_prefixes):
#     print(f"Processing {zip_file} with prefix {zip_prefix}")
#
#     # Create a unique output directory for each zip file
#     unique_output_path = os.path.join(data_to_process_folder, 'extracted_files', zip_prefix)
#     if not os.path.exists(unique_output_path):
#         os.makedirs(unique_output_path)
#
#         # Unzip the contents into the unique output directory
#         with zipfile.ZipFile(os.path.join(data_to_process_folder, zip_file), 'r') as zip_ref:
#             # Extract each member to the correct location
#             for member in zip_ref.namelist():
#                 # Construct the absolute path for each member
#                 target_path = os.path.join(unique_output_path, os.path.basename(member))
#                 # Extract the member only if it's a file (skip directories)
#                 if not member.endswith('/'):
#                     with zip_ref.open(member) as source, open(target_path, 'wb') as target:
#                         while chunk := source.read(chunk_size_bytes):  # use chunk size specified at the top, currently 8GB
#                             target.write(chunk)
#
#             # Process each .xml file in the unique output directory
#             for root, dirs, files in os.walk(unique_output_path):
#                 for file in files:
#                     if file.endswith('.pep.xml'):
#                         file_path = os.path.join(root, file)
#                         convert.convert(input=file_path, file_prefix=zip_prefix)
#
#
# # Now run FDR thresholding
# # first get list of all files we just converted to .tsv
# extension_to_find3 = '.pep.tsv'
# tsvs = LF.getfiles(folder=data_to_process_folder, extension=extension_to_find3, get_prefix=False)
# print(f"There are {len(tsvs)} tsv files detected.")
#
# FDR_thresh = 0.01
# print(f"Thresholding based on FDR q threshold set to {FDR_thresh}")
# for tsv in tsvs:
#     print(f"Working on {tsv}...")
#     # make sure absolute paths are used and convert from pepxml to tsv
#     full_file_path = os.path.join(data_to_process_folder, tsv)
#     Fdr.calculateFDR(results_file=full_file_path, decoy_string="DECOY", q_thresh=FDR_thresh, prob_column_header="pp_prob")

# Now Recalibrate the data using the thresholded data only
print("Recalibrating data...")
extension_to_find4 = "thresholded.tsv"
tsvs_to_calibrate = LF.getfiles(folder=data_to_process_folder, extension=extension_to_find4, get_prefix=False)
for tsv in tsvs_to_calibrate:
    full_file_path = os.path.join(data_to_process_folder, tsv)
    Cal.calibrate_errors(full_file_path)
