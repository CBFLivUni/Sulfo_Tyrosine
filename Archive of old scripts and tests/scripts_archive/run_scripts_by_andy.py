# import andy's scripts for xml file conversion to csv and data preprocessing
import os
from TPP_Process_Python import Convert_pepXML_toCSV as convert

# from TPP_Process_Python import CalculateFDR_and_threshold as threshold
# from TPP_Process_Python import Calibrate_errors as calibrate


# test by running on one file
current_folder = os.getcwd()
print(current_folder)
# xml_file_path = "/CPTAC_S013/interact-ipro-ptm.pep.xml"

# Run the conversion Python script on one file
# convert.convert(input = xml_file_path, file_prefix = "CPTAC_S013")
# getting error no such file or directory:
# FileNotFoundError: [Errno 2] No such file or directory: '/CPTAC_S013/CPTAC_S013_interact-ipro-ptm.pep.tsv'
# perhaps I need full path and not relative path? run a few tests to see if the path exists:


# print(f"testing if path exists for {xml_file_path}")
# if os.path.exists(xml_file_path):
#     print("path OK")
# else:
#     print("Path incorrect, does not exist")
# path did not exist, had to remove / from it, the one below works
xml_file_path = "PXD000222/PXD000222/U2OS_PLKinh-30m_TiO2-Diemthyl/comet/interact-ipro.pep.xml"

# now running the converison function
convert.convert(input = xml_file_path, file_prefix = "PXD000222")
# inspected tsv, looks good! now automate for multiple files.