import ListFilesFun as Lf
import os
from TPP_Process_Python import Convert_pepXML_toCSV as Convert

# Update the folder path to point to the 'data_to_process' subdirectory, relative to the 'scripts' directory
folder_path = os.path.join(os.getcwd(), '../data_to_process')
print(folder_path)

extension_to_find = '.pep.xml'
# Use the modified function to get a list of all pep.xml files that need to be converted
pep_xmls, prefixes = Lf.getfiles(folder=folder_path, extension=extension_to_find, get_prefix=True)

print(pep_xmls)
print(prefixes)

for file, prefix in zip(pep_xmls, prefixes):
    print(f"Converting {file} with prefix {prefix}")
    # The file path needs to be combined with the folder_path to get the absolute path
    absolute_file_path = os.path.join(folder_path, file)
    Convert.convert(input=absolute_file_path, file_prefix=prefix)
