import ListFilesFun as LF
import os
from TPP_Process_Python import Convert_pepXML_toCSV as convert # run andy's conversion script

# use function to get a list of all pep.xml files that need to be converted
folder_path = os.getcwd()  # Get the current directory
print(folder_path)
extension_to_find = '.pep.xml'
pep_xmls, prefixes = LF.getfiles(folder = folder_path, extension = extension_to_find, get_prefix = True)
# print(len(pep_xmls))
# print(prefixes)
for file, prefix in zip(pep_xmls, prefixes):
    print(f"Converting {file} with prefix {prefix}")
    convert.convert(input = file, file_prefix = prefix)


