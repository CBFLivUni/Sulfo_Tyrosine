import ListFilesFun as LF
import os
import zipfile
from TPP_Process_Python import Convert_pepXML_toCSV as convert

folder_path = os.path.join(os.getcwd(), '../data_to_process')
print(folder_path)

extension_to_find = '.zip'
zip_files, prefixes = LF.getfiles(folder=folder_path, extension=extension_to_find, get_prefix=True)

for zip_file, prefix in zip(zip_files, prefixes):
    print(f"Processing {zip_file} with prefix {prefix}")

    # Create a unique output directory for each zip file
    unique_output_path = os.path.join(folder_path, 'extracted_files', prefix)
    if not os.path.exists(unique_output_path):
        os.makedirs(unique_output_path)

    # Unzip the contents into the unique output directory
    with zipfile.ZipFile(os.path.join(folder_path, zip_file), 'r') as zip_ref:
        # Extract each member to the correct location
        for member in zip_ref.namelist():
            # Construct the absolute path for each member
            target_path = os.path.join(unique_output_path, os.path.basename(member))
            # Extract the member only if it's a file (skip directories)
            if not member.endswith('/'):
                with zip_ref.open(member) as source, open(target_path, 'wb') as target:
                    target.write(source.read())

        # Process each .xml file in the unique output directory
        for root, dirs, files in os.walk(unique_output_path):
            for file in files:
                if file.endswith('.pep.xml'):
                    file_path = os.path.join(root, file)
                    convert.convert(input=file_path, file_prefix=prefix)
