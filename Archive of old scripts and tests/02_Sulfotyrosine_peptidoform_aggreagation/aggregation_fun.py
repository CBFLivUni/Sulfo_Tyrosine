import ListFilesFun as LF
import os
import pandas as pd

# Step 1: Get all CSV file paths
folder_path = os.getcwd()
extension_to_find = '_peptidoforms.csv'
#
csv_files = LF.getfiles(folder = folder_path, extension = extension_to_find, get_prefix = False)
print(csv_files)

# Step 2: Create an empty dictionary to store the data
peptidoform_errors = {}

# Step 3: Process each CSV file
for file in csv_files:
    df = pd.read_csv(file)
    print(f"working on {file}")
    for index, row in df.iterrows():
        peptidoform_id = row['peptidoform_id']
        calibrated_error = row['calibrated_error']

        if peptidoform_id not in peptidoform_errors:
            peptidoform_errors[peptidoform_id] = []
        peptidoform_errors[peptidoform_id].append(calibrated_error)

# Step 4: Convert the dictionary to a df
data = []
for peptidoform_id, errors in peptidoform_errors.items():
    for error in errors:
        data.append([peptidoform_id, error])
result_df = pd.DataFrame(data, columns=['peptidoform_id', 'calibrated_error'])


# Step 5: Save the results as a CSV file
# this is n
result_df.to_csv('master_peptidoform_errors_v2.csv', index=False)