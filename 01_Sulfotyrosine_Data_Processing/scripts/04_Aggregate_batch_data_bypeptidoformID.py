import ListFilesFun as LF
import os
import pandas as pd

# Step 1: Get all CSV file paths
batch = 'new_peptidoform_IDs'
# Get the directory where peptidoform_id data is stored
data_to_process_folder = os.path.abspath(os.path.join(os.getcwd(), '..', 'out', batch))

print(data_to_process_folder)
extension_to_find = 'peptidoforms.csv'
# extension_to_find = '_ID.csv'
csv_files = LF.getfiles(folder=data_to_process_folder, extension=extension_to_find, get_prefix=False)
print(len(csv_files))
# we need absolute paths:
csv_file_paths = [os.path.join(data_to_process_folder, csv_file) for csv_file in csv_files]

# Step 2: Read and concatenate all CSV files
all_data = pd.concat([pd.read_csv(file) for file in csv_file_paths], ignore_index=True)

# Step 3: Group by 'peptidoform_id' and aggregate data; keep only unique values for all columns but calibrated errors
aggregated_data = all_data.groupby('peptidoform_id').agg({
    'calibrated_error': list,
    'dataset_ID': lambda x: list(set(x)),
    'protein': lambda x: list(set(x)),
    'alt_protein': lambda x: list(set(x))
}).reset_index()

# Step 4: Save the results as a CSV file
aggregated_data.to_csv(f'../out/{batch}_aggregated_by_peptidoform_ID.csv', index=False)
