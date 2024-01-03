import numpy as np
from sklearn.mixture import GaussianMixture
# from scipy.stats import norm
import csv
import pandas as pd
import ast   # For safely evaluating strings as Python literals
import matplotlib.pyplot as plt
# import random
import matplotlib.backends.backend_pdf
# import itertools

# Specify the CSV file path
csv_file = '../data/alldata_allbatches_aggregated_by_peptidoform_ID.csv'

# read in the data as a pandas .csv since this is where it originated as
aggregated_data = pd.read_csv(csv_file)
print(f"Total number of peptidoforms detected across all data: {len(set(aggregated_data['peptidoform_id']))}")

print(f"sanity check: are number of rows and number unique peptidoform IDs the same?\nnrows: {len(aggregated_data)}")
# Step 2: Convert 'dataset_ID' and 'calibrated_error' from string to list
aggregated_data['dataset_ID'] = aggregated_data['dataset_ID'].apply(ast.literal_eval)
aggregated_data['calibrated_error'] = aggregated_data['calibrated_error'].apply(ast.literal_eval)
aggregated_data['protein'] = aggregated_data['protein'].apply(ast.literal_eval)
aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(ast.literal_eval)
# Step 3: Filter data
# Keep rows with at least 3 unique dataset_IDs and at least 90 calibrated_error values
# 3 unique datasets suggested by Andy
# 50 was selected based on the following publication, which determines the number of values needed for
# different likelihood accuracies of GMMs (Table 1 and 2 of interest):
# https://www.sciencedirect.com/science/article/pii/S0031320319300445#:~:text=For%20instance%2C%20a%20sample%20size,only%202%20parameters%20are%20estimated).
threshold = 50
filtered_data = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= 3) &
                aggregated_data['calibrated_error'].apply(lambda x: len(x) >= threshold)]

print(f"Total number of peptidoforms post-filtering: {len(set(filtered_data['peptidoform_id']))}")

# Run a Gaussian mixed model - we expect 1-3 modes, but we might have an unexpected larger number.
# Therefore, we run the model with 1-3 components but if something weird pops up we might have to increase
# that number. However, larger number of components increases the risk of overfitting.
n_components_range = range(1, 3)
print(f"Fitting Gaussian mixed models with {n_components_range} components and selecting the best fit")

# we want to do the model fitting and selection for each peptidoform (key in filtered_dict)
# using the associated values as data to run the model on.

# first we will run the GMM for 1, 2, or 3 components, select the best model, then store that in a dictionary or list of models
#
# Convert string values to floats and reshape data for GMM
# for key in filtered_dict:
#     filtered_dict[key] = np.array(filtered_dict[key], dtype=float).reshape(-1, 1)

# # Print the names of the first 10 keys in filtered_dict
# first_10_keys = list(itertools.islice(filtered_dict.keys(), 10))
# print(first_10_keys)

# # Optionally, limit the dictionary to these 10 items for testing during further processing:
# filtered_dict = {key: filtered_dict[key] for key in first_10_keys}
# # print(filtered_dict)


# Dictionary to store the best model and its BIC for each peptidoform
best_models = {}

for peptidoform_id in filtered_data['peptidoform_id'].unique():
    # Extract calibrated errors for the current peptidoform_id
    errors = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].iloc[0]
    errors = np.array(errors, dtype=float).reshape(-1, 1)

    best_bic = np.inf
    best_gmm = None

    for n_components in n_components_range:
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(errors)
        bic = gmm.bic(errors)

        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm

    best_models[peptidoform_id] = {'model': best_gmm, 'bic': best_bic}

# now, best_models contains the best GMM for each peptidoform along with its BIC
# print(best_models)

# Filtering data to include peptidoforms with potential sulfonation
acceptable_range = (-0.015, -0.005)
sulfo_filtered_data = {}

print(f"Filtering data to only include peptiodoforms with potential sulfonation...")
## Select peptidoforms that exhibit a shift towards sulfonated.
for peptidoform_id, model_info in best_models.items():
    model = model_info['model']
    means = model.means_.ravel()

    if model.n_components == 1 and acceptable_range[0] <= means[0] <= acceptable_range[1]:
        sulfo_filtered_data[peptidoform_id] = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].iloc[0]
    elif model.n_components in [2, 3] and any(acceptable_range[0] <= mean <= acceptable_range[1] for mean in means):
        sulfo_filtered_data[peptidoform_id] = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].iloc[0]

# Now, sulfo_filtered_data contains the datasets that meet the GMM criteria

print(f"Writing results to sulfo_filtered_data.csv in the output folder...")
# Save sulfo_filtered_data to CSV
with open('../output/sulfo_filtered_data.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['PeptidoformID', 'Calibrated_mass_error'])  # Writing header
    for peptidoform_id, errors in sulfo_filtered_data.items():
        writer.writerow([peptidoform_id, ', '.join(map(str, errors))])

# # Randomly select 10 entries for histogram plotting just as a test
# random_keys = random.sample(list(sulfo_filtered_data.keys()), min(10, len(sulfo_filtered_data)))

print(f"Generating histograms for each peptidoform...")

# Create and save histograms in a PDF file
pdf = matplotlib.backends.backend_pdf.PdfPages(f"../output/sulfo_histograms_thresh{str(threshold)}.pdf")
for peptidoform_id, errors in sulfo_filtered_data.items():
    plt.figure()
    plt.hist(errors, bins=30, alpha=0.7)
    plt.title(f'Histogram for {peptidoform_id}')
    plt.xlabel('Calibrated_mass_error')
    plt.ylabel('Frequency')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
pdf.close()

print(f"Histograms saved in '../output/sulfo_histograms.pdf'")


######################## OLD CODE #########################


# # Open the CSV file and read the first n lines
# with open(csv_file, 'r') as file:
#     reader = csv.reader(file)
#
#     # Skip the header row - these are column names in our case
#     # next(reader, None)
#
#     # Read and process only the first 1000 lines for testing purposes
#     # n_lines = 1000
#     for row_num, row in enumerate(reader, start=1):
#         # print(row)
#
#         # if row_num > n_lines:
#         #     break  # Stop after the first n lines
# #
#         # column 1 contains unique IDs and column 2 contains calibrated mass error values
#         unique_id = row[0]
#         calibrated_errors = row[1]
#         dataset_ID = row[2]
#
#         # Check if the unique ID key already exists in the dictionary
#         if unique_id in data_dict:
#             # Append the value to the existing list
#             data_dict[unique_id].append(calibrated_errors)
#         else:
#             # Create a new key with the value as a list
#             data_dict[unique_id] = [calibrated_errors]

#
#