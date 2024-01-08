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
csv_file = '../data/clean_new_peptidoform_IDs_aggregated_by_peptidoform_ID.csv'

# read in the data as a pandas .csv since this is where it originated as
aggregated_data = pd.read_csv(csv_file)
print(f"Total number of peptidoforms detected across all data: {len(set(aggregated_data['peptidoform_id']))}")

print(f"sanity check: are number of rows and number unique peptidoform IDs the same?\nnrows: {len(aggregated_data)}")


# Step 2: Convert 'dataset_ID' and 'calibrated_error' from string to list
aggregated_data['dataset_ID'] = aggregated_data['dataset_ID'].apply(ast.literal_eval)
aggregated_data['calibrated_error'] = aggregated_data['calibrated_error'].apply(ast.literal_eval)
aggregated_data['protein'] = aggregated_data['protein'].apply(ast.literal_eval)

# for alt protein, sometimes there isn't one so we have nan(s); this throws an error when evaluated
# lets replace 'nan' and 'nan, ' with ''
strings_to_replace = ['nan, ', 'nan,', 'nan']  # quick and dirty
for string in strings_to_replace:
    aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(lambda x: x.replace(string, ""))

aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(ast.literal_eval)
# now this works fine.
print("Data read in complete.")
# Step 3: Filter data
# Keep rows with at least 3 unique dataset_IDs and at least 90 calibrated_error values
# 3 unique datasets suggested by Andy
# 90 was selected based on the following publication, which determines the number of values needed for
# different likelihood accuracies of GMMs (Table 1 and 2 of interest):
# https://www.sciencedirect.com/science/article/pii/S0031320319300445#:~:text=For%20instance%2C%20a%20sample%20size,only%202%20parameters%20are%20estimated).

# lets do it step by step and see how much data we 'lose'
n_datasets_threshold = 3
filtered_by_dataset_ID = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= n_datasets_threshold)]
print(f"Total number of peptidoforms post filtering by dataset ID: {len(set(filtered_by_dataset_ID['peptidoform_id']))}")

n_calerrors_threshold = 90
filtered_data = filtered_by_dataset_ID[filtered_by_dataset_ID['calibrated_error'].apply(lambda x: len(x) >= n_calerrors_threshold)]
print(f"Total number of peptidoforms post subsequent filtering by number of calibrated errors: {len(set(filtered_data['peptidoform_id']))}")
print(f"min_n_datasets: {n_datasets_threshold} \n min_n_masses: {n_calerrors_threshold}")
# test code up until here with all data.
# we started with 2191909 peptidoforms detected in total across all datasets
# Total number of peptidoforms post filtering by dataset ID: 421918 (~20%)
# Total number of peptidoforms post subsequent filtering by number of calibrated errors: 73672 (~3.36% of all)



# Run a Gaussian mixed model - we expect 1-3 modes, but we might have an unexpected larger number.
# Therefore, we run the model with 1-3 components but if something weird pops up we might have to increase
# that number. However, larger number of components increases the risk of overfitting.
n_components_range = range(1, 3)
print(f"Fitting Gaussian mixed models with {n_components_range} components and selecting the best fit")

# we want to do the model fitting and selection for each peptidoform
# using the associated values as data to run the model on.

# first we will run the GMM for 1, 2, or 3 components, select the best model, then store that in a dictionary or of models
# Dictionary to store the best model and its BIC for each peptidoform
best_models = {}
# for every peptidoform ID that has passed filtering
for peptidoform_id in filtered_data['peptidoform_id']:
    # Get all associated data for that peptidoform ID
    current_row = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]

    # Extract the calibrated errors values
    errors_list = current_row['calibrated_error'].values[0]

    # Reshape into a two-dimensional numpy array with one column - we need the data in this format to run GMMs
    errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

    # Set the best BIC value to infinite - lower values are better so we can update this for each mdoel and keep the best
    best_bic = np.inf
    best_gmm = None

    # For every number of components in our range
    for n_components in n_components_range:

        #fit a GMM model on the numpy array of calibrated errors and record the associated BIC value
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(errors_reformatted)
        bic = gmm.bic(errors_reformatted)

        # if it is smaller than the current lowest BIC, this model is better, update the best BIC and GMM
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm

    # add the best model and associated BIC score to the dictionary of best models
    best_models[peptidoform_id] = {'model': best_gmm, 'bic': best_bic}
###############################################

# TODO: add GMM bins and save tata to csv for each.

###############################################

# Filtering data to include peptidoforms with potential sulfonation + 4 bins to each side to
# assess false positive rates?
# List of acceptable ranges
acceptable_ranges = [(-0.013, -0.007), (-0.014, -0.006), (-0.0125, -0.0075)]

# Initialize the final list of bins
all_bins = []

# Loop through each acceptable range and create bins
for acceptable_range in acceptable_ranges:
    range_width = acceptable_range[1] - acceptable_range[0]  # calculate width of the range

    # Create bins: 3 to the left and 6 to the right of each specified range
    bins = []
    for i in range(-3, 7):  # 3 to the left, 1 original, 6 to the right
        new_range_start = acceptable_range[0] + i * range_width
        new_range_end = acceptable_range[1] + i * range_width
        bins.append((new_range_start, new_range_end))

    # Add the bins for the current acceptable range to the final list of all bins
    all_bins.extend(bins)

# Print out all the defined bins for all acceptable ranges
print("Defined bins: ", all_bins)

# for every created bin get the peptidoform IDs that fall within
for bin_index, bin_range in enumerate(all_bins):
    sulfo_filtered_peptidoform_ids = set()  # Using a set to avoid duplicate peptidoform IDs

    # Filter data for each bin
    for peptidoform_id, model_info in best_models.items():
        model = model_info['model']
        means = model.means_.ravel()

        if model.n_components == 1 and bin_range[0] <= means[0] <= bin_range[1]:
            sulfo_filtered_peptidoform_ids.add(peptidoform_id)
        elif model.n_components in [2, 3] and any(bin_range[0] <= mean <= bin_range[1] for mean in means):
            sulfo_filtered_peptidoform_ids.add(peptidoform_id)

    # Filter the 'filtered_data' DataFrame to include only rows with peptidoform IDs in 'sulfo_filtered_peptidoform_ids'
    df_tosave = filtered_data[filtered_data['peptidoform_id'].isin(sulfo_filtered_peptidoform_ids)]

    # Check if the resulting subset is empty
    if not df_tosave.empty:
        # If not empty, proceed to save the subset to a CSV file
        output_file = f'../output/postGMM_fitting_binrange_{bin_range[0]}_{bin_range[1]}.csv'

        df_tosave.to_csv(output_file, index=False)
        print(f"Written data for bin {bin_range} to {output_file}")
    else:
        # If the subset is empty, print a message indicating no data was saved for this bin
        print(f"No data matched criteria for bin {bin_range}. No file written.")

# # sulfo_filtered_data = {}
# #
# # print(f"Filtering data to only include peptiodoforms with potential sulfonation...")
# # ## Select peptidoforms that exhibit a shift towards sulfonated.
# # for peptidoform_id, model_info in best_models.items():
# #     model = model_info['model']
# #     means = model.means_.ravel()
# #
# #     if model.n_components == 1 and acceptable_range[0] <= means[0] <= acceptable_range[1]:
# #         sulfo_filtered_data[peptidoform_id] = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].iloc[0]
# #     elif model.n_components in [2, 3] and any(acceptable_range[0] <= mean <= acceptable_range[1] for mean in means):
# #         sulfo_filtered_data[peptidoform_id] = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].iloc[0]
# #
# # # Now, sulfo_filtered_data contains the datasets that meet the GMM criteria
# #
# # print(f"Writing results to sulfo_filtered_data.csv in the output folder...")
# # # Save sulfo_filtered_data to CSV
# # with open('../output/sulfo_filtered_data_newIDs.csv', 'w', newline='') as csvfile:
# #     writer = csv.writer(csvfile)
# #     writer.writerow(['PeptidoformID', 'Calibrated_mass_error'])  # Writing header
# #     for peptidoform_id, errors in sulfo_filtered_data.items():
# #         writer.writerow([peptidoform_id, ', '.join(map(str, errors))])
#
# # # Randomly select 10 entries for histogram plotting just as a test
# # random_keys = random.sample(list(sulfo_filtered_data.keys()), min(10, len(sulfo_filtered_data)))
#
# print(f"Generating histograms for each peptidoform...")
#
# # Create and save histograms in a PDF file
# pdf = matplotlib.backends.backend_pdf.PdfPages(f"../output/sulfo_histograms_thresh{str(n_calerrors_threshold)}_newIDs.pdf")
# for peptidoform_id, errors in sulfo_filtered_data.items():
#     plt.figure()
#     plt.hist(errors, bins=30, alpha=0.7)
#     plt.title(f'Histogram for {peptidoform_id}')
#     plt.xlabel('Calibrated_mass_error')
#     plt.ylabel('Frequency')
#     pdf.savefig()  # saves the current figure into a pdf page
#     plt.close()
# pdf.close()
#
# print(f"Histograms saved in '../output/sulfo_histograms_newIDs.pdf'")
#
#
