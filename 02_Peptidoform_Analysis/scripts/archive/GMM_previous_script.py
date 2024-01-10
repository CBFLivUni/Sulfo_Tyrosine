import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import csv
import sys, ctypes as ct
import matplotlib.pyplot as plt
import random
import matplotlib.backends.backend_pdf
import itertools


# Initialize an empty dictionary to store the data
calibrated_errors = {}

# Specify the CSV file path
csv_file = '../data/all_batches_peptidoforms_unique_aggregated_by_peptidoform_ID.csv'
# print(csv.field_size_limit())
# print(sys.maxsize)

# we have some fields of calibrated errors that exceed the limit
# increase the max field size limit of the csv reader, follwoing code taken from https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072
csv.field_size_limit(int(ct.c_ulong(-1).value // 2))
# print(csv.field_size_limit())


# Open the CSV file and read the first n lines
with open(csv_file, 'r') as file:
    reader = csv.reader(file)

    # Skip the header row - these are column names in our case
    next(reader, None)

    # Read and process only the first 1000 lines for testing purposes
    for row_num, row in enumerate(reader, start=1):

        # if row_num > n_lines:
        #     break  # Stop after the first n lines

        # column 1 with index 0 contains unique IDs and column 2 with index 1 contains calibrated mass error values
        unique_id = row[0]
        value = row[1]

        # Check if the unique ID key already exists in the dictionary
        if unique_id in calibrated_errors:
            # Append the value to the existing list
            calibrated_errors[unique_id].append(value)
        else:
            # Create a new key with the value as a list
            calibrated_errors[unique_id] = [value]
print(f"Total number of peptidofotms detected across all data: {len(calibrated_errors)}")


# Filter out keys with less than 50 values by keeping only keys with more than 50 values
# 50 was selected based on the following publication, which determines the number of values needed for
# different likelyhood accuracies of GMMs (Table 1 and 2 of interest):
# : https://www.sciencedirect.com/science/article/pii/S0031320319300445#:~:text=For%20instance%2C%20a%20sample%20size,only%202%20parameters%20are%20estimated).
filtered_errors = {key: values for key, values in calibrated_errors.items() if len(values) >= 2}
print(f"Total number of peptidofotms remaining post-filtering: {len(filtered_errors)}")




#
# # Run a Gaussian mixed model - we expect 1-3 modes, but we might have an unexpected larger number.
# # Therefore, we run the model with 1-3 components but if something weird pops up we might have to increase
# # that number. However, larger number of components increases the risk of overfitting.
# n_components_range = range(1, 3)
# print("Fitting Gaussian mixed models with {n_components_range} components and selecting the best fit")
#
# # we want to do the model fitting and selection for each peptidoform (key in filtered_dict)
# # using the associated values as data to run the model on.
#
# # first we will run the GMM for 1, 2, or 3 components, select the best model, then store that in a dictionary or list of models
# #
# # Convert string values to floats and reshape data for GMM
# for key in filtered_dict:
#     filtered_dict[key] = np.array(filtered_dict[key], dtype=float).reshape(-1, 1)
#
# # # Print the names of the first 10 keys in filtered_dict
# # first_10_keys = list(itertools.islice(filtered_dict.keys(), 10))
# # print(first_10_keys)
#
# # # Optionally, limit the dictionary to these 10 items for testing during further processing:
# # filtered_dict = {key: filtered_dict[key] for key in first_10_keys}
# # # print(filtered_dict)
#
#
# # Dictionary to store the best model and its BIC for each peptidofotm
# best_models = {}
#
# for key, values in filtered_dict.items():
#     best_bic = np.inf # low BIC is good so we start with an infinite number that can then be 'bested' when we run a model
#     best_gmm = None
#
#     for n_components in n_components_range:
#         # Fit the Gaussian Mixture Model
#         gmm = GaussianMixture(n_components=n_components, random_state=0)
#         gmm.fit(values)
#
#         # Calculate the BIC for the current model
#         bic = gmm.bic(values)
#
#         if bic < best_bic:  # update best BIC if new number of components results in a better model
#             best_bic = bic
#             best_gmm = gmm
#
#     # Store the best model and its BIC for the current peptidofotm
#     best_models[key] = {'model': best_gmm, 'bic': best_bic}
#
# # now, best_models contains the best GMM for each peptidofotm along with its BIC
# # print(best_models)
#
# print(f"Filtering data to only include peptiodoforms wiht potential sulfonation...")
# ## Select peptidoforms that exhibit a shift towards sulfonated.
#
# # Define the acceptable range for the mean - might need toying with.
# acceptable_range = (-0.015, -0.005)
#
# # Dictionary to store the final filtered data
# sulfo_filtered_data = {}
#
# for key, model_info in best_models.items():
#     model = model_info['model']
#     means = model.means_.ravel()  # Flatten the means array
#
#     # Check for unimodal distributions
#     if model.n_components == 1:
#         if acceptable_range[0] <= means[0] <= acceptable_range[1]:
#             sulfo_filtered_data[key] = filtered_dict[key]
#
#     # Check for bimodal or trimodal distributions
#     elif model.n_components in [2, 3]:
#         if any(acceptable_range[0] <= mean <= acceptable_range[1] for mean in means):
#             sulfo_filtered_data[key] = filtered_dict[key]
#
# # Now, sulfo_filtered_data contains the datasets that meet the criteria
#
# print(f"writing results to sulfo_filtered_data.csv in the output folder")
# # Save sulfo_filtered_data to CSV
# with open('../output/sulfo_filtered_data.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(['PeptidoformID', 'Calibrated_mass_error'])  # Writing header
#     for key, values in sulfo_filtered_data.items():
#         writer.writerow([key, ', '.join(map(str, values))])
#
# # # Randomly select 10 entries for histogram plotting just as a test
# # random_keys = random.sample(list(sulfo_filtered_data.keys()), min(10, len(sulfo_filtered_data)))
#
# print(f"Generating histograms for each peptidoform...")
#
# # Create and save histograms in a PDF file
# pdf = matplotlib.backends.backend_pdf.PdfPages("../output/sulfo_histograms.pdf")
# for key in sulfo_filtered_data:
#     plt.figure()
#     plt.hist(sulfo_filtered_data[key], bins=30, alpha=0.7)
#     plt.title(f'Histogram for {key}')
#     plt.xlabel('Value')
#     plt.ylabel('Frequency')
#     pdf.savefig()  # saves the current figure into a pdf page
#     plt.close()
# pdf.close()
#
# print(f"Histograms saved in '../output/sulfo_histograms.pdf'")