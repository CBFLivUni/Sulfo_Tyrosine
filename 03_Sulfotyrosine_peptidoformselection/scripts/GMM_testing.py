import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import csv
import matplotlib.pyplot as plt

# Initialize an empty dictionary to store the data
data_dict = {}

# Specify the CSV file path
csv_file = '../data/master_peptidoform_errors_06_12_23.csv'

# Open the CSV file and read the first 1000 lines
with open(csv_file, 'r') as file:
    reader = csv.reader(file)

    # Skip the header row - these are column names in our case
    next(reader, None)

    # Read and process only the first 1000 lines for testing purposes
    for row_num, row in enumerate(reader, start=1):
        # comment out to plot all data so far
        # if row_num > 1000:
        #     break  # Stop after the first 1000 lines

        # column 1 contains unique IDs and column 2 contains calibrated mass error values
        unique_id = row[0]
        value = row[1]

        # Check if the unique ID key already exists in the dictionary
        if unique_id in data_dict:
            # Append the value to the existing list
            data_dict[unique_id].append(value)
        else:
            # Create a new key with the value as a list
            data_dict[unique_id] = [value]
print(f"Total number of peptidofotms detected across all data: {len(data_dict)}")
# Now data_dict contains the data from the first 1000 lines of the CSV file
# As a sanity check we might want to plot a histogram of the number of values associated with each unique ID

# Extract the counts (number of values) for each key
# value_counts = [len(values) for values in data_dict.values()]

# Plot a histogram
# plt.hist(value_counts, bins=200, edgecolor='k', range=(1, 200))
# plt.xlabel('Number of Calibrated mass error Values')
# plt.ylabel('Frequency')
# plt.title('Histogram of Number of PSMs per Peptidoform')
# plt.grid(True)

# Save the plot as a PDF file - no need, better to show plot and be able to zoom in on it and save a few snapshots
# plt.savefig('../output/NumberOfPSMs_per_peptidoform_histogram_granular_.pdf', format='pdf')
#
# # Show the plot
# plt.show()
# we get lots of peptidoforms with only 1 PSM (~2.5*10^6) or < 10 PSMs,
# but quite a few with more than 10 out of just 31 datasets; so far the results are promising


# now lets filter the data so it includes only peptidoforms with sufficient number of PSMs to run a
# Gausian mixed model (GMM). e.g. lets start with total fo 10 PSMs minimum AND data comes from at leats 3 different datasets
# TODO: change the peptidoform ID collapse data funciton to include info on numbe rof datasets that the peptidoform was detected in

# Filter out keys with less than 50 values by keeping only keys with more than 50 values
# 50 was selected based on the following publication, which determines the number of values needed for
# different likelyhood accuracies of GMMs (Table 1 and 2 of interest):
# : https://www.sciencedirect.com/science/article/pii/S0031320319300445#:~:text=For%20instance%2C%20a%20sample%20size,only%202%20parameters%20are%20estimated).
filtered_dict = {key: values for key, values in data_dict.items() if len(values) >= 10}
print(f"Total number of peptidofotms remaining post-filtering: {len(filtered_dict)}")
# # Plot a histogram to see what we have got left
# # Extract the counts (number of values) for each key
# value_counts_postfilter = [len(values) for values in filtered_dict.values()]
# plt.hist(value_counts_postfilter, bins=200, edgecolor='k', range=(1, 200))
# plt.xlabel('Number of Calibrated mass error Values')
# plt.ylabel('Frequency')
# plt.title('Histogram of Number of PSMs per Peptidoform')
# plt.grid(True)
#
# # Save the plot as a PDF file - no need, better to show plot and be able to zoom in on it and save a few snapshots
# plt.savefig('../output/NumberOfPSMs_per_peptidoform_histogram_postfilter_threshold50.pdf', format='pdf')
#
# # Show the plot
# plt.show()


# Run a Gaussian mixed model - we expect 1-3 modes, but we might have an unexpected larger number.
# Therefore, we run the model with 1-3 components but if something weird pops up we might have to increase
# that number. However, larger number of components increases the risk of overfitting.
n_components_range = range(1, 3)
print("Fitting Gaussian mixed models with {n_components_range} components and selecting the best fit")

# we want to do this for each peptidoform
all_models = []


# we store the resulting models in a list
models = []
# and we keep track of their associated BIC scores for model selection purposes
bic_values = []

for n_components in n_components_range:
    gmm = GaussianMixture(n_components=n_components)
    gmm.fit(data)
    models.append(gmm)
    bic = gmm.bic(data)
    bic_values.append(bic)



filtered_data_dict = {}

# Define the range for centering check
# center_range = (0.02, 0.03)  #
#
# for dataset_name, dataset_info in data_dict.items():
#     values = dataset_info['values']
#     best_model = dataset_info['best_model']
#
#     # Check the number of components of the best model
#     num_components = best_model.n_components
#
#     if num_components == 1:  # Unimodal
#         # Check if the mean of the Gaussian component is within the specified range
#         mean = best_model.means_[0][0]
#         if center_range[0] <= mean <= center_range[1]:
#             filtered_data_dict[dataset_name] = values
#
#     elif num_components in (2, 3):  # Bimodal or trimodal
#         # Check if at least one cluster is centered around the specified range
#         means = best_model.means_.flatten()
#         if any(center_range[0] <= means <= center_range[1]):
#             filtered_data_dict[dataset_name] = values

# filtered_data_dict now contains the datasets that meet your criteria
