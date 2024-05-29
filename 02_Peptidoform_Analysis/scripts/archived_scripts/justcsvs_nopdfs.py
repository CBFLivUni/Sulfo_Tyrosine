import numpy as np
import pandas as pd
import ast  # For safely evaluating strings as Python literals
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from scipy.integrate import quad
import logging  # save log instead of printing

# Configure logging
logging.basicConfig(filename='my_script_log.log', level=logging.INFO,
                    format='%(asctime)s:%(levelname)s:%(message)s')

###### functions #######
def gaussian_auc(mean, covar, weight, start, end):
    # Function to calculate the weighted AUC for a single Gaussian component
    pdf = lambda x: weight * norm.pdf(x, mean, np.sqrt(covar))
    auc, _ = quad(pdf, start, end)
    return auc

# Function to calculate the weighted AUC for a single Gaussian component for a specific bin
def gaussian_auc_for_bin(mean, covar, weight, bin_range):
    return gaussian_auc(mean, covar, weight, bin_range[0], bin_range[1])

# specify the CSV file path to peptidoform data
csv_file = '../data/clean_new_peptidoform_IDs_aggregated_by_peptidoform_ID.csv'

# read in the data as a pandas .csv since this is what it originated as
aggregated_data = pd.read_csv(csv_file)

# Convert 'dataset_ID', 'calibrated_error', and 'protein' from string to list using literal eval
aggregated_data['dataset_ID'] = aggregated_data['dataset_ID'].apply(ast.literal_eval)
aggregated_data['calibrated_error'] = aggregated_data['calibrated_error'].apply(ast.literal_eval)
aggregated_data['protein'] = aggregated_data['protein'].apply(ast.literal_eval)

# Handle 'alt_protein' clean up before converting to list: replace 'nan', 'nan, ', and 'nan,' with ''
strings_to_replace = ['nan, ', 'nan,', 'nan']
for string in strings_to_replace:
    aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(lambda x: x.replace(string, ""))
aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(ast.literal_eval)

# Filter data based on number of dataset_IDs and number of calibrated_error values
n_datasets_threshold = 3
filtered_by_dataset_ID = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= n_datasets_threshold)]

n_calerrors_threshold = 90
filtered_data = filtered_by_dataset_ID[filtered_by_dataset_ID['calibrated_error'].apply(lambda x: len(x) >= n_calerrors_threshold)]

# Initiate empty dictionary to store best models
best_models = {}

# we need to fit the models for each peptidoform
for peptidoform_id in filtered_data['peptidoform_id']:
    errors_list = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].values[0]
    errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

    best_bic = np.inf
    for n_components in range(1, 4):
        gmm = GaussianMixture(n_components=n_components, random_state=0).fit(errors_reformatted)
        bic = gmm.bic(errors_reformatted)
        if bic < best_bic:
            best_bic = bic
            best_models[peptidoform_id] = {'model': gmm, 'bic': bic}

# Define AUC percentage thresholds
auc_percent_thresholds = [20, 25, 30, 45, 50]

# Define the range of bins of interest
acceptable_ranges = [(-0.0125, -0.0075)]
all_bins = []
for acceptable_range in acceptable_ranges:
    range_width = acceptable_range[1] - acceptable_range[0]
    bins = [(round(acceptable_range[0] + i * range_width, 4), round(acceptable_range[1] + i * range_width, 4)) for i in range(-3, 7)]
    all_bins.extend(bins)

# Iterate through all the bins and apply AUC filtering for each threshold
for bin_range in all_bins:
    for auc_percent_threshold in auc_percent_thresholds:
        filtered_peptidoform_ids = set()
        for peptidoform_id, model_data in best_models.items():
            model = model_data['model']
            errors_list = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].values[0]
            errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

            component_aucs = [gaussian_auc_for_bin(mean[0], covar[0][0], weight, bin_range)
                              for mean, covar, weight in zip(model.means_, model.covariances_, model.weights_)]
            component_aucs_percentage = [(auc / sum(model.weights_)) * 100 for auc in component_aucs]
            if any(auc_percentage >= auc_percent_threshold for auc_percentage in component_aucs_percentage):
                filtered_peptidoform_ids.add(peptidoform_id)

        # Filter the 'filtered_data' DataFrame and save to CSV
        if filtered_peptidoform_ids:
            df_tosave = filtered_data[filtered_data['peptidoform_id'].isin(filtered_peptidoform_ids)]
            output_file = f'../output/postGMM_fitting_threshold_{auc_percent_threshold}_binrange_{bin_range[0]}_{bin_range[1]}.csv'
            df_tosave.to_csv(output_file, index=False)
            logging.info(f"Data for bin {bin_range} with AUC >= {auc_percent_threshold}% written to {output_file}")
        else:
            logging.info(f"No data matched criteria for bin {bin_range} at AUC threshold of {auc_percent_threshold}%. No file written.")

logging.info('CSV files have been saved in output folder.')