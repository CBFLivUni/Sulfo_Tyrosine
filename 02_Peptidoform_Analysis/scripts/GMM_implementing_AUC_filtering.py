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
from scipy.stats import norm
from scipy.integrate import quad

# Specify the CSV file path
csv_file = '../data/clean_new_peptidoform_IDs_aggregated_by_peptidoform_ID.csv'

# read in the data as a pandas .csv since this is where it originated as
aggregated_data = pd.read_csv(csv_file, nrows= 50)
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
n_datasets_threshold = 1
filtered_by_dataset_ID = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= n_datasets_threshold)]
print(f"Total number of peptidoforms post filtering by dataset ID: {len(set(filtered_by_dataset_ID['peptidoform_id']))}")

n_calerrors_threshold = 10
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

    # Initialize the best BIC to infinity and the best model to None
    best_bic = np.inf
    best_gmm = None

    # Initialize the BIC of the simplest model (with the fewest components) to infinity
    simplest_bic = np.inf

    # For every number of components in our range
    for n_components in n_components_range:
        # Fit a GMM model on the numpy array of calibrated errors and record the associated BIC value
        gmm = GaussianMixture(n_components=n_components, random_state=0)
        gmm.fit(errors_reformatted)
        bic = gmm.bic(errors_reformatted)

        # Print the BIC for the current number of components
        print(f"BIC for {n_components} components: {bic}")

        # Update the simplest model's BIC if this is the first model (simplest one)
        if n_components == min(n_components_range):
            simplest_bic = bic

        # Update the best model if the BIC is lower than the best BIC by at least 10 points
        # or if this is the simplest model
        if bic < best_bic and (bic < simplest_bic - 10 or n_components == min(n_components_range)):
            best_bic = bic
            best_gmm = gmm
            print(f"Updated best model to {n_components} components with BIC: {best_bic}")

    # After the loop
    # best_gmm holds the best model
    # best_bic holds the BIC of the best model

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.stats import norm
    import numpy as np

    # Assuming best_models is your dictionary containing models and their BICs
    # and errors_reformatted is your data for the histogram

    # Prepare the PDF file
    pdf_pages = PdfPages('gmm_fit_results.pdf')

    # Plot for each peptidoform
    for peptidoform_id, model_data in best_models.items():
        model = model_data['model']

        # Start a new figure
        plt.figure(figsize=(10, 6))

        # Plot the histogram of the data
        n, bins, patches = plt.hist(errors_reformatted, bins=30, density=True, alpha=0.6, color='g')

        # Add a title and labels
        plt.title(f'Best Fit for Peptidoform ID: {peptidoform_id}')
        plt.xlabel('Calibrated Error')
        plt.ylabel('Density')

        # Plot each Gaussian component
        for i in range(model.n_components):
            mean = model.means_[i][0]
            var = model.covariances_[i][0]
            weight = model.weights_[i]

            # Create an array of values from the range of your data
            x = np.linspace(min(errors_reformatted), max(errors_reformatted), num=1000)

            # Compute the Gaussian function for these values
            gauss = weight * norm.pdf(x, mean, np.sqrt(var))

            # Plot the Gaussian function
            plt.plot(x, gauss, linewidth=2)

        # Combine the Gaussian components and plot the total density
        total_density = np.sum([weight * norm.pdf(bins, mean, np.sqrt(var)) for mean, var, weight in
                                zip(model.means_.flatten(), model.covariances_.flatten(), model.weights_)], axis=0)
        plt.plot(bins, total_density, color='red', linewidth=2, label='Total Density')

        # Add a legend
        plt.legend()

        # Save the current figure to the PDF
        pdf_pages.savefig()

    # Close the PDF and save the file
    pdf_pages.close()

    print('PDF file has been saved as gmm_fit_results.pdf')

    # add the best model and associated BIC score to the dictionary of best models
    best_models[peptidoform_id] = {'model': best_gmm, 'bic': best_bic}





    # Define the range of your m/z bin
    mz_bin_start = -0.0125
    mz_bin_end = -0.0075


    # Function to calculate the weighted AUC for a single Gaussian component
    def gaussian_auc(mean, covar, weight, start, end):
        # Define the probability density function (PDF) of the Gaussian component
        pdf = lambda x: weight * norm.pdf(x, mean, np.sqrt(covar))
        # Integrate the PDF over the bin range
        auc, _ = quad(pdf, start, end)
        return auc


    # Calculate the total AUC over the m/z bin for each best model
    for peptidoform_id, model_data in best_models.items():
        model = model_data['model']
        total_auc = sum(gaussian_auc(mean, covar, weight, mz_bin_start, mz_bin_end)
                        for mean, covar, weight in
                        zip(model.means_.flatten(), model.covariances_.flatten(), model.weights_))
        print(f"Peptidoform ID: {peptidoform_id}, AUC over m/z bin: {total_auc}")