import numpy as np
import pandas as pd
import ast  # For safely evaluating strings as Python literals
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns # for kernel density estimate histogram plots

# specify the CSV file path to peptidoform data
csv_file = '../data/clean_new_peptidoform_IDs_aggregated_by_peptidoform_ID.csv'

# read in the data as a pandas .csv since this is what it originated as
aggregated_data = pd.read_csv(csv_file, nrows=10000)

# sanity checks
print(f"Total number of peptidoforms detected across all data: {len(set(aggregated_data['peptidoform_id']))}")
print(f"Sanity check: Are number of rows and number unique peptidoform IDs the same?\nnrows: {len(aggregated_data)}")

# Convert 'dataset_ID', 'calibrated_error', and 'protein' from string to list using literal eval
aggregated_data['dataset_ID'] = aggregated_data['dataset_ID'].apply(ast.literal_eval)
aggregated_data['calibrated_error'] = aggregated_data['calibrated_error'].apply(ast.literal_eval)
aggregated_data['protein'] = aggregated_data['protein'].apply(ast.literal_eval)

# Handle 'alt_protein' clean up before converting to list: replace 'nan', 'nan, ', and 'nan,' with ''
strings_to_replace = ['nan, ', 'nan,', 'nan']  # quick and dirty fix
for string in strings_to_replace:
    aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(lambda x: x.replace(string, ""))
aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(ast.literal_eval)
print("Data read in complete.")

# Filter data based on number of dataset_IDs and number of calibrated_error values
n_datasets_threshold = 3
filtered_by_dataset_ID = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= n_datasets_threshold)]
print(f"Total number of peptidoforms post filtering by dataset ID: {len(set(filtered_by_dataset_ID['peptidoform_id']))}")

n_calerrors_threshold = 90
filtered_data = filtered_by_dataset_ID[filtered_by_dataset_ID['calibrated_error'].apply(lambda x: len(x) >= n_calerrors_threshold)]
print(f"Total number of peptidoforms post subsequent filtering by number of calibrated errors: {len(set(filtered_data['peptidoform_id']))}")
print(f"min_n_datasets: {n_datasets_threshold} \nmin_n_calerrors: {n_calerrors_threshold}")

# Selection of peptidoforms that fall in our bin of interest

# Step1: Fit gaussian mixture models with 1-4 components and select the optimal number of components for each peptidoform
n_components_range = range(1, 4)

# reason why we look for two componenents - we expect some sulfation, so centering around -0.01 and 0.
# reason for 3 - we have observed some data at + 0.01 in initial histograms

# select which model fits the data best depending on nuber of components, and BIC score
# we want the lowest BIC score but also we don't want to overfit - we have seen sometimes
# multiple components are selected based on BIC alone where eyeballing the histogram
# would have resulted in a lower number of components. Therefore, we prioritise a low score
# only if it's 'significantly' lower than the score of the smallest number of components with the best socre.
# rule of thumb - 10 points lower score is significant, but also test by plotting AUC

# initiate empty dictionary to store best models
best_models = {}

# we need to fit the models for each peptidoform
for peptidoform_id in filtered_data['peptidoform_id']:
    # print(f"Peptidoform:{peptidoform_id}")
    errors_list = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].values[0]
    # reshape data as required for GMM model
    errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

    # BIC scores are better when smaller.  initiate best score variable as infinitely large and replace with best score
    best_bic = np.inf
    # also track score for simplest model - e.g. 1 component in this case
    simplest_bic = np.inf

    for n_components in n_components_range:
        #fit model for every n of components
        gmm = GaussianMixture(n_components=n_components, random_state=0).fit(errors_reformatted)
        #get BIC score
        bic = gmm.bic(errors_reformatted)
        # print(f"BIC for {n_components} components: {bic}")  # Print the BIC for the current number of components
        # track bic score for simplest model
        if n_components == min(n_components_range):
            simplest_bic = bic

            # if new bic score is better AND 10 smaller than the simplest BIC or is the simplest bic
            # 10 is rule of thumb, might not really matter much for very small BIC scores anyway e.g. -900
        if bic < simplest_bic - 10 or n_components == min(n_components_range):
            # update the best bic
            best_bic = bic
            simplest_bic = bic #technically this is the simplest best bic now, as it might represent 2 components
            best_gmm = gmm
            best_models[peptidoform_id] = {'model': best_gmm, 'bic': best_bic}
            # print(f"Updated best model to {n_components} components with BIC: {best_bic}")

############### so far everything works as intended ###############
# now we wanto to visualise the data - e.g. histograms with best fit line
# NB: we get the histograms here during testing only! otherwise we save  them as pdfs for each bin?

# we want to also calculate the area under the curve for different m/z bins
# lets do it at this point and print it in the title of the plot to make sure everything works ok

# Step 2: AUC calculation and visualisation of model histograms
# Store AUC percentages and set threshold
peptidoform_aucs = {}
auc_percent_threshold = 15

# Prepare two PDFs to save the plots - one for all plots, one for just the peptidoforms that pass AUC filtering
pdf_pages_all = PdfPages('../output/gmm_models_histograms.pdf')
pdf_pages_filtered = PdfPages('../output/histograms_post_auc_filtering_20pc_mz15_05.pdf')

# Define the range of our m/z bin of interest
# NB: here I have used the narrow bins from our previous analysis as they showed highest pY proportions in bin of interest
mz_bin_start = -0.0125
mz_bin_end = -0.0075



# Color Universal Design (CUD) colors
colors = ['#0072B2',  # Blue
          '#E69F00',  # Orange
          '#009E73',  # Green
          '#CC79A7']  # Pink not reLLY needed here as only 3 components max

# Iterate through the best models and plot the histograms with best fit lines
for peptidoform_id, model_data in best_models.items():
    model = model_data['model']
    errors_list = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].values[0]
    errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

    # Start a new figure
    plt.figure(figsize=(10, 6))

    # Plot the histogram of the data
    n, bins, patches = plt.hist(errors_reformatted, bins=30, density=True, alpha=0.6, color='g')

    # Add a title and labels
    plt.title(f'Best Fit for Peptidoform ID: {peptidoform_id}')
    plt.xlabel('Calibrated Error')
    plt.ylabel('Density')

    # Plot each Gaussian component, its median line, and text for the median
    for i in range(model.n_components):
        mean = model.means_[i][0]
        var = model.covariances_[i][0]
        weight = model.weights_[i]
        x = np.linspace(min(errors_reformatted), max(errors_reformatted), 1000)
        gauss = weight * norm.pdf(x, mean, np.sqrt(var))
        current_color = colors[i % len(colors)]  # Use modulo to cycle through colors
        plt.plot(x, gauss, linewidth=2, color=current_color)

        # Add a dashed line at the median (mean for Gaussian), using the same color
        plt.axvline(x=mean, color=current_color, linestyle='--')

        # Add median text annotation
        plt.text(mean + plt.xlim()[1]/15, plt.ylim()[1] * 0.9, f'{mean:.4f}', color=current_color, horizontalalignment='center')

    # Calculate the total AUC over the m/z bin for the best model
    total_auc = sum(gaussian_auc(mean, covar, weight, mz_bin_start, mz_bin_end)
                    for mean, covar, weight in
                    zip(model.means_.flatten(), model.covariances_.flatten(), model.weights_))

    # Calculate AUC as a percentage and round to 2 decimals
    auc_percentage = round((total_auc / sum(weight for weight in model.weights_)) * 100, 2)
    # Store the AUC percentage for each peptidoform
    peptidoform_aucs[peptidoform_id] = auc_percentage

    # Display AUC value on the plot
    plt.text(0.05, 0.95, f'AUC: {auc_percentage}%', transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')

    # print(f"Peptidoform ID: {peptidoform_id}, AUC over m/z bin: {auc_percentage}%")

    # Save the current figure to the PDF
    pdf_pages_all.savefig()


    # further filter data based on the computed AUC for the bin of interest
    filtered_peptidoforms_by_auc = {peptidoform_id: auc for peptidoform_id, auc in peptidoform_aucs.items() if
                                    auc >= auc_percent_threshold}
    # Check if the AUC meets the threshold and save to the filtered PDF if it does
    if auc_percentage >= auc_percent_threshold:
        pdf_pages_filtered.savefig()
    # Close the current figure plot to avoid memory issues, otherwise all plots are stored in memory
    plt.close()

# Close the PDFs and save the files
pdf_pages_all.close()
pdf_pages_filtered.close()
print('PDF files have been saved in output folder.')

# Print the number of peptidoforms with AUC >= auc_percent_threshold
num_filtered_peptidoforms = len(filtered_peptidoforms_by_auc)
print(f"Number of peptidoforms with AUC >= {auc_percent_threshold}%: {num_filtered_peptidoforms}")

# plot auc percentages across all data
# gwt list of percentages
auc_percentages = list(peptidoform_aucs.values())

# Start a new figure
plt.figure(figsize=(10, 6))

# Create a histogram of the AUC percentages with kernel density estimate
sns.histplot(auc_percentages, kde=True, color='skyblue', bins=30)

# Add a title and labels
plt.title('Distribution of AUC Percentages')
plt.xlabel('AUC Percentage')
plt.ylabel('Frequency')

# save this figure to a file
plt.savefig(f'../output/auc_percentages_distribution_mz_{mz_bin_start}_{mz_bin_end}.png')

