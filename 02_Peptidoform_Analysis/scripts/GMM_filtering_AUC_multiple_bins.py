import numpy as np
import pandas as pd
import ast  # For safely evaluating strings as Python literals
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from scipy.integrate import quad
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns # for kernel density estimate histogram plots
import logging #save log instead of printing

# Configure logging
logging.basicConfig(filename='my_script_log.log', level=logging.INFO,
                    format='%(asctime)s:%(levelname)s:%(message)s')

###### functions #######
def gaussian_auc(mean, covar, weight, start, end):
    # Function to calculate the weighted AUC for a single Gaussian component
    # probability density function formula:
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

# sanity checks
logging.info(f"Total number of peptidoforms detected across all data: {len(set(aggregated_data['peptidoform_id']))}")
logging.info(f"Sanity check: Are number of rows and number unique peptidoform IDs the same?\nnrows: {len(aggregated_data)}")

# Convert 'dataset_ID', 'calibrated_error', and 'protein' from string to list using literal eval
aggregated_data['dataset_ID'] = aggregated_data['dataset_ID'].apply(ast.literal_eval)
aggregated_data['calibrated_error'] = aggregated_data['calibrated_error'].apply(ast.literal_eval)
aggregated_data['protein'] = aggregated_data['protein'].apply(ast.literal_eval)

# Handle 'alt_protein' clean up before converting to list: replace 'nan', 'nan, ', and 'nan,' with ''
strings_to_replace = ['nan, ', 'nan,', 'nan']  # quick and dirty fix
for string in strings_to_replace:
    aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(lambda x: x.replace(string, ""))
aggregated_data['alt_protein'] = aggregated_data['alt_protein'].apply(ast.literal_eval)
logging.info("Data read in complete.")

# Filter data based on number of dataset_IDs and number of calibrated_error values
n_datasets_threshold = 3
filtered_by_dataset_ID = aggregated_data[aggregated_data['dataset_ID'].apply(lambda x: len(set(x)) >= n_datasets_threshold)]
logging.info(f"Total number of peptidoforms post filtering by dataset ID: {len(set(filtered_by_dataset_ID['peptidoform_id']))}")

n_calerrors_threshold = 90
filtered_data = filtered_by_dataset_ID[filtered_by_dataset_ID['calibrated_error'].apply(lambda x: len(x) >= n_calerrors_threshold)]
logging.info(f"Total number of peptidoforms post subsequent filtering by number of calibrated errors: {len(set(filtered_data['peptidoform_id']))}")
logging.info(f"min_n_datasets: {n_datasets_threshold} \nmin_n_calerrors: {n_calerrors_threshold}")

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
    logging.info(f"Peptidoform:{peptidoform_id}")
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
        logging.info(f"BIC for {n_components} components: {bic}")  # Print the BIC for the current number of components
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
            logging.info(f"Updated best model to {n_components} components with BIC: {best_bic}")


# Step 2: AUC calculation and visualisation of model histograms
# Store AUC percentages and set threshold
peptidoform_aucs = {}
# auc_percent_threshold = 15
auc_percent_thresholds = [20, 25, 30, 45, 50]

# Color Universal Design (CUD) color palette
colors = ['#0072B2',  # Blue
          '#E69F00',  # Orange
          '#009E73',  # Green
          '#CC79A7']  # Pink not reLLY needed here as only 3 components max

# Define the range of our m/z bin of interest
# NB: here I have used the narrow bins from our previous analysis as they showed highest pY proportions in bin of interest
# mz_bin_start = -0.0125
# mz_bin_end = -0.0075

# List of acceptable ranges; can add or remove in future
acceptable_ranges = [(-0.0125, -0.0075)]

# Initialize the final list of bins
all_bins = []

# Loop through each acceptable range and create bins
for acceptable_range in acceptable_ranges:
    range_width = acceptable_range[1] - acceptable_range[0]  # calculate width of the range

    # Create bins: 3 to the left and 6 to the right of each specified range; make sure ranges are rounded to 4 decimals
    bins = []
    for i in range(-3, 7):  # 3 to the left, 1 original, 6 to the right
        new_range_start = round(acceptable_range[0] + i * range_width, 4)
        new_range_end = round(acceptable_range[1] + i * range_width, 4)
        bins.append((new_range_start, new_range_end))

    # Add the bins for the current acceptable range to the final list of all bins
    all_bins.extend(bins)

# now we need histograms for each bin so we create a pdf for each
pdf_pages_all_bins = {
    bin_range: PdfPages(f'../output/histograms_post_auc_filtering_bin_{bin_range[0]}_{bin_range[1]}.pdf') for
    bin_range in all_bins}
# also create a pdf to store all histograms of all peptidoforms - not for this script run
# pdf_pages_all = PdfPages('../output/gmm_models_histograms.pdf')

# Iterate through all the bins and apply AUC filtering
for bin_range in all_bins:
    for auc_percent_threshold in auc_percent_thresholds:
        filtered_peptidoform_ids = set()  # Using a set to avoid duplicate peptidoform IDs
        auc_percentages_bin = []  # Store AUC percentages for this bin


    # Iterate through the best models and apply AUC filtering based on the current bin
        for peptidoform_id, model_data in best_models.items():
            model = model_data['model']
            errors_list = filtered_data[filtered_data['peptidoform_id'] == peptidoform_id]['calibrated_error'].values[0]
            errors_reformatted = np.array(errors_list, dtype=float).reshape(-1, 1)

        # Calculate the AUC over the current bin for each component of the best model
            component_aucs = [gaussian_auc_for_bin(mean[0], covar[0][0], weight, bin_range)
                              for mean, covar, weight in zip(model.means_, model.covariances_, model.weights_)]

            # Calculate AUC as a percentage for each component and check if any meets the threshold
            component_aucs_percentage = [(auc / sum(model.weights_)) * 100 for auc in component_aucs]
            auc_meets_threshold = any(auc_percentage >= auc_percent_threshold for auc_percentage in component_aucs_percentage)

            # Check if the AUC meets the threshold
            if auc_meets_threshold:
                filtered_peptidoform_ids.add(peptidoform_id) # if so add the peptidoform to our filtered data

                # Start a new figure for each peptidoform
                plt.figure(figsize=(10, 6))

                # Plot the histogram of the data
                n, bins, patches = plt.hist(errors_reformatted, bins=30, density=True, alpha=0.6, color='g')

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
                    plt.text(mean + plt.xlim()[1] / 15, plt.ylim()[1] * 0.9, f'{mean:.4f}', color=current_color,
                             horizontalalignment='center')

                # add title and AUC to the figure
                plt.title(f'Best Fit for Peptidoform ID: {peptidoform_id} in Bin {bin_range}')
            # round AUCs first
                auc_percentages_to_plot = [round((auc / sum(model.weights_)) * 100, 2) for auc in component_aucs]
                plt.text(0.05, 0.95, f'AUCs: {auc_percentages_to_plot}%', transform=plt.gca().transAxes, fontsize=10,
                         verticalalignment='top')

                # Save the current figure to the PDF for the current bin
                pdf_pages_all_bins[bin_range].savefig()
                # pdf_pages_all.savefig()  # Save each histogram to the all-inclusive PDF, not in this run
                plt.close()

    # Print the number of peptidoforms with AUC >= auc_percent_threshold for the current bin
        num_filtered_peptidoforms_bin = len(filtered_peptidoform_ids)
        logging.info(f"Number of peptidoforms in bin {bin_range} with AUC >= {auc_percent_threshold}%: {num_filtered_peptidoforms_bin}")

        # Plot AUC percentages distribution for the current bin
        plt.figure(figsize=(10, 6))
        sns.histplot(auc_percentages_bin, kde=True, color='skyblue', bins=30)
        plt.title(f'Distribution of AUC Percentages for Bin {bin_range}')
        plt.xlabel('AUC Percentage')
        plt.ylabel('Frequency')
        plt.savefig(f'../output/auc_percentages_distribution_bin_{bin_range[0]}_{bin_range[1]}.png')
        plt.close()  # Close the plot to avoid memory issues


        # Filter the 'filtered_data' DataFrame to include only rows with peptidoform IDs in 'sulfo_filtered_peptidoform_ids'
        df_tosave = filtered_data[filtered_data['peptidoform_id'].isin(filtered_peptidoform_ids)]

        # Check if the resulting subset is empty
        if len(filtered_peptidoform_ids) > 0:
            # Filter the 'filtered_data' DataFrame
            df_tosave = filtered_data[filtered_data['peptidoform_id'].isin(filtered_peptidoform_ids)]

            # Save the filtered data to a CSV file
            output_file = f'../output/postGMM_fitting_threshold_{auc_percent_threshold}_binrange_{bin_range[0]}_{bin_range[1]}.csv'
            df_tosave.to_csv(output_file, index=False)
            logging.info(f"Data for bin {bin_range} with AUC >= {auc_percent_threshold}% written to {output_file}")
        else:
            # If the subset is empty, print a message indicating no data was saved for this bin
            logging.info(f"No data matched criteria for bin {bin_range} at AUC threshold of {auc_percent_threshold}%. No file written.")

# Close all PDFs
for pdf in pdf_pages_all_bins.values():
    pdf.close()

# pdf_pages_all.close() # we never opened this one
logging.info('PDF files have been saved in output folder.')
