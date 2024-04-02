######### libraries and directories #########
library(tidyverse)
library(janitor)
library(viridis)
library(ggplot2)


# set wd to script folder
setwd("C:/Users/jtzve/Desktop/Sufo_Tyrosine/12_simple_instrument_proportions/scripts/")
# instrument_counts_by_peptidoform_data is in the data dir
data_dir <- "C:/Users/jtzve/Desktop/Sufo_Tyrosine/10_confident_subset_BOI/out/total_PSM_counts_by_instrument/"



####### read in data ##############
# peptidoforms of interest metadata - subset of proteins  that could be sulfated based on CC
BOI_all<- read.csv("../in/BOI_peptidoforms_annotated.csv")
BOI_all_peptidoforms <- BOI_all$peptidoform_id
# stricter subset where histograms were convincing
BOI_strict <- read.csv("../in/strict_potential_sY_peptidoforms_BOI.csv")
BOI_likely_sY_peptidoforms <- BOI_strict$peptidoform_id

datasets_metadata <- read.csv("../in/human_phosphobuild_metadata.csv") %>% clean_names()

########### get instrument contribution to build ############
# we are interested in how many PSMs we had available per instrument to start with
# for this we will aggregate by instrument and takethe number of total PSMs id'd for each

# head(datasets_metadata)

instrument_data <- datasets_metadata[, c("instrument_name", "spectra_id_d")]

instrument_data$instrument_name <- as.factor(instrument_data$instrument_name)

# clean spectra_id_d column
instrument_data <- instrument_data %>%
  mutate(spectra_id_d = gsub(",", "", spectra_id_d), # remove commas
         spectra_id_d = as.numeric(spectra_id_d)) # convert to numeric

# aggregate by instrument and sum the values in spectra_id_d
PSMs_by_instrument <- instrument_data %>%
  group_by(instrument_name) %>%
  summarise(total_PSMs = sum(spectra_id_d, na.rm = TRUE), .groups = 'drop')

PSMs_by_instrument$PSMs_proportion <- PSMs_by_instrument$total_PSMs/sum(PSMs_by_instrument$total_PSMs)

PSMs_by_instrument$percent_contribution <- PSMs_by_instrument$PSMs_proportion*100

PSMs_by_instrument %>%
  arrange(desc(total_PSMs)) %>%
  write.csv(file = "../out/PSMs_IDd_by_instrument_in_phosphobuild.csv", row.names = FALSE)


########### get instrument contribution to BOI ############

# we want a matrix with insturments in rows and peptidoform ID in columns 
# generate a df with instruments as rows
instrument_names <- gsub(" ", "_", PSMs_by_instrument$instrument_name)
instruments_contingency_data <- data.frame(instrument = instrument_names,
                                           PSM_totals_BOI = 0, # preallocate column for total counts
                                           PSM_proportions_BOI = 0,
                                           PSM_percentages_BOI = 0,
                                           row.names = instrument_names)
# get list of all counts data files
counts_data_files  <- list.files(data_dir, full.names = TRUE)
# test gsub
# cleaned_string <- gsub(".*total_counts_by_instrument_for_(.*)\\.csv$", "\\1", counts_data_files[2])
file  = counts_data_files[1]

# read in each file one at a time
for (file in counts_data_files) {
 
  # read in actual instrument counts
  temp_file <- read.csv(file)
  # this has instrument in column 1 and # PSMs contributed in column 2
  
  # get peptidoform ID
  current_peptidoform_ID <- names(temp_file)[2]
  # set column for peptidoform id to have 0 counts for all instruments
  instruments_contingency_data[, current_peptidoform_ID] <- 0
  
  instruments_contingency_data[temp_file$instrument, current_peptidoform_ID] <- temp_file[,2]
}

# calculate the totals for each row 
instruments_contingency_data <- instruments_contingency_data %>%
  mutate(PSM_totals_BOI = rowSums(select(., -c(instrument, PSM_totals_BOI, PSM_proportions_BOI,PSM_percentages_BOI)), na.rm = TRUE))

# conver to proportions for the whole bin
instruments_contingency_data$PSM_proportions_BOI <- instruments_contingency_data$PSM_totals_BOI/sum(instruments_contingency_data$PSM_totals_BOI)
instruments_contingency_data$PSM_percentages_BOI <- instruments_contingency_data$PSM_proportions_BOI*100


# calculate totals for each peptidoform
total_counts <- colSums(instruments_contingency_data[, 5:ncol(instruments_contingency_data)], na.rm = TRUE)
# add row to data frame
instruments_contingency_data <- rbind(instruments_contingency_data, 
                                      c("Totals", sum(total_counts), 1, 100, total_counts))

write.csv(instruments_contingency_data, "../out/instruments_counts_by_peptidoform.csv", row.names = FALSE)

############# for likely sY ###############
stringent_contingency_data <- data.frame(instrument = instrument_names,
                                                                 PSM_totals_BOI = 0, # preallocate column for total counts
                                                                 PSM_proportions_BOI = 0,
                                                                 PSM_percentages_BOI = 0,
                                                                 row.names = instrument_names)

stringent_subset <- instruments_contingency_data[, c("instrument", BOI_strict$peptidoform_id)]
# str(stringent_subset)
# bind by common instrument name; in both dfs the instrument column should match,remainign cols should be kept
stringent_contingency_data <- left_join(stringent_contingency_data,
                                        stringent_subset,
                                        by = "instrument") 

# convertt all columns but instrument column to numeric 
stringent_contingency_data <- stringent_contingency_data %>%
  mutate(across(-instrument, ~suppressWarnings(as.numeric(.))))

# str(stringent_contingency_data)
# add totals and proportions as before
stringent_contingency_data <- stringent_contingency_data %>%
  mutate(PSM_totals_BOI = rowSums(select(., -c(instrument, PSM_totals_BOI, PSM_proportions_BOI ,PSM_percentages_BOI)), na.rm = TRUE))

# conver to proportions for the whole bin
stringent_contingency_data$PSM_proportions_BOI <- stringent_contingency_data$PSM_totals_BOI/sum(stringent_contingency_data$PSM_totals_BOI)
stringent_contingency_data$PSM_percentages_BOI <- stringent_contingency_data$PSM_proportions_BOI*100


# calculate totals for each peptidoform
total_counts <- colSums(stringent_contingency_data[, 5:ncol(stringent_contingency_data)], na.rm = TRUE)
# add row to data frame
stringent_contingency_data <- rbind(stringent_contingency_data, 
                                      c("Totals", sum(total_counts), 1, 100, total_counts))

write.csv(stringent_contingency_data, "../out/stringent_instruments_counts_by_peptidoform.csv", row.names = FALSE)

###########  for known sY ##############
# get list of detected peptidoforms
knonw_sYs <- read.csv("../in/knownsY_peptidoforms_annotated.csv")$peptidoform_id
# NB: we have detected Vitronextin and Secretogranin-1 as known sY

#  subset columns to only oinclude peptidoforms form them
known_sY_contingency_data <- data.frame(instrument = instrument_names,
                                      PSM_totals_BOI = 0, # preallocate column for total counts
                                      PSM_proportions_BOI = 0,
                                      PSM_percentages_BOI = 0,
                                      row.names = instrument_names)

#  generate contingency table
known_subset <- instruments_contingency_data[, c("instrument", knonw_sYs)]
# str(stringent_subset)
# bind by common instrument name; in both dfs the instrument column should match,remainign cols should be kept
known_sY_contingency_data <- left_join(known_sY_contingency_data,
                                       known_subset,
                                        by = "instrument") 

# convertt all columns but instrument column to numeric 
known_sY_contingency_data <- known_sY_contingency_data %>%
  mutate(across(-instrument, ~suppressWarnings(as.numeric(.))))

# str(stringent_contingency_data)
# add totals and proportions as before
known_sY_contingency_data <- known_sY_contingency_data %>%
  mutate(PSM_totals_BOI = rowSums(select(., -c(instrument, PSM_totals_BOI, PSM_proportions_BOI ,PSM_percentages_BOI)), na.rm = TRUE))

# conver to proportions for the whole bin
known_sY_contingency_data$PSM_proportions_BOI <- known_sY_contingency_data$PSM_totals_BOI/sum(known_sY_contingency_data$PSM_totals_BOI)
known_sY_contingency_data$PSM_percentages_BOI <- known_sY_contingency_data$PSM_proportions_BOI*100


# calculate totals for each peptidoform
total_counts <- colSums(known_sY_contingency_data[, 5:ncol(known_sY_contingency_data)], na.rm = TRUE)
# add row to data frame
known_sY_contingency_data <- rbind(known_sY_contingency_data, 
                                    c("Totals", sum(total_counts), 1, 100, total_counts))

write.csv(known_sY_contingency_data, "../out/known_sY_instruments_counts_by_peptidoform.csv", row.names = FALSE)

# generate baloon plots
colnames(PSMs_by_instrument)[c(1, 4) ] <- c("instrument", "PSM_percentages_all_data")
PSMs_by_instrument$instrument <- gsub (" ", "_", PSMs_by_instrument$instrument)
PSMs_by_instrument$PSM_percentages_all_data <- as.character(PSMs_by_instrument$PSM_percentages_all_data)

columns_to_include <- c("instrument", "PSM_percentages_BOI")

baloon_plot_percentages_data <- PSMs_by_instrument %>%
  left_join(instruments_contingency_data[-16, columns_to_include], by = "instrument", suffix = c("", "_all")) %>%
  left_join(stringent_contingency_data[-16, columns_to_include], by = "instrument", suffix = c("", "_stringent")) %>%
  left_join(known_sY_contingency_data[-16, columns_to_include], by = "instrument", suffix = c("", "_knownsY"))

# retain only columns to plot and conver to long format
baloon_plot_percentages_data_to_plot <- baloon_plot_percentages_data [ , c(1, 4:7)] %>%
  pivot_longer(cols = starts_with("PSM_percentages"),
               names_to = "category",
               values_to = "percentage") %>%
  mutate(percentage = as.numeric(percentage))

# plot baloon plot
# sort instruments descending by $percentages_PSM_all_data
PSMs_by_instrument$PSM_percentages_all_data <- as.numeric(PSMs_by_instrument$PSM_percentages_all_data)
instrument_plot_order <- PSMs_by_instrument %>%
  arrange(desc(PSM_percentages_all_data)) %>%
  pull(instrument)

category_order <- c("PSM_percentages_BOI_knownsY", "PSM_percentages_BOI_stringent","PSM_percentages_BOI", "PSM_percentages_all_data")
# apply the sorted order to the 'instrument' column for plotting
baloon_plot_percentages_data_to_plot$instrument <- factor(baloon_plot_percentages_data_to_plot$instrument, levels = instrument_plot_order)
baloon_plot_percentages_data_to_plot$category <- factor(baloon_plot_percentages_data_to_plot$category, levels = category_order)
baloon_plot_percentages_data_to_plot$percentage <- as.numeric(baloon_plot_percentages_data_to_plot$percentage)
# plot in that order form left to right

write.csv(baloon_plot_percentages_data, file = "../out/baloon_plot_wide_data.csv", row.names = FALSE)


baloon_plot <- ggplot(baloon_plot_percentages_data_to_plot, aes(x = instrument, y = category, size = percentage, color = category)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(1, 20)) +
  theme_minimal()+
  labs(title = "Percentage of PSMs contributed to a category by Instrument",
       x = "Instrument",
       y = "Category",
       size = "Percent Contribution"
       ) +
  guides(color = "none", # this removes the color legend
         size = guide_legend(title = "Percent Contribution")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  geom_text(aes(label = ifelse(percentage >= 20, paste0(sprintf("%.2f", percentage), "%"), ""), y = category),
            position = position_dodge(width = 0), vjust = 0.5, colour = "black",
            size = 3, check_overlap = FALSE)

ggsave(plot = baloon_plot, 
       filename = "../out/PSMs_Percent_contribution_by_category_baloonplot.png", 
       dpi = 600, width = 10, height = 8)

# generate bar charts

