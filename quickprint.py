acceptable_ranges = [(-0.013, -0.007)] #, (-0.014, -0.006), (-0.0125, -0.0075)]

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
def custom_round(number, precision):
    return round(number / precision) * precision

rounded_bins = [(custom_round(bin_value[0], 0.0004), custom_round(bin_value[1], 0.0004)) for bin_value in all_bins]
print("Defined bins: ", rounded_bins)
# Print out all the defined bins for all acceptable ranges
print("Defined bins: ", [round(bin_value) for bin_value in all_bins])
