import pandas as pd
import numpy as np
import logging
import sys
import os

# Set up logging
logging.basicConfig(filename='edu_seq_analysis.log', level=logging.INFO)

# Input variables
BIN_SIZE = 10000
SCALE_FACTOR = 1000  # Scaling factor for sigma calculation

# Load the input file names from command-line arguments
adjusted_counts_file = sys.argv[1]  # EduHU_HCT_Biotin_set2A_adjusted_sample_counts.txt
bin_counts_file = sys.argv[2]       # EduHU_HCT_Biotin_set2A_sample_bin_counts.txt
totalsheared_file = sys.argv[3]     # EduHU_HCT_TotalSheared_set2A_adjust.csv

# Extract the basename from the adjusted counts input file (without extension)
sample_basename = os.path.basename(adjusted_counts_file).split('_adjusted_sample_counts.txt')[0]

# Dynamically generate output file names based on the sample name
output_file = f"{sample_basename}_sigma_output.csv"
qual_counts_output = f"{sample_basename}_sigma_qual_counts.txt"
select_output_file = f"{sample_basename}_sigma_select_EU_0b.csv"
smoothed_output_file = f"{sample_basename}_sigma_all_EU_0b.csv"
qual_counts_eu_output = f"{sample_basename}_sigma_qual_counts_EU_0b.txt"

# Step 1: Read input files
try:
    logging.info(f"Loading input files: {adjusted_counts_file}, {bin_counts_file}, {totalsheared_file}")
    
    # Reading space-delimited .txt files
    adjusted_counts = pd.read_csv(adjusted_counts_file, delim_whitespace=True, names=["chromosome", "bin", "adjusted_1", "adjusted_2"])
    bin_counts = pd.read_csv(bin_counts_file, delim_whitespace=True, names=["chromosome", "bin", "bin_count_1", "bin_count_2"])
    
    # Reading CSV file (TotalSheared data)
    totalsheared = pd.read_csv(totalsheared_file, names=["chromosome", "bin", "sheared_counts"])
    
except Exception as e:
    logging.error(f"Error loading files: {e}")
    sys.exit(1)

# Merge the data on chromosome and bin position
merged_data = pd.merge(adjusted_counts, bin_counts, on=["chromosome", "bin"])
merged_data = pd.merge(merged_data, totalsheared, on=["chromosome", "bin"])

# Initialize columns for sigma calculations
merged_data['sigma'] = 0

# Calculate sigma values based on the original Perl logic
for index, row in merged_data.iterrows():
    if row['sheared_counts'] > 0:  # Avoid division by zero
        # Sigma calculation (adjusted bin counts / total sheared counts)
        sigma = (row['bin_count_1'] / row['sheared_counts']) * SCALE_FACTOR
        merged_data.at[index, 'sigma'] = sigma

# Save the sigma output to a new CSV file
logging.info(f"Saving sigma calculations to {output_file}")
merged_data.to_csv(output_file, index=False)

# Step 2: Percentile-based background noise calculation
def find_background_noise(data, column, low_percentile=9, high_percentile=99):
    """Find background noise based on percentiles."""
    non_zero_values = data[column][data[column] > 0]
    background_low = np.percentile(non_zero_values, low_percentile)
    background_high = np.percentile(non_zero_values, high_percentile)
    return background_low, background_high

# Calculate background noise percentiles
background_low, background_high = find_background_noise(merged_data, 'adjusted_1')

# Step 3: Sigma calculation with background adjustment
def calculate_sigma_with_background(data, background_low, background_high):
    """Calculate sigma with background adjustment."""
    data['sigma_mb'] = (data['adjusted_1'] - background_low) / (background_high - background_low)
    return data

merged_data = calculate_sigma_with_background(merged_data, background_low, background_high)

# Step 4: Fit power curve (log-log linear regression) to the data
def fit_power_curve(x, y):
    """Fit a power curve to sigma values."""
    log_x = np.log(x)
    log_y = np.log(y)
    
    # Perform linear regression
    slope, intercept = np.polyfit(log_x, log_y, 1)
    return slope, intercept

adjust_read_values = merged_data[merged_data['sheared_counts'] != 0]['sheared_counts'].tolist()
adjust_sd_values = merged_data[merged_data['sheared_counts'] != 0]['sigma'].tolist()

slope, intercept = fit_power_curve(adjust_read_values, adjust_sd_values)

# Apply power curve to calculate sigma based on adjusted reads
merged_data['fitted_sigma'] = np.exp(intercept) * np.power(merged_data['sheared_counts'], slope)

# Step 5: Smooth and trim sigma values based on percentile background
def smooth_trim_sigma(data, window_size=3, trim_factor=1.2):
    """Smooth and trim sigma values."""
    data['smoothed_sigma'] = data['fitted_sigma'].rolling(window=window_size, center=True).mean()
    
    # Trimming step: Apply trimming based on trim_factor
    condition = data['fitted_sigma'] > (data['fitted_sigma'].shift(1) * trim_factor)
    data.loc[condition, 'trimmed_sigma'] = data['fitted_sigma'].shift(1) * trim_factor
    data['trimmed_sigma'].fillna(data['fitted_sigma'], inplace=True)
    return data

merged_data = smooth_trim_sigma(merged_data)

# Step 6: Convert sigma values to log2 scale and adjust the baseline
def log2_convert_sigma(data, baseline_mean):
    """Convert sigma values to log2 scale and adjust baseline."""
    data['sigma_log2'] = np.log2(data['trimmed_sigma'] + baseline_mean)
    return data

baseline_mean = merged_data['trimmed_sigma'].mean()
merged_data = log2_convert_sigma(merged_data, baseline_mean)

# Step 7: Save final output files
logging.info(f"Saving final sigma calculations to {select_output_file}")
merged_data.to_csv(select_output_file, index=False)

# Save smoothed and trimmed sigma results
logging.info(f"Saving smoothed sigma results to {smoothed_output_file}")
merged_data[['smoothed_sigma', 'trimmed_sigma']].to_csv(smoothed_output_file, index=False)

# Write the final quality counts and power curve details to a separate file
with open(qual_counts_eu_output, 'w') as qual_file:
    qual_file.write(f"File Name: {output_file}\n")
    qual_file.write(f"Bin size: {BIN_SIZE}\n")
    qual_file.write(f"Background Noise (Low): {background_low}\n")
    qual_file.write(f"Background Noise (High): {background_high}\n")
    qual_file.write(f"Slope of Power Curve: {slope}\n")
    qual_file.write(f"Intercept of Power Curve: {intercept}\n")
    qual_file.write(f"Baseline Mean (log2 adjusted): {baseline_mean}\n")

logging.info(f"Sigma Post Processing Complete.")
