#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import logging
import os

# Set up logging
logging.basicConfig(filename='edu_seq_analysis.log', level=logging.INFO)

# Input variables
BIN_SIZE = 10000
CORRECTION_FACTOR = 1.0  # Placeholder for correction factor, can be adjusted later
SCALE_FACTOR = 1000  # Scaling factor added to sigma calculation

# Load the input file names from command-line arguments
adjusted_counts_file = sys.argv[1]  # EduHU_HCT_Biotin_set2A_adjusted_sample_counts.txt
bin_counts_file = sys.argv[2]       # EduHU_HCT_Biotin_set2A_sample_bin_counts.txt
totalsheared_file = sys.argv[3]     # EduHU_HCT_TotalSheared_set2A_adjust.csv

# Extract sample name from the input file
sample_name = os.path.basename(adjusted_counts_file).split('_')[0]

# Define the output filenames with sample name
output_file = f"results/sigma/{sample_name}_sigma_select_EU_0b.csv"
qual_counts_output = f"results/sigma/{sample_name}_sigma_qual_counts_EU_0b.txt"
smoothed_output_file = f"results/sigma/{sample_name}_sigma_all_EU_0b.csv"

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

# Step 2: Merge the data on chromosome and bin position
merged_data = pd.merge(adjusted_counts, bin_counts, on=["chromosome", "bin"])
merged_data = pd.merge(merged_data, totalsheared, on=["chromosome", "bin"])

# Step 3: Calculate sigma values with scaling
merged_data['sigma'] = 0
for index, row in merged_data.iterrows():
    if row['sheared_counts'] > 0:  # Avoid division by zero
        # Sigma calculation with scaling by 1000
        sigma = (row['bin_count_1'] / row['sheared_counts']) * SCALE_FACTOR
        merged_data.at[index, 'sigma'] = sigma

# Save the intermediate sigma output
logging.info(f"Saving intermediate sigma calculations to {output_file}")
merged_data.to_csv(output_file, index=False)

# Step 4: Sort adjusted bin values for background noise calculation
all_non0_adjbin = merged_data[merged_data['sheared_counts'] != 0]['adjusted_1'].tolist()
sort_all_non0_adjbin = sorted(all_non0_adjbin)

# Step 5: Define background noise (percentile-based calculation)
def find_background_noise(sorted_adjbin, low_percentile=9, high_percentile=99):
    low_index = int(len(sorted_adjbin) * (low_percentile / 100))
    high_index = int(len(sorted_adjbin) * (high_percentile / 100))
    background_low = sorted_adjbin[low_index]
    background_high = sorted_adjbin[high_index]
    return background_low, background_high

background_low, background_high = find_background_noise(sort_all_non0_adjbin)

# Step 6: Calculate sigma with background adjustment
def calculate_sigma_with_background(data, background_low, background_high):
    data['sigma_mb'] = (data['adjusted_1'] - background_low) / (background_high - background_low)
    return data

merged_data = calculate_sigma_with_background(merged_data, background_low, background_high)

# Step 7: Fit a power curve (log-log linear regression) to the data
def fit_power_curve(x, y):
    log_x = np.log(x)
    log_y = np.log(y)
    
    # Perform linear regression
    slope, intercept = np.polyfit(log_x, log_y, 1)
    return slope, intercept

adjust_read_values = merged_data[merged_data['sheared_counts'] != 0]['sheared_counts'].tolist()
adjust_sd_values = merged_data[merged_data['sheared_counts'] != 0]['sigma'].tolist()

slope, intercept = fit_power_curve(adjust_read_values, adjust_sd_values)

# Apply the power curve to calculate sigma based on adjusted reads
merged_data['fitted_sigma'] = np.exp(intercept) * np.power(merged_data['sheared_counts'], slope)

# Step 8: Smooth and trim sigma values
def smooth_trim_sigma(data, window_size=3, trim_factor=1.2):
    data['smoothed_sigma'] = data['fitted_sigma'].rolling(window=window_size, center=True).mean()
    
    # Trimming step: Apply trimming based on trim_factor
    condition = data['fitted_sigma'] > (data['fitted_sigma'].shift(1) * trim_factor)
    data.loc[condition, 'trimmed_sigma'] = data['fitted_sigma'].shift(1) * trim_factor
    data['trimmed_sigma'].fillna(data['fitted_sigma'], inplace=True)
    return data

merged_data = smooth_trim_sigma(merged_data)

# Step 9: Convert sigma values to log2 scale and adjust the baseline
def log2_convert_sigma(data, baseline_mean):
    data['sigma_log2'] = np.log2(data['trimmed_sigma'] + baseline_mean)
    return data

baseline_mean = merged_data['trimmed_sigma'].mean()
merged_data = log2_convert_sigma(merged_data, baseline_mean)

# Step 10: Save final output files
logging.info(f"Saving final sigma calculations to {output_file}")
merged_data.to_csv(output_file, index=False)

# Save smoothed and trimmed sigma results
logging.info(f"Saving smoothed sigma results to {smoothed_output_file}")
merged_data[['smoothed_sigma', 'trimmed_sigma']].to_csv(smoothed_output_file, index=False)

# Write the quality counts and power curve details to a separate file
total_sample_hits = merged_data['bin_count_1'].sum()
total_adjust_hits = merged_data['sheared_counts'].sum()

with open(qual_counts_output, 'w') as qual_file:
    qual_file.write(f"File Name: {output_file}\n")
    qual_file.write(f"Bin size: {BIN_SIZE}\n")
    qual_file.write(f"Background Noise (Low): {background_low}\n")
    qual_file.write(f"Background Noise (High): {background_high}\n")
    qual_file.write(f"Slope of Power Curve: {slope}\n")
    qual_file.write(f"Intercept of Power Curve: {intercept}\n")
    qual_file.write(f"Baseline Mean (log2 adjusted): {baseline_mean}\n")
    qual_file.write(f"Total sample hits: {total_sample_hits}\n")
    qual_file.write(f"Total adjusted hits: {total_adjust_hits}\n")
    qual_file.write(f"Correction factor: {CORRECTION_FACTOR}\n")
    qual_file.write(f"Scale factor: {SCALE_FACTOR}\n")

logging.info(f"Combined sigma calculation script complete.")
