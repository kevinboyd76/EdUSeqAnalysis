#!/usr/bin/env python

import pandas as pd
import numpy as np
import math

# Load the intermediate sigma output file (from the first script)
sigma_data = pd.read_csv("sigma_output.csv")

# Constants
number_split_adj_read = 200  # How many bins we will split the adjusted reads into
bin_size = 10000  # Assuming the bin size used earlier

# Step 1: Sort adjusted bin values
all_non0_adjbin = sigma_data[sigma_data['sheared_counts'] != 0]['adjusted_1'].tolist()
sort_all_non0_adjbin = sorted(all_non0_adjbin)

# Step 2: Define background noise
# Example: Find background noise level by analyzing adjbin percentiles
def find_background_noise(sorted_adjbin, low_percentile=9, high_percentile=99):
    low_index = int(len(sorted_adjbin) * (low_percentile / 100))
    high_index = int(len(sorted_adjbin) * (high_percentile / 100))
    background_low = sorted_adjbin[low_index]
    background_high = sorted_adjbin[high_index]
    return background_low, background_high

background_low, background_high = find_background_noise(sort_all_non0_adjbin)

# Step 3: Calculate SD for bins based on background noise
def calculate_sigma_with_background(data, background_low, background_high):
    data['sigma_mb'] = (data['adjusted_1'] - background_low) / (background_high - background_low)
    return data

sigma_data = calculate_sigma_with_background(sigma_data, background_low, background_high)

# Step 4: Power curve fitting (adjust reads vs adjbin_SD)
# Assume y = a * x^b and convert it to log-log linear form: log(y) = log(a) + b*log(x)
def fit_power_curve(x, y):
    log_x = np.log(x)
    log_y = np.log(y)
    
    # Perform linear regression
    slope, intercept = np.polyfit(log_x, log_y, 1)
    return slope, intercept

adjust_read_values = sigma_data[sigma_data['sheared_counts'] != 0]['sheared_counts'].tolist()
adjust_sd_values = sigma_data[sigma_data['sheared_counts'] != 0]['sigma'].tolist()

slope, intercept = fit_power_curve(adjust_read_values, adjust_sd_values)

# Apply power curve to calculate sigma based on adjusted reads
sigma_data['fitted_sigma'] = np.exp(intercept) * np.power(sigma_data['sheared_counts'], slope)

# Step 5: Smooth and Trim Sigma values
def smooth_trim_sigma(data, window_size=3, trim_factor=1.2):
    data['smoothed_sigma'] = data['fitted_sigma'].rolling(window=window_size, center=True).mean()
    
    # Trimming step: Apply trimming based on trim_factor
    condition = data['fitted_sigma'] > (data['fitted_sigma'].shift(1) * trim_factor)
    data.loc[condition, 'trimmed_sigma'] = data['fitted_sigma'].shift(1) * trim_factor
    data['trimmed_sigma'].fillna(data['fitted_sigma'], inplace=True)
    return data

sigma_data = smooth_trim_sigma(sigma_data)

# Step 6: Convert sigma values to log2 scale and adjust baseline
def log2_convert_sigma(data, baseline_mean):
    data['sigma_log2'] = np.log2(data['trimmed_sigma'] + baseline_mean)
    return data

baseline_mean = sigma_data['trimmed_sigma'].mean()
sigma_data = log2_convert_sigma(sigma_data, baseline_mean)

# Step 7: Save final output files
sigma_data.to_csv("sigma_select_EU_0b.csv", index=False)

# Optional: You can also write additional files with sigma smoothing and trimming results
sigma_data[['smoothed_sigma', 'trimmed_sigma']].to_csv("sigma_all_EU_0b.csv", index=False)

# Also write a quality control file
with open("sigma_qual_counts_EU_0b.txt", "w") as qual_file:
    qual_file.write(f"File Name: sigma_output.csv\n")
    qual_file.write(f"Bin size: {bin_size}\n")
    qual_file.write(f"Background Noise (Low): {background_low}\n")
    qual_file.write(f"Background Noise (High): {background_high}\n")
    qual_file.write(f"Slope of Power Curve: {slope}\n")
    qual_file.write(f"Intercept of Power Curve: {intercept}\n")
    qual_file.write(f"Baseline Mean (log2 adjusted): {baseline_mean}\n")

