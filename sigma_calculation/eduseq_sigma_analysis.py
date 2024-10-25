import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import sys
import os

# Set up logging
logging.basicConfig(filename='eduseq_sigma_analysis.log', level=logging.INFO)

# Constants
BIN_SIZE = 10000
SCALE_FACTOR = 1000  # Scaling factor added to sigma calculation
MIN_VALUE = 1e-9  # Small constant to avoid log issues

# Input variables from command-line arguments
adjusted_counts_file = sys.argv[1]  # EduHU_HCT_Biotin_set2A_adjusted_sample_counts.txt
bin_counts_file = sys.argv[2]       # EduHU_HCT_Biotin_set2A_sample_bin_counts.txt
totalsheared_file = sys.argv[3]     # EduHU_HCT_TotalSheared_set2A_adjust.csv
work_dir = sys.argv[4]              # Working directory for outputs
manual_max = float(sys.argv[5]) if len(sys.argv) > 5 else None  # Optional manual y-axis maximum

# Ensure the work directory exists
os.makedirs(work_dir, exist_ok=True)

# Extract the basename from the adjusted counts input file (without extension)
sample_basename = os.path.basename(adjusted_counts_file).split('_adjusted_sample_counts.txt')[0]

# Dynamically generate output file names based on the sample name
output_file = os.path.join(work_dir, f"{sample_basename}_sigma_output.csv")
qual_counts_output = os.path.join(work_dir, f"{sample_basename}_sigma_qual_counts.txt")
select_output_file = os.path.join(work_dir, f"{sample_basename}_sigma_select_EU_0b.csv")
smoothed_output_file = os.path.join(work_dir, f"{sample_basename}_sigma_all_EU_0b.csv")
qual_counts_eu_output = os.path.join(work_dir, f"{sample_basename}_sigma_qual_counts_EU_0b.txt")
output_plot_file = os.path.join(work_dir, f"{sample_basename}_sigma_plot.png")

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

# Calculate total reads in the sample and control to compute the correction factor
total_sample_reads = merged_data['bin_count_1'].sum()
total_control_reads = merged_data['sheared_counts'].sum()

# Step 2: Calculate or manually input the correction factor
if len(sys.argv) > 5:  # If the user provides a manual correction factor
    correction_factor = float(sys.argv[5])
else:
    if total_control_reads > 0:
        correction_factor = total_sample_reads / total_control_reads
    else:
        logging.error(f"Total control reads is zero, unable to calculate correction factor.")
        sys.exit(1)

logging.info(f"Correction factor calculated: {correction_factor}")

# Initialize columns for sigma calculations
merged_data['sigma'] = 0

# Step 3: Calculate sigma values, applying the correction factor
for index, row in merged_data.iterrows():
    if row['sheared_counts'] > 0:  # Avoid division by zero
        sigma = (row['bin_count_1'] / row['sheared_counts']) * SCALE_FACTOR * correction_factor
        merged_data.at[index, 'sigma'] = sigma

# Save the sigma output to a new CSV file
logging.info(f"Saving sigma calculations to {output_file}")
merged_data.to_csv(output_file, index=False)

# Write the quality counts to a separate file
with open(qual_counts_output, 'w') as qc_file:
    qc_file.write(f"Bin size: {BIN_SIZE}\n")
    qc_file.write(f"Total sample hits: {total_sample_reads}\n")
    qc_file.write(f"Total adjusted hits: {total_control_reads}\n")
    qc_file.write(f"Correction factor: {correction_factor}\n")

# Step 4: Percentile-based smoothing (background noise)
all_non0_adjbin = merged_data[merged_data['sheared_counts'] != 0]['adjusted_1'].tolist()
sort_all_non0_adjbin = sorted(all_non0_adjbin)

# Percentile background noise calculation
def find_background_noise(sorted_adjbin, low_percentile=9, high_percentile=99):
    low_index = int(len(sorted_adjbin) * (low_percentile / 100))
    high_index = int(len(sorted_adjbin) * (high_percentile / 100))
    background_low = sorted_adjbin[low_index]
    background_high = sorted_adjbin[high_index]
    return background_low, background_high

background_low, background_high = find_background_noise(sort_all_non0_adjbin)

# Step 5: Adjust sigma based on background noise
def calculate_sigma_with_background(data, background_low, background_high):
    data['sigma_mb'] = (data['adjusted_1'] - background_low) / (background_high - background_low)
    return data

merged_data = calculate_sigma_with_background(merged_data, background_low, background_high)

# Step 6: Smooth and trim sigma values
def smooth_trim_sigma(data, window_size=3, trim_factor=1.2):
    data['smoothed_sigma'] = data['sigma_mb'].rolling(window=window_size, center=True).mean()
    
    # Trimming step: Apply trimming based on trim_factor
    condition = data['sigma_mb'] > (data['sigma_mb'].shift(1) * trim_factor)
    data.loc[condition, 'trimmed_sigma'] = data['sigma_mb'].shift(1) * trim_factor
    data['trimmed_sigma'].fillna(data['sigma_mb'], inplace=True)
    return data

merged_data = smooth_trim_sigma(merged_data)

# Step 7: Convert sigma values to log2 scale and adjust the baseline
def log2_convert_sigma(data, baseline_mean, min_value=1e-9):
    # Add a small constant to avoid log2(0) or log2(negative numbers)
    data['sigma_log2'] = np.log2(data['trimmed_sigma'].clip(lower=min_value) + baseline_mean)
    
    # Log any remaining issues
    invalid_entries = data[data['sigma_log2'].isnull()]
    if not invalid_entries.empty:
        logging.warning(f"Invalid log2 values encountered for the following entries: {invalid_entries}")
    
    return data

baseline_mean = merged_data['trimmed_sigma'].mean()
merged_data = log2_convert_sigma(merged_data, baseline_mean)

# Step 8: Save final output files
logging.info(f"Saving final sigma calculations to {select_output_file}")
merged_data.to_csv(select_output_file, index=False)

# Save smoothed and trimmed sigma results
logging.info(f"Saving smoothed sigma results to {smoothed_output_file}")
merged_data[['smoothed_sigma', 'trimmed_sigma']].to_csv(smoothed_output_file, index=False)

# Write the final quality counts and details to a separate file
with open(qual_counts_eu_output, 'w') as qual_file:
    qual_file.write(f"File Name: {output_file}\n")
    qual_file.write(f"Bin size: {BIN_SIZE}\n")
    qual_file.write(f"Background Noise (Low): {background_low}\n")
    qual_file.write(f"Background Noise (High): {background_high}\n")
    qual_file.write(f"Baseline Mean (log2 adjusted): {baseline_mean}\n")
    qual_file.write(f"Total sample hits: {total_sample_reads}\n")
    qual_file.write(f"Total adjusted hits: {total_control_reads}\n")
    qual_file.write(f"Correction factor: {correction_factor}\n")

# Step 9: Global and Local Maximum Calculation
def calculate_global_local_max(df, num_bins=500):
    global_max = df['smoothed_sigma'].max()  # Use smoothed_sigma for scaling
    local_max = df['smoothed_sigma'][:num_bins].max()  # First `num_bins` for local scaling
    return global_max, local_max

global_max, local_max = calculate_global_local_max(merged_data)

# Step 10: Plotting Sigma Values as a Bar Plot
def plot_sigma(data, y_max, local_max=None):
    plt.figure(figsize=(12, 6))

    # Plot the smoothed and trimmed sigma values as bars
    plt.bar(data['bin'], data['smoothed_sigma'], label="Smoothed Sigma", color='blue', alpha=0.6)
    plt.bar(data['bin'], data['trimmed_sigma'], label="Trimmed Sigma", color='green', alpha=0.4)

    # Set plot limits based on the chosen maximum value
    plt.ylim(0, y_max)

    # If we have local max data, we can highlight this region
    if local_max:
        plt.axhline(local_max, color='red', linestyle=':', label=f'Local Max ({local_max:.2f})')

    # Plot titles and labels
    plt.title(f"Sigma Values for {sample_basename}")
    plt.xlabel('Bin')
    plt.ylabel('Sigma')
    plt.legend()

    # Show plot
    plt.tight_layout()
    plt.show()

# Step 11: Determine y_max for plot scaling
if manual_max:
    y_max = manual_max  # Use the manually provided maximum if specified
else:
    y_max = global_max  # Otherwise, default to the global maximum

# Plot sigma values as bars
plot_sigma(merged_data, y_max, local_max)

# Save the plot to a file
plt.savefig(output_plot_file)
logging.info(f"Plot saved to {output_plot_file}")
