#!/usr/bin/env python

import pandas as pd
import numpy as np
import logging
import sys

# Set up logging
logging.basicConfig(filename='edu_seq_analysis.log', level=logging.INFO)

# Input variables
BIN_SIZE = 10000
CORRECTION_FACTOR = 1.0  # Placeholder for correction factor, can be adjusted later

# Load the input file names from command-line arguments
adjusted_counts_file = sys.argv[1]  # EduHU_HCT_Biotin_set2A_adjusted_sample_counts.txt
bin_counts_file = sys.argv[2]       # EduHU_HCT_Biotin_set2A_sample_bin_counts.txt
totalsheared_file = sys.argv[3]     # EduHU_HCT_TotalSheared_set2A_adjust.csv

output_file = "sigma_output.csv"
qual_counts_output = "sigma_qual_counts.txt"

# Read the three input files
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
        # Sigma calculation (simplified): (adjusted bin counts / total sheared counts)
        sigma = row['bin_count_1'] / row['sheared_counts']
        merged_data.at[index, 'sigma'] = sigma

# Save the sigma output to a new CSV file
logging.info(f"Saving sigma calculations to {output_file}")
merged_data.to_csv(output_file, index=False)

# Write the quality counts to a separate file
total_sample_hits = merged_data['bin_count_1'].sum()
total_adjust_hits = merged_data['sheared_counts'].sum()

with open(qual_counts_output, 'w') as qc_file:
    qc_file.write(f"Bin size: {BIN_SIZE}\n")
    qc_file.write(f"Total sample hits: {total_sample_hits}\n")
    qc_file.write(f"Total adjusted hits: {total_adjust_hits}\n")
    qc_file.write(f"Correction factor: {CORRECTION_FACTOR}\n")

logging.info(f"Sigma calculations complete.")

