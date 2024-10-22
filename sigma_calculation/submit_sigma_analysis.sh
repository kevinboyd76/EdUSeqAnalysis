#!/bin/bash

#SBATCH --job-name=sigma_analysis
#SBATCH --cpus-per-task=28              # Number of CPU cores per task
#SBATCH --time=48:00:00                 # Max time for the job (48 hours here)
#SBATCH --mem=128G

# Load necessary modules (adjust based on your HPC environment)
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1

# Run the Python script with input files
python eduseq_sigma_analysis.py ../results/sigma/HCT116_EdU_HU_Aph_set1_adjusted_sample_counts.txt ../results/sigma/HCT116_EdU_HU_Aph_set1_sample_bin_counts.txt ../results/sigma/MTBP6Tir114_TotalSheared_set1_adjust.csv
