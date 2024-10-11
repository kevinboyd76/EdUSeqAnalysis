#!/bin/bash

#SBATCH --job-name=edu_seq_sigma
#SBATCH --cpus-per-task=28              # Number of CPU cores per task
#SBATCH --time=48:00:00                 # Max time for the job (48 hours here)
#SBATCH --mem=128G

# Load necessary modules (adjust based on your HPC environment)
module load python

# Run the Python script with input files
python sigmaCalc.py \
    EduHU_HCT_Biotin_set2A_adjusted_sample_counts.txt \
    EduHU_HCT_Biotin_set2A_sample_bin_counts.txt \
    EduHU_HCT_TotalSheared_set2A_adjust.csv

