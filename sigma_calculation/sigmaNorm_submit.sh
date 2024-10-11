#!/bin/bash

#SBATCH --job-name=edu_seq_sigma
#SBATCH --cpus-per-task=28              # Number of CPU cores per task
#SBATCH --time=48:00:00                 # Max time for the job (48 hours here)
#SBATCH --mem=128G

# Load necessary modules (adjust based on your HPC environment)
module load python

# Run the Python script with input files
python sigmaCalc.py \
    sigma_output.csv
