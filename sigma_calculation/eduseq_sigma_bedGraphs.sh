#!/bin/bash
#SBATCH --job-name=sigma_bedgraphs
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Check if the correct number of arguments are provided
if [[ $# -lt 3 ]]; then
    echo "Usage: sbatch create_bedgraphs.sh <input_csv> <bin_size> <output_dir>"
    exit 1
fi

# Set variables from positional arguments
INPUT_CSV=$1         # The input CSV file
BIN_SIZE=$2          # Bin size in base pairs
OUTPUT_DIR=$3        # Output directory

# Validate that the input CSV file exists
if [[ ! -f "$INPUT_CSV" ]]; then
    echo "Error: Input CSV file does not exist."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to create a bedGraph for a specific column
create_bedgraph() {
    local column=$1
    local column_name=$2
    local output_file="${OUTPUT_DIR}/$(basename "${INPUT_CSV%.csv}")_${column_name}.bedGraph"

    # Skip the header and process the CSV
    awk -v bin_size="$BIN_SIZE" -v col="$column" -F, 'NR > 1 {
        # Extract chromosome, bin, and the specific column value
        chrom=$1;
        bin=$2;
        value=$col;

        # Calculate start and end positions based on the bin index
        start = bin * bin_size;
        end = (bin + 1) * bin_size;

        # Output in bedGraph format: chromosome, start, end, value
        print chrom "\t" start "\t" end "\t" value;
    }' "$INPUT_CSV" > "$output_file"

    # Sort the bedGraph file
    sort -k1,1 -k2,2n "$output_file" > "${output_file%.bedGraph}_sorted.bedGraph"

    echo "bedGraph for ${column_name} created: ${output_file}"
    echo "Sorted bedGraph for ${column_name} created: ${output_file%.bedGraph}_sorted.bedGraph"
}

# Call the function for each column
create_bedgraph 8 "sigma"               # Column 8 is "sigma"
create_bedgraph 9 "sigma_mb"            # Column 9 is "sigma_mb"
create_bedgraph 10 "smoothed_sigma"     # Column 10 is "smoothed_sigma"
create_bedgraph 11 "trimmed_sigma"      # Column 11 is "trimmed_sigma"
#create_bedgraph 12 "sigma_log2"         # Column 12 is "sigma_log2"

echo "All bedGraphs created successfully."
