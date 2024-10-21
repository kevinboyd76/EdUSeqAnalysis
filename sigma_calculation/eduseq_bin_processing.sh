#!/bin/bash
#SBATCH --job-name=eduseq_bin_processing
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Load necessary modules
module load samtools
module load bedtools
module load deeptools

# Define variables passed from Snakemake
BIN_SIZE=10000  # Bin size in base pairs
QUALITY=60      # Minimum mapping quality
CHROM_SIZES="/hpc-prj/cobre-dev-bio/boydk/reference_info/hg38/hg38.chrom.sizes"  # Update this to the actual location

# Input SAM files
ADJUST_SAM=$1  # Total sheared DNA SAM file
SAMPLE_SAM=$2  # Edu-labeled SAM file
WORK_DIR=$3    # Working directory passed from Snakemake

# Check if the input files are provided
if [[ $# -lt 3 ]]; then
    echo "Usage: sbatch eduseq_pipeline.sh <adjust_file.sam> <sample_file.sam> <work_dir>"
    exit 1
fi

# Validate that the input files exist
if [[ ! -f "$ADJUST_SAM" || ! -f "$SAMPLE_SAM" ]]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Generate file prefixes for output
ADJUST_PREFIX=$(basename "$ADJUST_SAM" .sam)
SAMPLE_PREFIX=$(basename "$SAMPLE_SAM" .sam)

# Step 1: Initial input and variable setup
echo "Step 1 complete. Input files and variables are correctly set."

# Step 2: Generate the Adjust File from Total Sheared DNA Sample
echo "Step 2: Generating Adjust File from $ADJUST_SAM..."
samtools view -h -q $QUALITY "$ADJUST_SAM" | \
    awk -v bin_size="$BIN_SIZE" '
    BEGIN { FS = "\t"; }
    !/^@/ {
        chr = $3;
        bin = int($4 / bin_size);
        if (($2 % 32) >= 16) {
            reverse_hits[chr][bin]++;
        } else {
            forward_hits[chr][bin]++;
        }
    }
    END {
        for (chr in forward_hits) {
            for (bin in forward_hits[chr]) {
                print chr, bin, forward_hits[chr][bin], reverse_hits[chr][bin];
            }
        }
    }' > "${WORK_DIR}/${ADJUST_PREFIX}_adjust_bin_counts.txt"
echo "Step 2 complete. Adjust file counts saved to ${WORK_DIR}/${ADJUST_PREFIX}_adjust_bin_counts.txt."

# Step 3: Convert bin counts to CSV format
echo "Step 3: Converting Adjust bin counts to CSV format..."
awk '
{
    adjust_hits = $3 + $4;  # sum of forward and reverse hits
    print $1 "," $2 "," adjust_hits;
}' "${WORK_DIR}/${ADJUST_PREFIX}_adjust_bin_counts.txt" > "${WORK_DIR}/${ADJUST_PREFIX}_adjust.csv"
echo "Step 3 complete. Adjust CSV file created: ${WORK_DIR}/${ADJUST_PREFIX}_adjust.csv."

# Step 4: Generate bin counts for Edu-labeled sample
echo "Step 4: Generating bin counts for Edu-labeled sample $SAMPLE_SAM..."
samtools view -h -q $QUALITY "$SAMPLE_SAM" | \
    awk -v bin_size="$BIN_SIZE" '
    BEGIN { FS = "\t"; }
    !/^@/ {
        chr = $3;
        bin = int($4 / bin_size);
        if (($2 % 32) >= 16) {
            reverse_hits[chr][bin]++;
        } else {
            forward_hits[chr][bin]++;
        }
    }
    END {
        for (chr in forward_hits) {
            for (bin in forward_hits[chr]) {
                print chr, bin, forward_hits[chr][bin], reverse_hits[chr][bin];
            }
        }
    }' > "${WORK_DIR}/${SAMPLE_PREFIX}_sample_bin_counts.txt"
echo "Step 4 complete. Sample bin counts saved to ${WORK_DIR}/${SAMPLE_PREFIX}_sample_bin_counts.txt."

# Step 5: Adjust sample bin counts using the adjust file
echo "Step 5: Adjusting sample bin counts using the adjust file..."
awk -F'[ ,]' 'NR==FNR { adjust[$1","$2] = $3; next }
{
    bin = $1","$2;
    adjust_hits = (bin in adjust) ? adjust[bin] : 0.01;  # Small constant instead of 1 to avoid division by zero
    adjbin_f = int($3 * 1000 / adjust_hits + 0.5);
    adjbin_r = int($4 * 1000 / adjust_hits + 0.5);
    print $1, $2, adjbin_f, adjbin_r;
}' "${WORK_DIR}/${ADJUST_PREFIX}_adjust.csv" "${WORK_DIR}/${SAMPLE_PREFIX}_sample_bin_counts.txt" > "${WORK_DIR}/${SAMPLE_PREFIX}_adjusted_sample_counts.txt"
echo "Step 5 complete. Adjusted sample counts saved to ${WORK_DIR}/${SAMPLE_PREFIX}_adjusted_sample_counts.txt."

# Step 6: Create a Bed file from adjusted counts
echo "Step 6: Combining forward and reverse hits to create Bed file..."
awk -v bin_size="$BIN_SIZE" '
{
    total_hits = $3 + $4;  # Combine forward and reverse hits
    print $1 "\t" ($2 * bin_size) "\t" (($2 + 1) * bin_size) "\t" total_hits;
}' "${WORK_DIR}/${SAMPLE_PREFIX}_adjusted_sample_counts.txt" > "${WORK_DIR}/${SAMPLE_PREFIX}.bed"

# Sort the BED file
sort -k1,1 -k2,2n "${WORK_DIR}/${SAMPLE_PREFIX}.bed" > "${WORK_DIR}/${SAMPLE_PREFIX}_sorted.bed"
echo "Step 6 complete. Bed files created: ${WORK_DIR}/${SAMPLE_PREFIX}.bed, ${WORK_DIR}/${SAMPLE_PREFIX}_sorted.bed."

# Step 7: Create a bedGraph file from the adjusted counts
echo "Step 7: Creating bedGraph file from adjusted counts..."

awk -v bin_size="$BIN_SIZE" '
{
    total_hits = $3 + $4;  # Combine forward and reverse hits
    print $1 "\t" ($2 * bin_size) "\t" (($2 + 1) * bin_size) "\t" total_hits;
}' "${WORK_DIR}/${SAMPLE_PREFIX}_adjusted_sample_counts.txt" > "${WORK_DIR}/${SAMPLE_PREFIX}.bedGraph"

# Sort the bedGraph
sort -k1,1 -k2,2n "${WORK_DIR}/${SAMPLE_PREFIX}.bedGraph" > "${WORK_DIR}/${SAMPLE_PREFIX}_sorted.bedGraph"

echo "Step 7 complete. bedGraph created from adjusted counts."

# Step 8: Convert sorted bedGraph to BigWig
echo "Step 8: Creating BigWig file from sorted bedGraph..."

bedGraphToBigWig "${WORK_DIR}/${SAMPLE_PREFIX}_sorted.bedGraph" "$CHROM_SIZES" "${WORK_DIR}/${SAMPLE_PREFIX}.bw"

if [[ -f "${WORK_DIR}/${SAMPLE_PREFIX}.bw" ]]; then
    echo "Step 8 complete. BigWig file created: ${WORK_DIR}/${SAMPLE_PREFIX}.bw"
else
    echo "Error: BigWig file creation failed."
fi

# Final message
echo "EduSeq analysis pipeline complete."
