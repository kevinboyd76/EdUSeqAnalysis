# EdUSeqAnalysis
 
## 1. Clone repository
```
git clone https://github.com/SansamLab/EdUSeqAnalysis.git
```
## 2. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
## 3. Modify Samples file
```
vim samples.csv
```
## 4. Dry Run
```
snakemake -npr
```
## 5. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 30 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

# Explanation of Final Output
{sample}_sigma_select_EU_0b.csv
- Columns
  - chromosome,bin,adjusted_1,adjusted_2,bin_count_1,bin_count_2,sheared_counts,sigma,sigma_mb,smoothed_sigma,trimmed_sigma,sigma_log2
    •	Chromosome: The chromosome identifier
    •	Bin: Bin number to describe the specific genomic location
    •	Adjusted_1, Adjusted_2: These represent the adjusted hit counts from the Edu-labeled sample, in both forward and reverse directions, for the corresponding bin.
    •	Bin_count_1, Bin_count_2: These represent the unadjusted counts from the Edu-labeled sample for forward and reverse strands, before normalization with the control (total sheared).
    •	Sheared_counts: These are the counts from the total sheared control sample for that bin.
    •	Sigma: This is the initial sigma value, calculated as the ratio of the Edu-labeled sample counts to the total sheared counts, scaled by a factor (SCALE_FACTOR).
    •	Sigma_mb: The sigma value adjusted for background noise, which is normalized using the background noise low and high values.
    •	Smoothed_sigma: The sigma value after smoothing, where a rolling window was applied to the sigma values to remove noise and get a more consistent signal.
    •	Trimmed_sigma: After smoothing, a trimming step was applied to remove extreme outliers (using a trim factor) from the sigma values.
    •	Sigma_log2: The final sigma value transformed to the log2 scale for easier visualization and comparison. The very negative values, such as -29.897352853986263, indicate that the corresponding bin has very low adjusted sigma values.

