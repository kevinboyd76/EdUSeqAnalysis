# EdUSeqAnalysis
![EdU](/images/EdUseq.png)

# 1) Project Description
EdUSeqAnalysis Pipeline is a Snakemake pipeline designed to analyze whole-genome sequencing data from Edu-labeled DNA samples. This pipeline generates trimmed FASTQ files, genome alignments (hg38), coverage files, and background-adjusted sigma values, providing a robust framework for analyzing Edu-DNA incorporation. The sigma values were calculated, trimmed, and smoothed based on methologies adapted from Macheret and Halazonetis (Nature, 2018). Each processing step is clearly defined, with dependencies managed through Snakemake, and execution automated via module environments.

The pipeline utilizes a control sample to normalize EdU-DNA counts, producing files compatible with genome visualization and quantitative analysis of DNA replication and synthesis. To facilitate quick testing, we include a compact dataset within the repository. Additionally, a detailed example is provided to demonstrate how to run this pipeline. This workflow is inspired by and extends the protocols provided by the Sansam Lab and Macheret and Halazonetis, particularly through modifications in genome alignment, sigma calculation, background subtraction, and smoothing techniques aligned to the hg38 genome assembly.


# 2) Instructions to Run Pipeline on Slurm Managed HPC
## 2A. Clone repository
```
git clone https://github.com/SansamLab/EdUSeqAnalysis.git
```
## 2B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
## 2C. Modify Samples file
```
vim samples.csv
```
## 2D. Dry Run
```
snakemake -npr
```
## 2E. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --latency-wait 30 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```

# 3) Explanation of samples.csv
Note. Make sure to check sample.csv before each run

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files. Columns 2 and 3 are the paths to the paired fastq files. Column 4 is the sample type (either "treatment" or "control"). Column 5 is the name of the corresponding Control sample for each treated sample (use "NA" if the sample is a control).

| sample      | fastq1              | fastq2              | sampleType | Control   |
|-------------|---------------------|---------------------|------------|-----------|
| testSample  | sample_R1.fastq.gz  | sample_R2.fastq.gz  | treatment  | testInput |
| testInput   | sample2_R1.fastq.gz | sample2_R2.fastq.gz | control    | NA        |


# 4) Explanation of Final Output
{sample}_sigma_select_EU_0b.csv
- Columns: chromosome, bin, adjusted_1, adjusted_2, bin_count_1, bin_count_2, sheared_counts, sigma, sigma_mb, smoothed_sigma, trimmed_sigma, sigma_log2
    +	Chromosome: The chromosome identifier for each bin, aligned to hg38
    +	Bin: Bin number to describe the specific genomic location, set at 10,000 base pairs
    +	Adjusted_1, Adjusted_2: These represent the adjusted read counts for the Edu-labeled sample in the forward and reverse directions, respectively, for each bin. These values are derived by normalizing the original counts from the Edu sample against the control sample (total sheared DNA).
    +	Bin_count_1, Bin_count_2: The raw, unadjusted counts from the Edu-labeled sample in forward and reverse directions before any normalization. These counts give an initial measure of signal intensity for each strand in each bin.
    +	Sheared_counts: The counts from the total sheared control sample for each bin, representing background or baseline DNA levels for comparison.
    +	Sigma: The initial sigma value is calculated as the ratio of Edu-labeled sample counts to total sheared counts, adjusted by a scaling factor (SCALE_FACTOR). This value reflects the relative enrichment of EdU-labeled DNA in each bin before further background correction.
    +	Sigma_mb: The background-adjusted sigma value, which is normalized using the background noise thresholds calculated from the low and high percentiles. This adjustment helps standardize the sigma values by reducing the impact of noisy bins with low background counts.
    +	Smoothed_sigma: The sigma value after percentile-based smoothing, where outliers and background noise are reduced based on selected percentiles. This percentile-based approach yields a stable and consistent signal.
    +	Trimmed_sigma: Post-smoothing, a trimming step is applied to further reduce extreme outliers, using a trim factor to cap extreme deviations.
    +	Sigma_log2: The final sigma value transformed to the log2 scale for better visualization and comparison. Very negative values indicate bins with low or near-zero adjusted sigma values.


# 5) Citations
Macheret, M., & Halazonetis, T. D. (2018). Intragenic origins due to short G1 phases underlie oncogene-induced DNA replication stress. Nature, 555(7694), 112–116. https://doi.org/10.1038/nature25507

