# EdUSeqAnalysis
 
## clone repository
```
git clone https://github.com/SansamLab/EdUSeqAnalysis.git
```
## load modules
```
module purge
module load slurm python pandas numpy bwa fastqc fastp samtools sicer
```
## run on hpc
```
sbatch --wrap="snakemake -j 999 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output}'"
```
