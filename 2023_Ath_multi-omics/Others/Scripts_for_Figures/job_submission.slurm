#!/bin/bash --login 
#SBATCH --time=5:00:00                         # Allocated run time
#SBATCH --nodes=1                               # Number of nodes
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=2                       # Number of CPUs per task
#SBATCH --mem-per-cpu=50G                       # Allocate CPU memory
#SBATCH --job-name map_snps_to_TAIR10_genes     # Job name
#SBATCH --output=%x_%j                          # Output file in the format <job_name_slurm_job_id>

### Load relevant modules ###
module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  Python/3.6.4

### Directory to data ###
cd /mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Data

### Code to run ###
python /mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/multi-omics/Scripts_Figures/map_SNPs_to_TAIR10_genes.py -geno SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv -gff Athaliana_167_TAIR10.gene.gff3 -save TAIR10_genes_of_SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv

scontrol show job $SLURM_JOB_ID # print job information
