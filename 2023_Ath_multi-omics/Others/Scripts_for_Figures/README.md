# Descriptions and Usage of Scripts
__feature_selection_plots.ipynb__<br />
Jupyter notebook containing code for generating scatterplots of feature importance scores (RF) or feature coefficients (rrBLUP) for genomic, transcriptomic, and methylomic features.
***
__job_submission.slurm__<br />
SLURM script to submit code as jobs. Modify it depending on the code you want to run, it's and the input data location, and resources you need.
***
__map_SNPs_to_TAIR10_genes.py__<br />
To map the bi-allelic SNPs to the *Arabidopsis thaliana* TAIR10 reference genome.<br />
Input:
- genotype matrix: SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv 
- gff file: Athaliana_167_TAIR10.gene.gff3
- save name: "TAIR10_genes_of_SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv"

Output:
- all_genes_TAIR10.txt
    - TAIR10 reference genome genes with chromosomal location, gene name, start and stop basepair location, and strand (+ or -)
- TAIR10_genes_of_SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv
    - TAIR10 reference genes that are represented in the genotype matrix. Has columns containing SNP ID, chromosome number, basepair position, and gene name
```
# on command line:
python match_SNPs_to_TAIR10_genes.py -geno SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv -gff Athaliana_167_TAIR10.gene.gff3 -save TAIR10_genes_of_SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv 

# or submit a slurm job:
sbatch job_submission.slurm
squeue -u userid # check job status
```
***
__plot_RF_imp_between_omics.R & plot_rrBLUP_coefs_between_omics.R__<br />
R script for making density scatterplots comparing features from genome, transcriptome, and gene body methylome for flowering time (383 accessions) prediction models (Random Forest [RF] and rrBLUP). The actual feature coefficients or importance scores are plotted and the 95% and 99% quantile lines are marked.
***
