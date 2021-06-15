#!/bin/sh --login

#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name rrBLUP.sh
#SBATCH -e rrBLUP.sh.e
#SBATCH -o rrBLUP.sh.o
cd /mnt/home/peipeiw/Documents/Genome_selection/CV/Anthesis_Date_8632_all_top_markers/

module load R
module load Python/3.6.4
#Rscript /mnt/home/peipeiw/Documents/Genome_selection/CV/11_split_geno_pheno_fread.r geno.csv pheno.csv Test.txt
#Rscript /mnt/home/peipeiw/Documents/Genome_selection/CV/09_rrBLUP_fread.r geno_training.csv pheno_training.csv all Anthesis_Date_8632 5 10 CVFs.csv Anthesis_Date_8632_all_top_markers_training_geno
#python /mnt/home/peipeiw/Documents/Genome_selection/CV/12_select_markers_according_to_abs_coef.py -coef Coef_Anthesis_Date_8632_all_top_markers_training_geno_Anthesis_Date_8632.csv -start 250 -stop 18750 -step 250
Rscript /mnt/home/peipeiw/Documents/Genome_selection/CV/09_rrBLUP_fread.r genoMarkers_top18500.txt.csv pheno_training.csv all Anthesis_Date_8632 5 10 CVFs.csv Anthesis_Date_8632_18500_test
