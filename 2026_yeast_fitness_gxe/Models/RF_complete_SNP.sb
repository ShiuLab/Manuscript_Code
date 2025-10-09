#!/bin/bash --login 
#SBATCH --array=33-34
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=6G
#SBATCH --job-name complete_rf_snp
#SBATCH --output=%x_%j.out

cd /mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/baseline

module purge
conda activate /mnt/home/seguraab/miniconda3/envs/ml-pipeline
TRAITS=($(head -n 1 ~/Shiu_Lab/Project/Data/Peter_2018/pheno.csv | sed 's/,/ /g'))

echo "${SLURM_ARRAY_TASK_ID} ; ${TRAITS[${SLURM_ARRAY_TASK_ID}]}"

path1=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018
python ${path1}/ML_regression.py -df ${data}/geno.csv -df2 ${data}/pheno.csv -y_name ${TRAITS[${SLURM_ARRAY_TASK_ID}]} -sep , -test ${data}/Test.txt -alg RF -n_jobs 12 -n 20 -cv_num 5 -save ${TRAITS[${SLURM_ARRAY_TASK_ID}]}_rf_complete -plots t

conda deactivate

scontrol show job $SLURM_JOB_ID



