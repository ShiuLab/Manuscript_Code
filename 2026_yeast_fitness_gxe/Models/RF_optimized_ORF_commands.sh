#!/usr/bin/bash

conda activate ml-pipeline

traits=($(head -n 1 /mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv | sed 's/,/ /g'))
fs_path=/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF
models=/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results

############### Create feature selection subsets for each model ################
for i in {1..35}; do
    ## PAV RF Models
    python ${fs_path}/1_feature_selection.py \
        -i ${models}/baseline/${traits[$i]}_pav_baseline_imp \
        -start 0 \
        -stop 7500 \
        -step 250 \
        -o ${models}/fs/features_pav_${traits[$i]} -d n
    
    python ${fs_path}/1_feature_selection.py \
        -i ${models}/baseline/${traits[$i]}_pav_baseline_imp \
        -start 1 \
        -stop 10 \
        -base 2 \
        -o ${models}/fs/features_pav_${traits[$i]} -d n

    ## CNV RF Models
    python ${fs_path}/1_feature_selection.py \
        -i ${models}/baseline/${traits[$i]}_cnv_baseline_imp \
        -start 0 \
        -stop 7500 \
        -step 250 \
        -o ${models}/fs/features_cnv_${traits[$i]} -d n
    
    python ${fs_path}/1_feature_selection.py \
        -i ${models}/baseline/${traits[$i]}_cnv_baseline_imp \
        -start 1 \
        -stop 10 \
        -base 2 \
        -o ${models}/fs/features_cnv_${traits[$i]} -d n
done

########### Create file to call RF pipeline for each feature subset ############
pipe=/mnt/home/seguraab/Shiu_Lab/Project/External_software/ML-Pipeline
data=/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018

for j in {1..35}; do
    ## PAVs non-exponential features
    feat_pav=($(ls ${models}/fs/features_pav_${traits[${j}]}*))
    for feat in ${feat_pav[@]}; do
        feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
        echo "python ${pipe}/ML_regression.py -df ${data}/ORFs_pres_abs.csv -df2 ${data}/pheno.csv -y_name ${traits[${j}]} -sep , -feat ${feat} -test ${data}/Test.txt -alg RF -n_jobs 12 -n 20 -cv_num 5 -save ${traits[${j}]}_pav_top_${feat_num} -plots t" >> RF_FS_runs.txt
    done

    ## CNVs non-exponential features
    feat_cnv=($(ls ${models}/fs/features_cnv_${traits[${j}]}*))
    for feat in ${feat_cnv[@]}; do
        feat_num=$(echo ${feat} | grep -P [0-9]+$ -o)
        echo "python ${pipe}/ML_regression.py -df ${data}/ORFs_no_NA.csv -df2 ${data}/pheno.csv -y_name ${traits[${j}]} -sep , -feat ${feat} -test ${data}/Test.txt -alg RF -n_jobs 12 -n 20 -cv_num 5 -save ${traits[${j}]}_cnv_top_${feat_num} -plots t" >> RF_FS_runs.txt
    done
done