import pandas as pd
import os
import sys
importance_file = sys.argv[1]
#list of the files you want to subset your dataframe for different features
list_files = "whole_matrix.txt", "test_enriched_kmers_matrix.txt", "train_enriched_kmers_matrix.txt"
#feature importance file
feature_imp = pd.read_csv(importance_file,sep = "\t", header=0)

important_fet_num = []
x = 10

while x < feature_imp.shape[0]:
	important_fet_num.append(x)
	x = x + 20
for index in important_fet_num:
	imp_fet_mat = feature_imp[["Feature"]].head(index)
	mat_name = "top_features_"+str(index)+".txt"
	imp_fet_mat.to_csv(mat_name, sep = "\t",index = False, header = None)
	#print(imp_fet_mat)
	script_name = str(index)+"_run_model.sh"
	script = open(script_name,"w")
	script.write(f'#!/bin/sh --login\n#SBATCH --time=3:50:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=20G\n#SBATCH --job-name {mat_name}\ncd ./\n')
	script.write(f'export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH\nmodule load Python/3.6.4\n')
	script.write(f'python ~/02_switchgrass_CRE_evolution/Motif_discovary/13_RF_holdout_before_enrichment_draw_two_AUC_SMOTE_upsampling_val_02_feature_selection.py -file train_enriched_kmers_matrix.txt -smote y -cv_num 5 -feat {mat_name}')
	script.close()
	os.system(f'sbatch {script_name}')	
