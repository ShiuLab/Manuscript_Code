import os
out = open("feature_selection_results.txt","w")
out.write(f'features\tF1_CV\n')
for files in os.listdir("./"):
	if files.startswith("train_enriched_kmers_matrix.txt_RF_cv5_score_roc_auc_y_top_features") and files.endswith(".txt"):
		score_file = open(files,"r").readlines()
		for lines in score_file:
			if lines.startswith("F1_CV:"):
				F1_cv = lines.split(":")[1].lstrip().rstrip()
				feature = files.replace("train_enriched_kmers_matrix.txt_RF_cv5_score_roc_auc_y_top_features_","").strip("")
				feature =feature.replace(".txt","")
				#print(F1_cv,feature)
				out.write(f'{feature}\t{F1_cv}\n')
out.close()	
