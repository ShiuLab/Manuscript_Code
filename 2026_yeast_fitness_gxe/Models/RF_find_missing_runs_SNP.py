#!/usr/bin/env python3

import os
import re
import pandas as pd

os.chdir('/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs')

# results summary file
res = pd.read_csv("RESULTS_reg.txt", sep="\t")
res = res.loc[res.DateTime.str.contains("2024-06")|res.DateTime.str.contains("2024-07"),:]

# get the list of feature selection files that did not run
path = '/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/SNPs_as_Feat/Feature_Selection/RF/RF_FS_runs.txt'
runs = open(path, 'r').read().split('\n')

with open(path.replace('RF_FS_runs.txt','RF_FS_missing_runs.txt'), 'w') as f:
    for run in runs[:1400]:
        feat_file = run.split(' ')[11].split('/fs/')[1]
        # check if the _imp file exists (since it's created after all training reps)
        feat_num = re.search(r'top_[0-9]+', feat_file).group()
        env = re.search(r'(?<=_)[A-Z0-9]+(?=_)', feat_file).group()
        # sometimes an imp file is created but the results don't get written to results file
        if (not os.path.isfile(f'{env}_rf_{feat_num}_imp')) | \
            (not f"{env}_rf_{feat_num}" in res.ID.values): # check if the results summary was written to res
            print(run)
            f.write(f'{run}\n')
        del run, feat_file, feat_num, env

# July 9: 7 missing runs; already submitted
# July 14: No missing runs!