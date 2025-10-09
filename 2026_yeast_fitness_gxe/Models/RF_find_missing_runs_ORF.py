#!/usr/bin/env python3

import os
import re
import pandas as pd

os.chdir('/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/fs')

# check to make sure no old files are present
os.system('ls -l | grep -P Aug')

# results summary file
res = pd.read_csv("RESULTS_reg.txt", sep="\t")

# get the list of feature selection files that did not run
path = '/mnt/home/seguraab/Shiu_Lab/Project/Job_Submission_Scripts/ORFs_as_Feat/Feature_Selection/RF/RF_FS_runs.txt'
runs = open(path, 'r').read().split('\n')

with open(path.replace('RF_FS_runs.txt','RF_FS_missing_runs.txt'), 'w') as f:
    for run in runs[:2800]:
        try:
            feat_file = run.split(' ')[11].split('/fs/')[1]
        except IndexError:
            try:
                feat_file = run.split(' ')[12].split('/fs/')[1]
            except IndexError:
                feat_file = run.split(' ')[13].split('/fs/')[1]
        # check if the _imp file exists (since it's created after all training reps)
        feat_num = re.search(r'top_[0-9]+', feat_file).group()
        env = re.search(r'(?<=_)[A-Z0-9]+(?=_)', feat_file).group()
        dat_type = re.search(r'(?<=features_)[a-z]+', feat_file).group()
        # sometimes an imp file is created but the results don't get written to results file
        if (not os.path.isfile(f'{env}_{dat_type}_{feat_num}_imp')) | \
            (not f"{env}_{dat_type}_{feat_num}" in res.ID.values): # check if the results summary was written to res
            print(run)
            f.write(f'{run}\n')
        del run, feat_file, feat_num, env, dat_type

# June 11. there are 3 missing runs
# July 17. there are no missing pav/cnv runs!