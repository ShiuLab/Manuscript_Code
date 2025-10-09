"""Get the best model of the n training repetitions based on the validation set performance"""

import joblib
import os
import pandas as pd
from sklearn.metrics import r2_score

mapping = {"YPACETATE": "YP Acetate 2%", "YPD14": "YPD 14ºC", "YPD40": "YPD 40ºC",
           "YPD42": "YPD 42ºC", "YPD6AU": "YPD 6-Azauracile 600 µg/ml",
           "YPDANISO10": "YPD Anisomycin 10 µg/ml", "YPDANISO20": "YPD Anisomycin 20 µg/ml",
           "YPDANISO50": "YPD Anisomycin 50 µg/ml", "YPDBENOMYL200": "YPD Benomyl 200 µg/ml",
           "YPDBENOMYL500": "YPD Benomyl 500 µg/ml", "YPDCAFEIN40": "YPD Caffeine 40 mM",
           "YPDCAFEIN50": "YPD Caffeine 50 mM", "YPDCHX05": "YPD Cycloheximide 0.5 µg/ml",
           "YPDCHX1": "YPD Cycloheximide 1 µg/ml", "YPDCUSO410MM": "YPD CuSO4 10 mM",
           "YPDDMSO": "YPD DMSO 6%", "YPDETOH": "YPD Ethanol 15%",
           "YPDFLUCONAZOLE": "YPD Fluconazole 20 µg/ml", "YPDFORMAMIDE4": "YPD Formamide 4%",
           "YPDFORMAMIDE5": "YPD Formamide 5%", "YPDHU": "YPD Hydroxyurea 30 mg/ml",
           "YPDKCL2M": "YPD KCL 2 M", "YPDLICL250MM": "YPD LiCl 250 mM",
           "YPDMV": "YPD Methylviologen 20 mM", "YPDNACL15M": "YPD NaCl 1.5 M",
           "YPDNACL1M": "YPD NaCl 1 M", "YPDNYSTATIN": "YPD Nystatin 10 µg/ml",
           "YPDSDS": "YPD SDS 0.2%", "YPDSODIUMMETAARSENITE": "YPD Sodium metaarsenite 2.5 mM",
           "YPETHANOL": "YP Ethanol 2%", "YPGALACTOSE": "YP Galactose 2%",
           "YPRIBOSE": "YP Ribose 2%", "YPGLYCEROL": "YP Glycerol 2%",
           "YPXYLOSE": "YP Xylose 2%", "YPSORBITOL": "YP Sorbitol 2%"}

test_set = pd.read_csv(
    "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=None)

########################## PC/SNP RF complete models ###########################
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project"
# list of best models
snp_out = open(f"{d}/SNP_yeast_RF_results/baseline/snp_best_models.txt", "a")
pc_out = open(f"{d}/PC_yeast_RF_results/pc_best_models.txt", "a")
for env in mapping.keys():
    try:
        # read in predicted labels for each training rep
        snp = pd.read_csv(
            f"{d}/SNP_yeast_RF_results/baseline/{env}_rf_baseline_scores.txt", sep="\t", index_col=0)
        pc = pd.read_csv(
            f"{d}/PC_yeast_RF_results/{env}_PCs_tassel_scores.txt", sep="\t", index_col=0)
        # select only test instances
        snp_val = snp.loc[~snp.index.isin(test_set[0]), :]
        pc_val = pc.loc[~pc.index.isin(test_set[0]), :]
        # Calculate r2 performances
        r2_snp_val = snp_val.iloc[:, 3:].apply(
            lambda x: r2_score(snp_val.Y, x), axis=0)
        r2_pc_val = pc_val.iloc[:, 3:].apply(
            lambda x: r2_score(pc_val.Y, x), axis=0)
        # Determine which model has the maximum performance on the training data
        snp_model = int(r2_snp_val.idxmax().split("_")[1])
        pc_model = int(r2_pc_val.idxmax().split("_")[1])
        snp_out.write(f"{env}_rf_baseline_models_rep_{snp_model-1}.pkl\n")
        pc_out.write(
            f"{env}_PCs_tassel_baseline_models_rep_{pc_model-1}.pkl\n")
        del snp, pc
    except:
        print(env)
        continue

snp_out.close()
pc_out.close()

########################## PAV/CNV RF complete models ##########################
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/baseline"
pav_out = open(f"{d}/pav_best_models.txt", "a")  # list of best models
cnv_out = open(f"{d}/cnv_best_models.txt", "a")
for env in mapping.keys():
    # read in predicted labels for each training rep
    pav = pd.read_csv(f"{d}/{env}_pav_baseline_scores.txt",
                      sep="\t", index_col=0)
    cnv = pd.read_csv(f"{d}/{env}_cnv_baseline_scores.txt",
                      sep="\t", index_col=0)
    # select only validation instances
    pav_val = pav.loc[~pav.index.isin(test_set[0]), :]
    cnv_val = cnv.loc[~cnv.index.isin(test_set[0]), :]
    # Calculate r2 performances
    r2_pav_val = pav_val.iloc[:, 3:].apply(
        lambda x: r2_score(pav_val.Y, x), axis=0)
    r2_cnv_val = cnv_val.iloc[:, 3:].apply(
        lambda x: r2_score(cnv_val.Y, x), axis=0)
    # Determine which model has the maximum performance on the training data
    pav_model = int(r2_pav_val.idxmax().split("_")[1])
    cnv_model = int(r2_cnv_val.idxmax().split("_")[1])
    pav_out.write(f"{env}_pav_baseline_models_rep_{pav_model-1}.pkl\n")
    cnv_out.write(f"{env}_cnv_baseline_models_rep{cnv_model-1}.pkl\n")
    del pav, cnv

pav_out.close()
cnv_out.close()

####################### SNP RF feature selection models ########################
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/fs"
snp_out = open(f"{d}/snp_best_fs_models.txt", "a")  # list of best models
res = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t")
for env in mapping.keys():
    try:
        res_env = res.loc[res.Y == env, :].sort_values(by="FeatureNum")
        # res_env = res_env.loc[res_env.FeatureNum <= 30000,:]
        print(env, res_env.shape)
        # Get the exact model repetition
        id = res_env.r2_val.idxmax()
        opt = res_env.loc[res_env.index == id, "FeatureNum"]
        mod_preds = pd.read_csv(
            f"{d}/{env}_rf_top_{str(opt.values[0])}_scores.txt", sep="\t", index_col=0)
        mod_preds = mod_preds.loc[~mod_preds.index.isin(test_set[0]), :]
        r2_val = mod_preds.iloc[:, 3:].apply(
            lambda x: r2_score(mod_preds.Y, x), axis=0)
        mod_id = int(r2_val.idxmax().split("_")[1])
        snp_out.write(
            f"{env}_rf_top_{str(opt.values[0])}_models_rep_{mod_id-1}.pkl\n")
    except:
        print(env)
        continue

snp_out.close()

####################### PAV/CNV RF feature selection models ########################
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/fs"
res = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t")
res = res.iloc[res.ID.drop_duplicates(keep="last").index]

missing = []  # missing runs
for file in os.listdir(d):
    if file.startswith('features'):
        dattype = file.split('_')[1]
        env = file.split('_')[2]
        featnum = '_'.join(file.split('_')[3:5])
        if f'{env}_{dattype}_{featnum}' not in res.ID.values:
            missing.append(file)

len(missing)

res.ID.nunique()

for data in ["pav", "cnv"]:
    # for env in mapping.keys(): # check which envs have missing models
    #     res_env = res.loc[res.ID.str.contains(f'{env}_{data}'),:].sort_values(by="FeatureNum")
    #     print(data, env, res_env.shape)
    out = open(f"{d}/{data}_best_fs_models.txt", "a")  # list of best models
    for i, env in enumerate(mapping.keys()):
        # for i,env in enumerate(['YPDCAFEIN40', 'YPDCAFEIN50', 'YPDCUSO410MM', 'YPDSODIUMMETAARSENITE', 'YPDBENOMYL500']):
        # Generate feature selection curve
        res_env = res.loc[res.Y == env, :].sort_values(by="FeatureNum")
        res_env = res_env.loc[res_env.ID.str.contains(data), :]
        print(env, res_env.shape)
        # Get the exact model repetition
        id = res_env.r2_val.idxmax()
        opt = res_env.loc[res_env.index == id, "FeatureNum"]
        mod_preds = pd.read_csv(
            f"{d}/{env}_{data}_top_{str(opt.values[0])}_scores.txt", sep="\t", index_col=0)
        mod_preds = mod_preds.loc[~mod_preds.index.isin(test_set[0]), :]
        r2_val = mod_preds.iloc[:, 3:].apply(
            lambda x: r2_score(mod_preds.Y, x), axis=0)
        mod_id = int(r2_val.idxmax().split("_")[1])
        out.write(
            f"{env}_{data}_top_{str(opt.values[0])}_models_rep_{mod_id-1}.pkl\n")
    out.close()
