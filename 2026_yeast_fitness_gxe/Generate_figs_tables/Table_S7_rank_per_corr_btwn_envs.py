#!/usr/bin/env python3
################################################################################
# TABLE S7: Spearman's rho of feature ranks between environments
################################################################################
import os
import pandas as pd
import datatable as dt
import numpy as np
from scipy.stats import spearmanr
from itertools import combinations

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Isolate growth condition labels; will be used throughout the script
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

for i, mod_type in enumerate(["complete", "optimized"]):
    for j, data_type in enumerate(["snp", "pav", "cnv"]):
        for k, imp_type in enumerate(["gini", "shap"]):
            print(
                f"Calculating spearman's rho for Scripts/Data_Vis/Section_3/RF_{mod_type}_{imp_type}_{data_type}.tsv")
            df = dt.fread(
                f"Scripts/Data_Vis/Section_3/RF_{mod_type}_{imp_type}_{data_type}.tsv").to_pandas()
            df.set_index(df.columns[0], inplace=True)
            #
            # calculate the ranks for each feature in each env
            df[df == 0] = np.nan  # set 0 importance to NAs
            df = df[df.index != ""]  # remove the row with empty feature name
            if data_type == "snp":
                df.iloc[:, 2:] = df.iloc[:, 2:].rank(axis=0, method="average",
                                                     pct=False, ascending=False)
                # calculate spearman's rho
                rho = df.iloc[:, 2:].corr(
                    method=lambda x, y: spearmanr(x, y).statistic)
                pval = df.iloc[:, 2:].corr(
                    method=lambda x, y: spearmanr(x, y).pvalue)
                env_pairs = list(combinations(df.columns[2:], 2))
            else:
                df.iloc[:, 1:] = df.iloc[:, 1:].rank(axis=0, method="average",
                                                     pct=False, ascending=False)
                # calculate spearman's rho
                rho = df.iloc[:, 1:].corr(
                    method=lambda x, y: spearmanr(x, y).statistic)
                pval = df.iloc[:, 1:].corr(
                    method=lambda x, y: spearmanr(x, y).pvalue)
                env_pairs = list(combinations(df.columns[1:], 2))
            #
            # get the number of common features between envs
            shared_feats = {}
            for env1, env2 in env_pairs:
                shared_feats[(env1, env2)] = df[[env1, env2]].dropna().shape[0]
            #
            # prepare table to save
            print(rho.shape, pval.shape)
            out = rho.where(np.triu(rho, k=0).astype(bool)).stack()
            out2 = pval.where(np.triu(pval, k=0).astype(bool)).stack()
            out = pd.concat([out, out2], ignore_index=False, axis=1)
            out.columns = ["rho", "pval"]
            # 0.0 pvalues become NaNs with astype(bool)
            out.index.set_names(["Env1", "Env2"], inplace=True)
            print(out.shape)
            if (i == 0 and j == 0 and k == 0):
                res = out.copy(deep=True)
                res.columns = [f"rho_{imp_type}_{data_type}_{mod_type}",
                               f"pval_{imp_type}_{data_type}_{mod_type}"]
                shared_feats_dict = {
                    f"NumSharedFeats_{imp_type}_{data_type}_{mod_type}": shared_feats}
            else:
                out.columns = [f"rho_{imp_type}_{data_type}_{mod_type}",
                               f"pval_{imp_type}_{data_type}_{mod_type}"]
                res = pd.concat([res, out], ignore_index=False, axis=1)
                shared_feats_dict[f"NumSharedFeats_{imp_type}_{data_type}_{mod_type}"] = shared_feats
            del df, rho, pval, out

# Order environments based on the Figure 1A clustering and save
env_order = ["YPDCHX05", "YPDCHX1", "YPDANISO50", "YPDANISO10", "YPDANISO20",
             "YPDDMSO", "YPDMV", "YPDSDS", "YPD40", "YPD42", "YPDKCL2M",
             "YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL200", "YPDBENOMYL500",
             "YPDETOH", "YPDNYSTATIN", "YPACETATE", "YPXYLOSE", "YPRIBOSE",
             "YPSORBITOL", "YPGLYCEROL", "YPETHANOL", "YPGALACTOSE",
             "YPDLICL250MM", "YPDNACL15M", "YPDNACL1M", "YPDFORMAMIDE4",
             "YPDFORMAMIDE5", "YPDHU", "YPD14", "YPDFLUCONAZOLE",
             "YPDSODIUMMETAARSENITE", "YPD6AU", "YPDCUSO410MM"]

res.reset_index(inplace=True)
res = res.set_index("Env2").loc[env_order]
res = res.reset_index().set_index("Env1").loc[env_order]
res.reset_index(inplace=True)
res.to_csv(
    "Scripts/Data_Vis/Section_3/Table_S7_rank_per_corr_btwn_envs.tsv", sep="\t")

feats = pd.DataFrame(shared_feats_dict)
feats = pd.concat([feats, pd.DataFrame([[np.nan] * 12],
                                       index=[("YPSORBITOL", "YPACETATE")],
                                       columns=feats.columns)], axis=0)
feats.index.names = ["Env1", "Env2"]
feats = feats.reset_index().set_index("Env2").loc[env_order]
feats = feats.reset_index().set_index("Env1").loc[env_order]
feats.to_csv(
    "Scripts/Data_Vis/Section_3/Table_S7_rank_per_corr_btwn_envs_shared_feat_counts.tsv", sep="\t")


''' Summary statistics:
>>> res.iloc[:,2:].describe().T
                         count      mean       std            min            25%            50%           75%  max
rho_gini_snp_complete    630.0  0.211503  0.235630  -8.336781e-02   6.939274e-02   1.441885e-01  3.028398e-01  1.0
pval_gini_snp_complete   349.0  0.142834  0.325466  6.888430e-308  3.692431e-121   1.057203e-23  2.074783e-03  1.0
rho_shap_snp_complete    630.0  0.099006  0.236031  -3.776224e-01  -1.517539e-03   4.360521e-02  1.030239e-01  1.0
pval_shap_snp_complete   627.0  0.242337  0.334107  4.271594e-281   1.521418e-06   4.286226e-02  4.276967e-01  1.0
rho_gini_pav_complete    630.0  0.727196  0.134374   2.569054e-01   6.637737e-01   7.517486e-01  8.030471e-01  1.0
pval_gini_pav_complete   230.0  0.152174  0.359973  5.928788e-323  4.713557e-232  2.511940e-146  1.404615e-37  1.0
rho_shap_pav_complete    630.0  0.566188  0.204012   6.129416e-02   4.476659e-01   6.002723e-01  6.980137e-01  1.0
pval_shap_pav_complete   598.0  0.060033  0.234917  3.952525e-323  1.659305e-162   1.438386e-84  7.851175e-28  1.0
rho_gini_cnv_complete    630.0  0.622925  0.150051   2.525476e-01   5.355634e-01   6.359437e-01  7.008104e-01  1.0
pval_gini_cnv_complete   170.0  0.205882  0.405539  1.088845e-310  2.161353e-206   1.734329e-88  5.065225e-59  1.0
rho_shap_cnv_complete    630.0  0.451203  0.196598   3.883042e-02   3.239675e-01   4.715941e-01  5.449895e-01  1.0
pval_shap_cnv_complete   580.0  0.065873  0.241197  9.021639e-320  1.205849e-150   8.097407e-75  8.607170e-12  1.0
rho_gini_snp_optimized   414.0  0.223253  0.564980  -1.000000e+00  -1.019552e-01   2.500000e-01  6.236347e-01  1.0
pval_gini_snp_optimized  352.0  0.484194  0.338397   1.611489e-93   1.820725e-01   4.927262e-01  7.735994e-01  1.0
rho_shap_snp_optimized   415.0  0.282608  0.550708  -1.000000e+00  -8.902610e-03   3.000000e-01  7.402930e-01  1.0
pval_shap_snp_optimized  351.0  0.475156  0.340359   4.910174e-29   1.573290e-01   4.600303e-01  7.674146e-01  1.0
rho_gini_pav_optimized   585.0  0.259274  0.502604  -1.000000e+00  -3.333333e-02   3.081731e-01  5.970588e-01  1.0
pval_gini_pav_optimized  532.0  0.357559  0.337273   2.646329e-26   4.533488e-02   2.423175e-01  6.354808e-01  1.0
rho_shap_pav_optimized   582.0  0.356543  0.443828  -1.000000e+00   1.193580e-01   3.794236e-01  6.494944e-01  1.0
pval_shap_pav_optimized  534.0  0.365235  0.343619   1.404265e-24   3.531056e-02   2.633922e-01  6.666667e-01  1.0
rho_gini_cnv_optimized   557.0  0.457245  0.467581  -1.000000e+00   2.721805e-01   5.127820e-01  8.000000e-01  1.0
pval_gini_cnv_optimized  487.0  0.291513  0.347116   3.435652e-30   3.242129e-03   1.107872e-01  6.000000e-01  1.0
rho_shap_cnv_optimized   559.0  0.241913  0.512219  -1.000000e+00  -7.197802e-02   2.797555e-01  5.609890e-01  1.0
pval_shap_cnv_optimized  493.0  0.415177  0.336918   7.887333e-14   9.186805e-02   3.819039e-01  6.806381e-01  1.0
'''

res = res.loc[res.Env1 != res.Env2]
''' More summary statistics
>>> res[['rho_shap_snp_optimized', 'pval_shap_snp_optimized']].dropna().describe().T
                         count      mean       std           min       25%       50%       75%       max
rho_shap_snp_optimized   313.0  0.195790  0.378526 -8.000000e-01 -0.017781  0.203297  0.475256  0.942857
pval_shap_snp_optimized  313.0  0.411437  0.303859  4.910174e-29  0.115714  0.396501  0.666667  0.991178
>>> res[['rho_shap_pav_optimized', 'pval_shap_pav_optimized']].dropna().describe().T
                         count      mean       std           min       25%       50%       75%       max
rho_shap_pav_optimized   494.0  0.314794  0.323143 -8.000000e-01  0.113319  0.335347  0.538929  1.000000
pval_shap_pav_optimized  494.0  0.313837  0.303834  1.404265e-24  0.026224  0.208000  0.575114  0.986743
>>> res[['rho_shap_cnv_optimized', 'pval_shap_cnv_optimized']].dropna().describe().T
                         count      mean       std           min       25%       50%       75%       max
rho_shap_cnv_optimized   457.0  0.168993  0.374822 -9.642857e-01 -0.084211  0.202941  0.450688  0.964286
pval_shap_cnv_optimized  457.0  0.369108  0.305527  7.887333e-14  0.064359  0.298089  0.649437  0.988506
'''
