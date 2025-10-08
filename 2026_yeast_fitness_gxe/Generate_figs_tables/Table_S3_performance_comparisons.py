#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
# TABLE S3
################################################################################
# Correlations between performance of models trained on different data types
# Model results files:
res_base = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_baseline.txt", sep="\t")
res_fs = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_ALL_ALG_FS.txt", sep="\t")

############################ COMPLETE MODELS FIRST #############################
# Compare between data types
pcc_btwn_data = pd.DataFrame()
mwu_btwn_data = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_base.loc[res_base.Alg == alg, ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_data = pd.concat([pcc_btwn_data, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_data = pd.concat([pcc_btwn_data, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")
    # Calculate the Mann-Whitney U test between data types
    for data1 in ["PCs_sklearn", "SNP", "PAV", "CNV"]:
        for data2 in ["PCs_sklearn", "SNP", "PAV", "CNV"]:
            if data1 != data2:
                u_stat, p_value = mannwhitneyu(
                    df[data1], df[data2], alternative="greater")
                mwu_btwn_data = pd.concat([mwu_btwn_data, pd.DataFrame(
                    {f"{alg}.{data1} > {data2}.U": [u_stat], f"{alg}.{data1} > {data2}.p": [p_value]})], axis=1)
                del u_stat, p_value
    del df

pcc_btwn_data.columns = columns
pcc_btwn_data.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_baseline_performance_correlations_between_data_types.tsv", index=True, sep="\t")

mwu_btwn_data = mwu_btwn_data.T
mwu_btwn_data.index = pd.MultiIndex.from_tuples(
    mwu_btwn_data.index.str.split(".", expand=True))
mwu_btwn_data.reset_index(inplace=True)
mwu_btwn_data.columns = ["Algorithm", "Comparison", "Statistic", "Value"]
mwu_btwn_data = mwu_btwn_data.pivot_table(
    index=["Algorithm", "Comparison"], columns="Statistic", values="Value").reset_index()
mwu_btwn_data.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_baseline_performance_mannwhitneyu_between_data_types.tsv", index=False, sep="\t")


# Compare between algorithms
pcc_btwn_alg = pd.DataFrame()
columns = []
for data_type in ["PCs_sklearn", "SNP", "PAV", "CNV"]:
    df = res_base.loc[res_base.Data == data_type, ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_alg = pd.concat([pcc_btwn_alg, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_alg = pd.concat([pcc_btwn_alg, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_alg.columns = columns
pcc_btwn_alg.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_baseline_performance_correlations_between_algorithms.tsv", index=True, sep="\t")

# compare between environments across algorithms
pcc_btwn_env = pd.DataFrame()
columns = []
for data_type in ["PCs_sklearn", "SNP", "PAV", "CNV"]:
    df = res_base.loc[(res_base.Data == data_type),
                      ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_baseline_performance_correlations_between_envs_across_data_types.tsv", index=True, sep="\t")

# compare between environments across data types
pcc_btwn_env = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_base.loc[(res_base.Alg == alg), ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_baseline_performance_correlations_between_envs_across_algorithms.tsv", index=True, sep="\t")

############################### OPTIMIZED MODELS ###############################
# Compare between data types
pcc_btwn_data_fs = pd.DataFrame()
mwu_btwn_data = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_fs.loc[res_fs.Alg == alg, ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_data_fs = pd.concat([pcc_btwn_data_fs, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_data_fs = pd.concat([pcc_btwn_data_fs, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")
    # Calculate the Mann-Whitney U test between data types
    for data1 in ["SNP", "PAV", "CNV"]:
        for data2 in ["SNP", "PAV", "CNV"]:
            if data1 != data2:
                u_stat, p_value = mannwhitneyu(
                    df[data1], df[data2], alternative="greater")
                mwu_btwn_data = pd.concat([mwu_btwn_data, pd.DataFrame(
                    {f"{alg}.{data1} > {data2}.U": [u_stat], f"{alg}.{data1} > {data2}.p": [p_value]})], axis=1)
                del u_stat, p_value

pcc_btwn_data_fs.columns = columns
pcc_btwn_data_fs.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_FS_performance_correlations_between_data_types.tsv", index=True, sep="\t")

mwu_btwn_data = mwu_btwn_data.T
mwu_btwn_data.index = pd.MultiIndex.from_tuples(
    mwu_btwn_data.index.str.split(".", expand=True))
mwu_btwn_data.reset_index(inplace=True)
mwu_btwn_data.columns = ["Algorithm", "Comparison", "Statistic", "Value"]
mwu_btwn_data = mwu_btwn_data.pivot_table(
    index=["Algorithm", "Comparison"], columns="Statistic", values="Value").reset_index()
mwu_btwn_data.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_optimized_performance_mannwhitneyu_between_data_types.tsv", index=False, sep="\t")

# Compare between algorithms
pcc_btwn_alg_fs = pd.DataFrame()
columns = []
for data_type in ["SNP", "PAV", "CNV"]:
    df = res_fs.loc[res_fs.Data == data_type, ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_alg_fs = pd.concat([pcc_btwn_alg_fs, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_alg_fs = pd.concat([pcc_btwn_alg_fs, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_alg_fs.columns = columns
pcc_btwn_alg_fs.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_FS_performance_correlations_between_algorithms.tsv", index=True, sep="\t")

# compare between environments across algorithms
pcc_btwn_env = pd.DataFrame()
columns = []
for data_type in ["SNP", "PAV", "CNV"]:
    df = res_fs.loc[(res_fs.Data == data_type), ["Alg", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Alg", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{data_type}_r")
    columns.append(f"{data_type}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_FS_performance_correlations_between_envs_across_data_types.tsv", index=True, sep="\t")

# compare between environments across data types
pcc_btwn_env = pd.DataFrame()
columns = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    df = res_fs.loc[(res_fs.Alg == alg), ["Data", "Trait", "r2_test"]]
    df = df.pivot(index="Trait", columns="Data", values="r2_test")
    # Calculate pearson correlation r and p-value
    rho = df.T.corr(method=lambda x, y: pearsonr(x, y).statistic)
    pval = df.T.corr(method=lambda x, y: pearsonr(x, y).pvalue)
    pcc_btwn_env = pd.concat([pcc_btwn_env, rho.where(
        np.triu(rho, k=1).astype(bool)).stack()], axis=1)
    pcc_btwn_env = pd.concat([pcc_btwn_env, pval.where(
        np.triu(pval, k=1).astype(bool)).stack()], axis=1)
    columns.append(f"{alg}_r")
    columns.append(f"{alg}_p")

pcc_btwn_env.columns = columns
pcc_btwn_env.to_csv(
    "Scripts/Data_Vis/Section_2/Table_S3_FS_performance_correlations_between_envs_across_algorithms.tsv", index=True, sep="\t")

# Which algorithm results in better performance overall?
rf_greater = []
for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    for alg2 in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
        for data in ["SNP", "PAV", "CNV"]:
            if alg != alg2:
                dist1 = res_fs.loc[(res_fs.Alg == alg) & (res_fs.Data == data), [
                    "new_cond", "r2_test"]].set_index("new_cond")
                dist2 = res_fs.loc[(res_fs.Alg == alg2) & (res_fs.Data == data), [
                    "new_cond", "r2_test"]].set_index("new_cond")
                stat, pval = mannwhitneyu(dist1, dist2, alternative="greater")
                rf_greater.append(
                    [f"{alg} > {alg2}?", data, "FS", stat[0], pval[0]])
                del dist1, dist2, stat, pval


for alg in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
    for alg2 in ["RF", "XGBoost", "rrBLUP", "Bayesian LASSO", "BayesC"]:
        for data in ["PCs_sklearn", "SNP", "PAV", "CNV"]:
            if alg != alg2:
                dist1 = res_base.loc[(res_base.Alg == alg) & (res_base.Data == data), [
                    "new_cond", "r2_test"]].set_index("new_cond")
                dist2 = res_base.loc[(res_base.Alg == alg2) & (res_base.Data == data), [
                    "new_cond", "r2_test"]].set_index("new_cond")
                stat, pval = mannwhitneyu(dist1, dist2, alternative="greater")
                rf_greater.append(
                    [f"{alg} > {alg2}?", data, "Baseline", stat[0], pval[0]])
                del dist1, dist2, stat, pval

rf_greater = pd.DataFrame(
    rf_greater, columns=["Comparison", "Data", "Model Type", "U", "p-value"])

_, q_values, _, _ = multipletests(
    rf_greater["p-value"].values, alpha=0.05, method='fdr_bh')

rf_greater.insert(5, "q-value", q_values)
rf_greater.\
    to_csv("Scripts/Data_Vis/Section_2/Table_S3_RF_vs_other_algs_mannwhitneyu_greater.tsv",
           index=False, sep="\t")
