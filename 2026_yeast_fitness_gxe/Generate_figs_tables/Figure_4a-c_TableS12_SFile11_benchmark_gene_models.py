#!/usr/bin/env python3
'''This script does the following:
1. Create the feature lists for training three types of models:
    a. Only benchmark genes (from the baseline models)
    b. Only important non-benchmark genes (from the optimized models)
    c. Benchmark genes + important non-benchmark genes (combined feature sets)
    Modeling performances are saved as Supplementary file 11.
2. Plot the performances of the three types of models (Fig. 4A-C)
3. Regress model performance on the number of features used to train models (Table S12)
'''

import glob
import os
import re
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, chisquare, ks_2samp, linregress

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
# 1. Feature lists for the "only benchmark genes", "only important non-benchmark
# genes", and "benchmark plus important non-benchmark genes" models
################################################################################
# feature to gene maps
map_snps = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t", index_col=0)
map_orfs.index = map_orfs.index.str.replace(
    "-", ".", regex=False)  # replace . with - in the index names
# remove X from the beginning of the index names
map_orfs.index = map_orfs.apply(lambda x: "X" + x.name, axis=1)
map_orfs.reset_index(inplace=True)
# rename index column to orf
map_orfs.rename(columns={"index": "orf"}, inplace=True)

# benchmark genes to features maps
ben_snp = map_snps.loc[map_snps.Benomyl == 1, ["snp", "gene"]]
ben_orf = map_orfs.loc[map_orfs.Benomyl == 1, ["orf", "gene"]]
caf_snp = map_snps.loc[map_snps.Caffeine == 1, ["snp", "gene"]]
caf_orf = map_orfs.loc[map_orfs.Caffeine == 1, ["orf", "gene"]]
cu_snp = map_snps.loc[map_snps.CuSO4 == 1, ["snp", "gene"]]
cu_orf = map_orfs.loc[map_orfs.CuSO4 == 1, ["orf", "gene"]]
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, ["snp", "gene"]]
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, ["orf", "gene"]]

# feature importance scores from RF models built with the complete feature sets
snp_shap = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_complete_shap_snp.tsv").to_pandas()
pav_shap = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_complete_shap_pav.tsv").to_pandas()
cnv_shap = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_complete_shap_cnv.tsv").to_pandas()
pav_shap = pav_shap.iloc[1:, :]  # remove the first row of NAs
cnv_shap = cnv_shap.iloc[1:, :]
snp_shap.set_index("snp", inplace=True)
pav_shap.set_index("orf", inplace=True)
cnv_shap.set_index("orf", inplace=True)

# feature importance scores from RF models built with the optimized feature sets
snp_shap_opt = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_snp.tsv").to_pandas()
pav_shap_opt = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_pav.tsv").to_pandas()
cnv_shap_opt = dt.fread(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_cnv.tsv").to_pandas()
pav_shap_opt = pav_shap_opt.iloc[1:, :]  # remove the first row of NAs
cnv_shap_opt = cnv_shap_opt.iloc[1:, :]
snp_shap_opt.set_index("snp", inplace=True)
pav_shap_opt.set_index("orf", inplace=True)
cnv_shap_opt.set_index("orf", inplace=True)

# a. Create feature tables that only contain benchmark genes (one variant per gene)
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
snp_bench_list = [caf_snp, caf_snp, ben_snp, cu_snp, sma_snp]
orf_bench_list = [caf_orf, caf_orf, ben_orf, cu_orf, sma_orf]
for i, env in enumerate(target_envs):
    for data_type in ["snp", "pav", "cnv"]:
        # subset the benchmark genes from the feature importance tables
        if data_type == "snp":
            shap_env = snp_shap.loc[
                snp_shap.index.isin(snp_bench_list[i].snp), ["gene", env]]
            shap_opt_env = snp_shap_opt.loc[
                ~snp_shap_opt.index.isin(snp_bench_list[i].snp), ["gene", env]].dropna()  # this still has intergenic snps
        #
        elif data_type == "pav":
            shap_env = pav_shap.loc[
                pav_shap.index.isin(orf_bench_list[i].orf), ["gene", env]]
            shap_opt_env = pav_shap_opt.loc[
                ~pav_shap_opt.index.isin(orf_bench_list[i].orf), ["gene", env]].dropna()
        #
        elif data_type == "cnv":
            shap_env = cnv_shap.loc[cnv_shap.index.isin(
                orf_bench_list[i].orf), ["gene", env]]
            shap_opt_env = cnv_shap_opt.loc[
                ~cnv_shap_opt.index.isin(orf_bench_list[i].orf), ["gene", env]].dropna()
        #
        # Are their missing shap importance values?
        assert shap_env[env].isna().sum(
        ) == 0, f"Missing {data_type} shap importance values for benchmark genes in {env}!"
        assert shap_opt_env[env].isna().sum(
        ) == 0, f"Missing {data_type} shap importance values for important non-benchmark genes in {env}!"
        #
        # Keep one feature per gene
        best_feat = shap_env.groupby("gene")[env].idxmax()
        best_opt_feat = shap_opt_env.groupby("gene")[env].idxmax()
        #
        # Drop intergenic SNPs
        if data_type == "snp":
            best_opt_feat = best_opt_feat[best_opt_feat.index != "intergenic"]
        #
        # Ensure there are no duplicate features (those that mapped to multiple genes)
        assert best_feat.nunique() == len(
            best_feat), f"Duplicate {data_type} features found for {env}!"
        assert best_opt_feat.nunique() == len(
            best_opt_feat), f"Duplicate {data_type} features found in optimized for {env}!"
        #
        # Write the feature lists to files
        best_feat.to_csv(
            f"{path}/Features_complete_only_benchmark_genes_{env}_{data_type}.txt", index=False, header=False)
        best_opt_feat.to_csv(
            f"{path}/Features_optimized_only_important_non_bench_genes_{env}_{data_type}.txt", index=False, header=False)
        #
        combined = pd.concat([best_feat, best_opt_feat], axis=0)
        assert combined.nunique() == len(
            combined), f"Duplicate {data_type} features found in combined for {env}!"
        combined.to_csv(
            f"{path}/Features_bench_plus_important_non_bench_{env}_{data_type}.txt", index=False, header=False)
        #
        del best_feat, best_opt_feat, shap_env, shap_opt_env, combined


# No assertion errors were triggered.

################################################################################
# 2. Plot the performances of the three types of models (Figure 4A-C, S11 File)
################################################################################
# Model performance figure
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF/bench_gene_models"
res = pd.read_csv(f"{d}/RESULTS_reg.txt", sep="\t")
res_lit = res[res.Tag.str.contains("only_benchmark")]
res_lit = res_lit[~(res_lit.Tag.str.startswith("snp") &
                    res_lit.DateTime.str.contains("2025-06-05"))]
res_nlit = res[res.Tag.str.contains("only_important_non_bench")]
res_nlit = res_nlit[~(res_nlit.Tag.str.startswith("snp") &
                      res_nlit.DateTime.str.contains("2025-06-05"))]
res_combined = res[res.Tag.str.contains("plus")]
res_combined = res_combined[~(res_combined.Tag.str.startswith("snp") &
                              res_combined.DateTime.str.contains("2025-06-07"))]

res_combined.insert(2, "Data", res.Tag.str.split(
    "_", expand=True)[0])  # add data column
res_combined.insert(3, "Env", res.Tag.str.split(
    "_", expand=True)[1])  # add env column
res_nlit.insert(2, "Data", res_nlit.Tag.str.split("_", expand=True)[0])
res_nlit.insert(3, "Env", res_nlit.Tag.str.split("_", expand=True)[1])
res_lit.insert(2, "Data", res_lit.Tag.str.split("_", expand=True)[0])
res_lit.insert(3, "Env", res_lit.Tag.str.split("_", expand=True)[1])

# Plot performances of all models
fig, ax = plt.subplots(nrows=3, ncols=3, sharex=True,
                       sharey=True, figsize=(7, 8))
for i, data_type in enumerate(["snp", "pav", "cnv"]):
    # Figures for important non-benchmark genes + benchmark genes
    sns.barplot(res_combined[res_combined.Data == data_type].sort_values(by="Tag"),
                x="Y", y="r2_test", hue="Y", ax=ax[0][i])
    ax[0][i].errorbar(
        x=res_combined.loc[res_combined.Data == data_type].sort_values(by="Tag")[
            "Y"],
        y=res_combined.loc[res_combined.Data == data_type].sort_values(by="Tag")[
            "r2_test"],
        yerr=res_combined.loc[res_combined.Data == data_type].sort_values(by="Tag")[
            "r2_test_sd"],
        fmt='o', color='black', capsize=5)
    ax[0][i].set_title(
        f"{data_type} RF models bench + important non-bench genes")
    ax[0][i].set_xticklabels(
        res_combined[res_combined.Data == data_type].sort_values(by="Tag").Y,
        size=7, rotation=50, ha="right")
    #
    # Figures for only important non-benchmark genes from the optimized RF models
    sns.barplot(res_nlit[res_nlit.Data == data_type].sort_values(by="Tag"),
                x="Y", y="r2_test", hue="Y", ax=ax[1][i])
    ax[1][i].errorbar(
        x=res_nlit.loc[res_nlit.Data == data_type].sort_values(by="Tag")["Y"],
        y=res_nlit.loc[res_nlit.Data == data_type].sort_values(by="Tag")[
            "r2_test"],
        yerr=res_nlit.loc[res_nlit.Data == data_type].sort_values(by="Tag")[
            "r2_test_sd"],
        fmt='o', color='black', capsize=5)
    ax[1][i].set_title(f"{data_type} RF important non-bench genes only")
    ax[1][i].set_xticklabels(res_nlit[res_nlit.Data == data_type].sort_values(by="Tag").Y,
                             size=7, rotation=50, ha="right")
    #
    # Figures for benchmark genes only from the baseline models
    sns.barplot(res_lit[res_lit.Data == data_type].sort_values(by="Tag"),
                x="Y", y="r2_test", hue="Y", ax=ax[2][i])
    ax[2][i].errorbar(
        x=res_lit.loc[res_lit.Data == data_type].sort_values(by="Tag")["Y"],
        y=res_lit.loc[res_lit.Data == data_type].sort_values(by="Tag")[
            "r2_test"],
        yerr=res_lit.loc[res_lit.Data == data_type].sort_values(by="Tag")[
            "r2_test_sd"],
        fmt='o', color='black', capsize=5)
    ax[2][i].set_title(f"{data_type}: only bench genes from complete RF")
    ax[2][i].set_xticklabels(res_lit[res_lit.Data == data_type].sort_values(by="Tag").Y,
                             size=7, rotation=50, ha="right")


plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_4/Figure_4a_benchmark_gene_model_performances_not_RF_FS.pdf")
plt.close("all")

# What is the difference in model performances between the three types of models?
'''>>> res_nlit.reset_index().sort_values(by=["Data", "Env"]).reset_index().r2_test - \
... res_lit.reset_index().sort_values(by=["Data", "Env"]).reset_index().r2_test
0.311960 <-- cnv, YPDBENOMYL500
0.186067 <-- cnv, YPDCAFEIN40
0.179139 <-- cnv, YPDCAFEIN50
0.418395 <-- cnv, YPDCUSO410MM
0.070128 <-- cnv, YPDSODIUMMETAARSENITE
0.620563 <-- pav, YPDBENOMYL500
0.461988 <-- pav, YPDCAFEIN40
0.505226 <-- pav, YPDCAFEIN50
0.298687 <-- pav, YPDCUSO410MM
0.232443 <-- pav, YPDSODIUMMETAARSENITE
0.006514 <-- snp, YPDBENOMYL500
0.015958 <-- snp, YPDCAFEIN40
-0.000193 <-- snp, YPDCAFEIN50
0.052007 <-- snp, YPDCUSO410MM
0.036731 <-- snp, YPDSODIUMMETAARSENITE
>>> res_combined.reset_index().sort_values(by=["Data", "Env"]).reset_index().r2_test - \
... res_lit.reset_index().sort_values(by=["Data", "Env"]).reset_index().r2_test
0.310364  <-- cnv, YPDBENOMYL500
0.213597  <-- cnv, YPDCAFEIN40
0.205563  <-- cnv, YPDCAFEIN50
0.423909  <-- cnv, YPDCUSO410MM
0.080054  <-- cnv, YPDSODIUMMETAARSENITE
0.608730  <-- pav, YPDBENOMYL500
0.481588  <-- pav, YPDCAFEIN40
0.509258  <-- pav, YPDCAFEIN50
0.289978  <-- pav, YPDCUSO410MM
0.226765  <-- pav, YPDSODIUMMETAARSENITE
0.009382  <-- snp, YPDBENOMYL500
0.016142  <-- snp, YPDCAFEIN40
0.014707  <-- snp, YPDCAFEIN50
0.055453  <-- snp, YPDCUSO410MM
0.045679  <-- snp, YPDSODIUMMETAARSENITE
'''

################################################################################
# 3. Regress model performance on the number of features used to train models (Table S12)
################################################################################
linreg_res = {'benchmark gene models': {},
              'important non-benchmark gene models': {},
              'combined models': {}}
for data_type in ["snp", "pav", "cnv"]:
    m, b, rval, pval, se = linregress(
        res_lit[res_lit.Data == data_type]["FeatureNum"],
        res_lit[res_lit.Data == data_type]["r2_test"],
        alternative="two-sided")
    linreg_res['benchmark gene models'][data_type] = {
        "slope": m, "intercept": b, "r": rval, "p-value": pval, "se": se}
    m, b, rval, pval, se = linregress(
        res_nlit[res_nlit.Data == data_type]["FeatureNum"],
        res_nlit[res_nlit.Data == data_type]["r2_test"],
        alternative="two-sided")
    linreg_res['important non-benchmark gene models'][data_type] = {
        "slope": m, "intercept": b, "r": rval, "p-value": pval, "se": se}
    m, b, rval, pval, se = linregress(
        res_combined[res_combined.Data == data_type]["FeatureNum"],
        res_combined[res_combined.Data == data_type]["r2_test"],
        alternative="two-sided")
    linreg_res['combined models'][data_type] = {
        "slope": m, "intercept": b, "r": rval, "p-value": pval, "se": se}

linreg_res_df = pd.json_normalize(linreg_res).transpose()
linreg_res_df.index = pd.MultiIndex.from_tuples(linreg_res_df.index.str.split(".").map(tuple),
                                                names=["Model Type", "Data", "Metric"])
out = linreg_res_df.pivot_table(index=["Model Type", "Data"], columns="Metric")

# combine with the results if all the models are included instead of separated by data type
m, b, r, p, se = linregress(
    res_lit.FeatureNum, res_lit.r2_test, alternative="two-sided")
m1, b1, r1, p1, se1 = linregress(
    res_nlit.FeatureNum, res_nlit.r2_test, alternative="two-sided")
m2, b2, r2, p2, se2 = linregress(
    res_combined.FeatureNum, res_combined.r2_test, alternative="two-sided")

out.reset_index(inplace=True)
out.columns = out.columns.droplevel(0)
out.columns = ["Model Type", "Data",
               "intercept", "p-value", "r", "se", "slope"]
pd.concat([out, pd.DataFrame.from_dict({
    "benchmark gene models": {"Data": "all", "intercept": b, "p-value": p,
                              "r": r, "se": se, "slope": m},
    "important non-benchmark gene models": {"Data": "all", "intercept": b1, "p-value": p1,
                                            "r": r1, "se": se1, "slope": m1},
    "combined models": {"Data": "all", "intercept": b2, "p-value": p2,
                        "r": r2, "se": se2, "slope": m2}}).transpose().reset_index()],
          axis=0, ignore_index=True).to_csv(
              "Scripts/Data_Vis/Section_4/Table_S12_benchmark_gene_model_performances_vs_feature_num.csv")

linregress(res_lit.FeatureNum, res_lit.r2_test)
# LinregressResult(slope=0.00017299985470615398, intercept=0.2710279232187656, rvalue=0.18788040604989847, pvalue=0.4276454100569397, stderr=0.00021316914804781604, intercept_stderr=0.1055754074416109)
linregress(res_nlit.FeatureNum, res_nlit.r2_test)
# LinregressResult(slope=4.516410514222374e-05, intercept=0.501827342962556, rvalue=0.32488150775165187, pvalue=0.23740645506868835, stderr=3.646491980499188e-05, intercept_stderr=0.026167808447933198)
linregress(res_combined.FeatureNum, res_combined.r2_test)
# LinregressResult(slope=4.764530204557169e-05, intercept=0.486583257565177, rvalue=0.3702542760751369, pvalue=0.174318524382342, stderr=3.3153663261356074e-05, intercept_stderr=0.034807811959239855)
