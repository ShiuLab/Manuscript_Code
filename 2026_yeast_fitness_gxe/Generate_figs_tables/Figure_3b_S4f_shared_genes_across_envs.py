#!/usr/bin/env python3
################################################################################
# This script is used to generate:
# 1. Figure 3B: Counts of shared genes across n environments (using SHAP)
# 2. Figure S4F: Counts of shared genes across n environments (using Gini)
# The figures were only made based on the optimized RF model feature importances.
################################################################################

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp
from tqdm import tqdm

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

# Feature to gene mappings; will be used throughout the script
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")


def make_fig3b(df, data_type, imp_type, ax, ax_ridx, ax_cidx):
    """Make Figures 3B & S4F: Distribution of the number of shared genes across n
    number of environments"""
    #
    # Randomize environment columns individually
    n = 10000  # number of iterations to randomize counts
    # empty dictionary to store new env counts for each gene after randomization
    counts_rand = {str(i): None for i in range(n)}
    ks_res = {str(i): None for i in range(n)}  # KS test results
    #
    # original number of environments per gene
    orig_env_counts = df.sum(axis=1)
    #
    for i in tqdm(range(n)):
        df_rand = {}
        for j, env in enumerate(df.columns):
            # randomized sample of if gene is a top gene or not (1/0)
            df_rand[env] = np.random.default_rng(
                seed=i * 100 + j).permutation(df[env])
        #
        df_rand = pd.DataFrame(df_rand)
        # Recalculate number of environments a gene is a top feature for
        env_counts = df_rand.sum(axis=1)
        counts_rand[str(i)] = env_counts  # Save results
        # 2-sample KS test (orig_env_counts cdf is greater than random chance?)
        ks_res[str(i)] = ks_2samp(orig_env_counts,
                                  env_counts, alternative="greater")
        #
    # Plot KS test results
    ks_df = pd.DataFrame.from_dict(ks_res, orient="index")
    f = ks_df.hist(figsize=(8.5, 4), bins=50)
    f[0][1].set_xlim(0, 0.05)
    plt.savefig(
        f"Scripts/Data_Vis/Section_3/Env_shared_genes_{data_type}_{imp_type}_ks2samp_rand_{n}permutations.pdf")
    plt.close()
    ks_df.to_csv(
        f"Scripts/Data_Vis/Section_3/Env_shared_genes_{data_type}_{imp_type}_ks2samp_rand_{n}permutations.csv")
    del f
    #
    # Plot histograms of counts
    # Calculate the median number of environments across all randomizations and the original counts
    # get frequency of env counts at each randomization iteration
    counts_rand_freq = pd.DataFrame(counts_rand).apply(pd.value_counts)
    counts_rand_freq_median = counts_rand_freq.median(
        axis=1).astype("int")  # calculate median counts for each bin
    # KS test between original and median randomized counts
    ks_test = ks_2samp(
        orig_env_counts, counts_rand_freq_median, alternative="greater")
    print(f"KS test between original and median randomized counts: {ks_test}")
    # counted environments per gene
    counted_envs = df.apply(lambda x: " // ".join(x.index), axis=1)
    #
    # Plot the distribution of median counts for each random iteration
    orig_env_counts.plot.hist(bins=35, range=[1, 35], alpha=0.5, color="#009CDC",
                              label=imp_type, ax=ax[ax_ridx][ax_cidx])
    pd.concat([orig_env_counts, counted_envs], ignore_index=False, axis=1).to_csv(
        f"Scripts/Data_Vis/Figure_3f_{data_type}_{imp_type}_data_CORRECTED.csv")
    ax[ax_ridx][ax_cidx].axvline(orig_env_counts.median(
    ), ls="--", color="#FF0000", linewidth=1)  # median env count
    counts_rand_freq_median.plot.line(alpha=0.5, color="#000000",
                                      label="median randomized counts",
                                      ax=ax[ax_ridx][ax_cidx])  # randomized counts
    ax[ax_ridx][ax_cidx].set_title(f"{data_type}: {imp_type}")
    ax[ax_ridx][ax_cidx].legend()
    #
    return ks_df, ax


# Calculate the Gini feature importance rank percentiles
snp_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_snp.tsv", sep="\t", index_col=0)
pav_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_pav.tsv", sep="\t", index_col=0)
cnv_gini = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_cnv.tsv", sep="\t", index_col=0)

snp_gini_env_rank_out = pd.DataFrame()
pav_gini_env_rank_out = pd.DataFrame()
cnv_gini_env_rank_out = pd.DataFrame()
for env in mapping.keys():
    snp_gini_env = snp_gini[[env]].dropna()
    pav_gini_env = pav_gini[[env]].dropna()
    cnv_gini_env = cnv_gini[[env]].dropna()
    #
    # drop features with 0 importance
    snp_gini_env = snp_gini_env[snp_gini_env[env] > 0]
    pav_gini_env = pav_gini_env[pav_gini_env[env] > 0]
    cnv_gini_env = cnv_gini_env[cnv_gini_env[env] > 0]
    #
    # calculate the ranks
    snp_gini_env_rank = snp_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    pav_gini_env_rank = pav_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    cnv_gini_env_rank = cnv_gini_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    #
    # concatenate all env ranks
    snp_gini_env_rank_out = pd.concat(
        [snp_gini_env_rank_out, snp_gini_env_rank], axis=1)
    pav_gini_env_rank_out = pd.concat(
        [pav_gini_env_rank_out, pav_gini_env_rank], axis=1)
    cnv_gini_env_rank_out = pd.concat(
        [cnv_gini_env_rank_out, cnv_gini_env_rank], axis=1)


# Calculate the average absolute SHAP feature importance rank percentiles
snp_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_snp.tsv", sep="\t", index_col=0)  # these are average absolute shap values
pav_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_pav.tsv", sep="\t", index_col=0)
cnv_shap = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_shap_cnv.tsv", sep="\t", index_col=0)

# Calculate the ranks
snp_shap_env_rank_out = pd.DataFrame()
pav_shap_env_rank_out = pd.DataFrame()
cnv_shap_env_rank_out = pd.DataFrame()
for env in mapping.keys():
    snp_shap_env = snp_shap[[env]].dropna()
    pav_shap_env = pav_shap[[env]].dropna()
    cnv_shap_env = cnv_shap[[env]].dropna()
    #
    # drop features with 0 importance
    snp_shap_env = snp_shap_env[snp_shap_env[env] != 0]
    pav_shap_env = pav_shap_env[pav_shap_env[env] != 0]
    cnv_shap_env = cnv_shap_env[cnv_shap_env[env] != 0]
    #
    # calculate the ranks based on the absolute average SHAP values
    snp_shap_env_rank = snp_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    pav_shap_env_rank = pav_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    cnv_shap_env_rank = cnv_shap_env.rank(
        method="average", axis=0, numeric_only=True, ascending=False)
    #
    # concatenate all env ranks
    snp_shap_env_rank_out = pd.concat(
        [snp_shap_env_rank_out, snp_shap_env_rank], axis=1)
    pav_shap_env_rank_out = pd.concat(
        [pav_shap_env_rank_out, pav_shap_env_rank], axis=1)
    cnv_shap_env_rank_out = pd.concat(
        [cnv_shap_env_rank_out, cnv_shap_env_rank], axis=1)


# Map ranked features to genes
snp_gini_rank = map_snps[["snp", "gene"]].merge(
    snp_gini_env_rank_out, left_on="snp", right_index=True, how="right")
pav_gini_rank = map_orfs[["orf", "gene"]].merge(
    pav_gini_env_rank_out, left_on="orf", right_index=True, how="right")
cnv_gini_rank = map_orfs[["orf", "gene"]].merge(
    cnv_gini_env_rank_out, left_on="orf", right_index=True, how="right")

snp_shap_rank = map_snps[["snp", "gene"]].merge(
    snp_shap_env_rank_out, left_on="snp", right_index=True, how="right")
pav_shap_rank = map_orfs[["orf", "gene"]].merge(
    pav_shap_env_rank_out, left_on="orf", right_index=True, how="right")
cnv_shap_rank = map_orfs[["orf", "gene"]].merge(
    cnv_shap_env_rank_out, left_on="orf", right_index=True, how="right")

# Replace gene names with ORFs that aren't in S288C
pav_gini_rank["gene"] = pav_gini_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
cnv_gini_rank["gene"] = cnv_gini_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
pav_shap_rank["gene"] = pav_shap_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
cnv_shap_rank["gene"] = cnv_shap_rank.apply(
    lambda x: x["orf"] if pd.isna(x["gene"]) else x["gene"], axis=1)
pav_gini_rank = pav_gini_rank.drop(columns="orf").set_index("gene")
cnv_gini_rank = cnv_gini_rank.drop(columns="orf").set_index("gene")
pav_shap_rank = pav_shap_rank.drop(columns="orf").set_index("gene")
cnv_shap_rank = cnv_shap_rank.drop(columns="orf").set_index("gene")

# Get the SNP feature with the best rank for each gene
snp_gini_rank_max = snp_gini_rank.groupby("gene").min()  # rank 1 is best
snp_shap_rank_max = snp_shap_rank.groupby("gene").min()
snp_gini_rank_max.drop(index="intergenic", inplace=True)
snp_shap_rank_max.drop(index="intergenic", inplace=True)
snp_gini_rank_max = snp_gini_rank_max.loc[~snp_gini_rank_max.index.str.contains(
    ","), :]  # drop snps that mapped to multiple genes
snp_shap_rank_max = snp_shap_rank_max.loc[~snp_shap_rank_max.index.str.contains(
    ","), :]
snp_gini_rank_max.drop(columns="snp", inplace=True)
snp_shap_rank_max.drop(columns="snp", inplace=True)

# Get the counts of how many environments each gene is in the optimized feature sets for
snp_gini_counts = snp_gini_rank_max.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 2
pav_gini_counts = pav_gini_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 3
cnv_gini_counts = cnv_gini_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 1
snp_shap_counts = snp_shap_rank_max.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 2
pav_shap_counts = pav_shap_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 3
cnv_shap_counts = cnv_shap_rank.applymap(
    lambda x: 1 if x > 0 else 0).sum(axis=1)  # median = 1

# Make Figure 3B and S4F
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(
    12, 7.5), sharex=True)  # env 1-35 x-axis
ks_df, ax = make_fig3b(
    snp_gini_rank_max.applymap(lambda x: 1 if x > 0 else 0), "snp", "gini", ax, 0, 0)
ks_df1, ax = make_fig3b(
    pav_gini_rank.applymap(lambda x: 1 if x > 0 else 0), "pav", "gini", ax, 0, 1)
ks_df2, ax = make_fig3b(
    cnv_gini_rank.applymap(lambda x: 1 if x > 0 else 0), "cnv", "gini", ax, 0, 2)

ks_df3, ax = make_fig3b(
    snp_shap_rank_max.applymap(lambda x: 1 if x > 0 else 0), "snp", "shap", ax, 1, 0)
ks_df4, ax = make_fig3b(
    pav_shap_rank.applymap(lambda x: 1 if x > 0 else 0), "pav", "shap", ax, 1, 1)
ks_df5, ax = make_fig3b(
    cnv_shap_rank.applymap(lambda x: 1 if x > 0 else 0), "cnv", "shap", ax, 1, 2)

plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_3/Figure_3b_S4f_gini_shap_optimized_RF_genes_shared_across_envs_with_rand.pdf")
plt.close('all')

# Outputs:
# snp, gini: KS test between original and median randomized counts: KstestResult(statistic=0.5711833531843341, pvalue=3.902908784435636e-05, statistic_location=13, statistic_sign=1)
# pav, gini: KS test between original and median randomized counts: KstestResult(statistic=0.4444444444444444, pvalue=0.0005590074707047819, statistic_location=18, statistic_sign=1)
# cnv, gini: KS test between original and median randomized counts: KstestResult(statistic=0.5593355481727574, pvalue=6.704828589356284e-05, statistic_location=19, statistic_sign=1)
# snp, shap: KS test between original and median randomized counts: KstestResult(statistic=0.6151385523924894, pvalue=1.5115991629540632e-05, statistic_location=13, statistic_sign=1)
# pav, shap: KS test between original and median randomized counts: KstestResult(statistic=0.4444444444444444, pvalue=0.0005590074707047819, statistic_location=18, statistic_sign=1)
# cnv, shap: KS test between original and median randomized counts: KstestResult(statistic=0.5212403100775194, pvalue=0.00014579537229090375, statistic_location=19, statistic_sign=1)

# Median p-values from the KS tests across all 10,000 randomizations
ks_df = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_snp_gini_ks2samp_rand_10000permutations.csv", index_col=0)
ks_df1 = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_pav_gini_ks2samp_rand_10000permutations.csv", index_col=0)
ks_df2 = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_cnv_gini_ks2samp_rand_10000permutations.csv", index_col=0)
ks_df3 = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_snp_shap_ks2samp_rand_10000permutations.csv", index_col=0)
ks_df4 = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_pav_shap_ks2samp_rand_10000permutations.csv", index_col=0)
ks_df5 = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Env_shared_genes_cnv_shap_ks2samp_rand_10000permutations.csv", index_col=0)

ks_df.pvalue.median()  # 5.944713732943106e-28
ks_df1.pvalue.median()  # 1.524324440450238e-31
ks_df2.pvalue.median()  # 4.707639945685864e-71
ks_df3.pvalue.median()  # 6.140974005168578e-28
ks_df4.pvalue.median()  # 1.524324440450238e-31
ks_df5.pvalue.median()  # 2.0975573145092256e-71
