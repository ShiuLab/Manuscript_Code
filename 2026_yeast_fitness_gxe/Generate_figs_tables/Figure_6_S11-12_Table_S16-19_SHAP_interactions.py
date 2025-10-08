#! /usr/bin/env python3
"""Analysis for the section on genetic interactions underlying benomyl stress in
yeast. This script does the following:
1. Generate feature tables for SNP, PAV, CNV, SNP+PAV, SNP+CNV, PAV+CNV, and
   SNP+PAV+CNV benomyl benchmark gene RF models. (Model performances in Fig. 6A)
2. Plot the performances (Table S16) of the best RF benomyl models out of all
   the training repetitions. These models will be used to generate SHAP
   interaction scores (Table S17).
3. Exploratory analysis of SHAP interaction scores (Table S17):
   - Determine the number of unique variant-variant and gene-gene interactions
     (Fig. S12: UpSet plot of the unique gene-gene interactions identified by
     the SHAP interactions from different models)
   - Determine what variant-variant types identified the same gene-gene interactions
     (Table S18)
4. Enrichment analysis of experimentally verified GIs among the SHAP interactions
   - Supplementary file 10: Experimentally verified genetic interactions from
     BioGRID & Costanzo et al. 2021
   - Figure S11: Venn diagram of overlap between experimentally verified GIs
   - Table S19: Enrichment results
5. SHAP interaction plots for the top SHAP interactions from each model
   - Figures 6B-F
"""

import pickletools
import joblib
import os
import re
import datatable as dt
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import numpy as np
from sklearn.metrics import r2_score
from statannotations.Annotator import Annotator
from itertools import combinations
from scipy.stats import mannwhitneyu
from collections import defaultdict
from upsetplot import from_memberships
from upsetplot import UpSet
from venn import venn
from tqdm import tqdm
from scipy.stats import false_discovery_control
from scipy.stats import fisher_exact

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
""" 1. Generate feature tables for SNP, PAV, CNV, SNP+PAV, SNP+CNV, PAV+CNV, and
       SNP+PAV+CNV benomyl benchmark gene RF models. Each gene will be
       represented by one feature."""
################################################################################
# Path to save feature tables to
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_rf"

# Benomyl benchmark genes
ben_genes = pd.read_csv(
    "Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")

# Generate sets of benomyl benchmark gene SNP, PAV, and CNV features
for data_type in ["snp", "pav", "cnv"]:
    # Read in average absolute SHAP values from RF models trained on complete feature sets
    df = dt.fread(
        f"Scripts/Data_Vis/Section_3/RF_complete_shap_{data_type}.tsv").to_pandas()
    print(df.head())
    #
    if data_type == "snp":
        df.set_index(["snp", "gene"], inplace=True)
        # remove "gene_with_intergenic" column; ranking will not include intergenic snps
        df = df.iloc[:, 1:]
        df = df.loc[df.index.get_level_values(
            "gene") != "intergenic", :]  # drop intergenic snps
    else:
        df = df.loc[df.orf != "", :]  # drop the extra row with no orf name
        df.set_index(["orf", "gene"], inplace=True)
    #
    # Get all variants per gene
    if data_type == "snp":
        ben_feat = df.loc[df.index.get_level_values("gene").
                          isin(ben_genes['Gene Systematic Name']), :].\
            index.get_level_values("snp")
    else:
        ben_feat = df.loc[df.index.get_level_values("gene").
                          isin(ben_genes['Gene Systematic Name']), :].\
            index.get_level_values("orf")
    ben_feat = pd.Series(ben_feat)
    assert ben_feat.nunique() == len(ben_feat), "Duplicate features found"
    #
    # ben_feat.to_csv(
    #     f"{d}/Features_all_variants_per_gene_benomyl_500ugml_{data_type}.txt",
    #     index=False, header=False)
    #
    # Get only the most important variant per benomyl benchmark gene
    # Take the max average absolute SHAP value per gene
    if data_type == "snp":
        env_imp = df.loc[df.index.get_level_values("snp").isin(ben_feat),
                         "YPDBENOMYL500"].dropna().sort_values(ascending=False)
        # remove features with 0 gini importance
        # env_imp = env_imp.loc[env_imp != 0.0, :] (need representation of all benchmark genes)
        max_features = pd.concat([
            env_imp.groupby("gene").idxmax().apply(lambda x: x[0]),
            env_imp.groupby("gene").max()], ignore_index=False, axis=1)
        max_features.columns = ["snp", "max_imp"]
        max_features["snp"].to_csv(
            f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp.txt",
            sep="\t", header=None, index=False)
    else:
        env_imp = df.loc[df.index.get_level_values("orf").isin(ben_feat),
                         "YPDBENOMYL500"].dropna().sort_values(ascending=False)
        max_features = pd.concat([
            env_imp.groupby("gene").idxmax().apply(lambda x: x[0]),
            env_imp.groupby("gene").max()], ignore_index=False, axis=1)
        max_features.columns = ["orf", "max_imp"]
        max_features["orf"].to_csv(
            f"{path}/Features_one_variant_per_gene_benomyl_500ugml_{data_type}.txt",
            sep="\t", header=None, index=False)
    del df, ben_feat, env_imp, max_features


# Read in the benomyl benchmark gene SNP, PAV, CNV feature tables
snp_1feat = pd.read_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp.txt",
    sep="\t", header=None, names=["Feature"])
pav_1feat = pd.read_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_pav.txt",
    sep="\t", header=None, names=["Feature"])
cnv_1feat = pd.read_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_cnv.txt",
    sep="\t", header=None, names=["Feature"])

# pav_1feat["Feature"] = pav_1feat.apply(lambda x: "X" + x.Feature, axis=1)
# pav_1feat["Feature"] = pav_1feat.apply(
#     lambda x: re.sub("-", ".", x.Feature), axis=1)
# cnv_1feat["Feature"] = cnv_1feat.apply(lambda x: "X" + x.Feature, axis=1)
# cnv_1feat["Feature"] = cnv_1feat.apply(
#     lambda x: re.sub("-", ".", x.Feature), axis=1)
# pav_1feat.to_csv(
#     f"{path}/Features_one_variant_per_gene_benomyl_500ugml_pav_fixed.txt",
#     sep="\t", header=None, index=False)
# cnv_1feat.to_csv(
#     f"{path}/Features_one_variant_per_gene_benomyl_500ugml_cnv_fixed.txt",
#     sep="\t", header=None, index=False)

# PAVs and CNVs have the same ORFs
set(pav_1feat["Feature"]) == set(cnv_1feat["Feature"])  # True

# Read in the SNP, PAV, and CNV feature tables
snp = dt.fread(
    "Data/Peter_2018/geno_corrected.csv").to_pandas().set_index("ID")
pav = dt.fread(
    "Data/Peter_2018/ORFs_pres_abs.csv").to_pandas().set_index("ID").astype(int)
cnv = dt.fread(
    "Data/Peter_2018/ORFs_no_NA.csv").to_pandas().set_index("ID")


def fix_orf_ids(df):
    """ Fix the ORF IDs in the dataframe by removing the prefix 'X' and
    replacing '.' with '-' """
    df.columns = df.apply(lambda x: re.sub(
        "^X", "", x.name), axis=0)  # fix orf IDs
    df.columns = df.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    return df


pav = fix_orf_ids(pav)
cnv = fix_orf_ids(cnv)

pav = pav.loc[snp.index, :]  # make sure the rows are in the same order
cnv = cnv.loc[snp.index, :]

# integrate the features and keep record from which dataset the ORFs are coming from
integrated_spc = pd.concat([
    snp.loc[:, snp_1feat.Feature],
    pav.loc[:, pav_1feat.Feature].rename(columns=lambda x: f"{x};PAV"),
    cnv.loc[:, cnv_1feat.Feature].rename(columns=lambda x: f"{x};CNV")],
    axis=1, ignore_index=False)
integrated_spc.to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav_cnv.txt")

integrated_sp = pd.concat([
    snp.loc[:, snp_1feat.Feature],
    pav.loc[:, pav_1feat.Feature].rename(columns=lambda x: f"{x};PAV")],
    axis=1, ignore_index=False)
integrated_sp.to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav.txt")

integrated_sc = pd.concat([
    snp.loc[:, snp_1feat.Feature],
    cnv.loc[:, cnv_1feat.Feature].rename(columns=lambda x: f"{x};CNV")],
    axis=1, ignore_index=False)
integrated_sc.to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_snp_cnv.txt")

integrated_pc = pd.concat([
    pav.loc[:, pav_1feat.Feature].rename(columns=lambda x: f"{x};PAV"),
    cnv.loc[:, cnv_1feat.Feature].rename(columns=lambda x: f"{x};CNV")],
    axis=1, ignore_index=False)
integrated_pc.to_csv(
    f"{path}/Features_one_variant_per_gene_benomyl_500ugml_pav_cnv.txt")

################################################################################
# 2. Determine the best RF benomyl models out of all training repetitions, plot
#    their performances. These models will be used to generate SHAP interaction
#    scores.
################################################################################
# Read in the model performance scores
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_rf"

snp_scores = pd.read_csv(
    f"{path}/snp_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
pav_scores = pd.read_csv(
    f"{path}/pav_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
cnv_scores = pd.read_csv(
    f"{path}/cnv_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
snp_pav_cnv_scores = pd.read_csv(
    f"{path}/integrated_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
snp_pav_scores = pd.read_csv(
    f"{path}/integrated_snp_pav_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
snp_cnv_scores = pd.read_csv(
    f"{path}/integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)
pav_cnv_scores = pd.read_csv(
    f"{path}/integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml_scores.txt", sep="\t", index_col=0)

# get the best model based on the validation R2 score
test = pd.read_csv(
    "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/Test.txt", sep="\t", header=None)

save_names = ["snp_one_variant_per_gene_benomyl_500ugml",
              "pav_one_variant_per_gene_benomyl_500ugml",
              "cnv_one_variant_per_gene_benomyl_500ugml",
              "integrated_one_variant_per_gene_benomyl_500ugml",
              "integrated_snp_pav_one_variant_per_gene_benomyl_500ugml",
              "integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml",
              "integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml"]

with open(f"{path}/best_benomyl_models.txt", "w") as f:
    for i, scores in enumerate([snp_scores, pav_scores, cnv_scores,
                               snp_pav_cnv_scores, snp_pav_scores,
                               snp_cnv_scores, pav_cnv_scores]):
        train_scores = scores.loc[~scores.index.isin(test[0]), :]
        best_rep = train_scores.drop(columns=["Mean", "stdev"]).apply(
            lambda c: r2_score(train_scores["Y"], train_scores[c.name]
                               if c.name != "Y" else c), axis=0).iloc[1:].idxmax()
        f.write(
            f"{save_names[i]}_models_rep_{int(best_rep.split('_')[1]) - 1}.pkl\n")

### Figure 6A and Table S16. Model performances -------------------------------#
rf = pd.read_csv(f"{path}/RESULTS_reg.txt", sep="\t")

# for plotting purposes, set x-axis labels
rf["Model"] = ["SNP", "SNP + PAV + CNV", "SNP + PAV",
               "SNP + CNV", "PAV", "CNV", "PAV + CNV"]
rf.to_excel(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S16_benomyl_model_performances_RF.xlsx",
    index=False)

ax = sns.barplot(data=rf, x="Model", y="r2_test")
plt.errorbar(data=rf, x="Model", y="r2_test",
             yerr="r2_test_sd", fmt="o", color="black", capsize=5)
plt.ylabel("RF R2 test")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_6/shap_interaction/Figure_6a_benomyl_model_performances_RF.pdf")
plt.close()
### ---------------------------------------------------------------------------#

################################################################################
# 3. Exploratory analysis of SHAP interaction scores
################################################################################
os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# SHAP interaction scores from SNP, PAV, CNV, and Integrated models
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_rf"
snp_shap_int = pd.read_csv(
    f"{path}/snp/shap_interaction_scores_snp_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
pav_shap_int = pd.read_csv(
    f"{path}/pav/shap_interaction_scores_pav_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
cnv_shap_int = pd.read_csv(
    f"{path}/cnv/shap_interaction_scores_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
spc_shap_int = pd.read_csv(
    f"{path}/snp_pav_cnv/shap_interaction_scores_integrated_snp_pav_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
sp_shap_int = pd.read_csv(
    f"{path}/snp_pav/shap_interaction_scores_integrated_snp_pav_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
sc_shap_int = pd.read_csv(
    f"{path}/snp_cnv/shap_interaction_scores_integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")
pc_shap_int = pd.read_csv(
    f"{path}/pav_cnv/shap_interaction_scores_integrated_pav_cnv_one_variant_per_gene_benomyl_500ugml_summed.txt", sep="\t")

# Add a column to indicate which model the SHAP interaction scores come from
snp_shap_int.insert(3, "Model", "SNP")
pav_shap_int.insert(3, "Model", "PAV")
cnv_shap_int.insert(3, "Model", "CNV")
spc_shap_int.insert(3, "Model", "SNP + PAV + CNV")
sp_shap_int.insert(3, "Model", "SNP + PAV")
sc_shap_int.insert(3, "Model", "SNP + CNV")
pc_shap_int.insert(3, "Model", "PAV + CNV")

### SHAP interaction score violin plots (PROBABLY DELETE) -------------------------------#
violin_data = pd.concat([snp_shap_int, pav_shap_int, cnv_shap_int,
                         spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int], axis=0,
                        ignore_index=True)

models = violin_data["Model"].unique()
pairs = list(combinations(models, 2))  # all pairwise combinations

ax = sns.violinplot(data=violin_data, inner="box", fill=False, x="Model", log_scale=10,
                    y="Interaction", hue="Model")
plt.axhline(violin_data["Interaction"].quantile(
    0.95), color="red", linestyle="--")
plt.axhline(violin_data["Interaction"].quantile(
    0.99), color="blue", linestyle="--")
annotator = Annotator(ax, pairs, data=violin_data,
                      x="Model", y="Interaction", order=models)
annotator.configure(test="Mann-Whitney", text_format="star",
                    loc="outside", comparisons_correction="bonferroni")
annotator.apply_and_annotate()
'''Annotator results:
p-value annotation legend:
      ns: 5.00e-02 < p <= 1.00e+00
       *: 1.00e-02 < p <= 5.00e-02
      **: 1.00e-03 < p <= 1.00e-02
     ***: 1.00e-04 < p <= 1.00e-03
    ****: p <= 1.00e-04

SNP vs. PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:2.030e-56 U_stat=1.332e+07
PAV vs. CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.000e+00 U_stat=3.293e+06
CNV vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:2.277e-188 U_stat=1.256e+09
SNP + PAV + CNV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.270e+09
SNP + PAV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.141e+09
SNP + CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=5.632e+08
SNP vs. CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=4.354e+08
PAV vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:5.093e-06 U_stat=3.933e+07
CNV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.555e+08
SNP + PAV + CNV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=2.241e+09
SNP + PAV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=4.959e+08
SNP vs. SNP + PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=4.952e+09
PAV vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:6.959e-49 U_stat=5.262e+06
CNV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=2.614e+08
SNP + PAV + CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.278e-71 U_stat=1.291e+09
SNP vs. SNP + PAV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:4.568e-15 U_stat=8.202e+08
PAV vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.889e-16 U_stat=8.532e+06
CNV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.017e-271 U_stat=1.395e+08
SNP vs. SNP + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=1.230e+09
PAV vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:1.088e-13 U_stat=4.327e+06
SNP vs. PAV + CNV: Mann-Whitney-Wilcoxon test two-sided with Bonferroni correction, P_val:0.000e+00 U_stat=5.289e+08
'''
plt.ylabel("SHAP interaction score")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_6/shap_interaction/violin_log10_summed_shap_interaction.pdf")
plt.close()

rows_above_95th = violin_data[violin_data["Interaction"]
                              > violin_data["Interaction"].quantile(0.95)]
rows_above_99th = violin_data[violin_data["Interaction"]
                              > violin_data["Interaction"].quantile(0.99)]
rows_above_95th["Model"].value_counts()
# Model
# SNP                4490
# SNP + PAV          4236
# SNP + CNV          4127
# SNP + PAV + CNV    1857
# CNV                 533
# PAV + CNV           274
# PAV                  24
rows_above_99th["Model"].value_counts()
# SNP + CNV          985
# SNP + PAV          889
# SNP                786
# SNP + PAV + CNV    213
# CNV                150
# PAV + CNV           80
# PAV                  6

### ---------------------------------------------------------------------------#
"""To determine the number of unique genetic interactions, there are different
levels to consider: variant-to-variant(might have duplicate interactions at the
gene level) and gene-to-gene (this is the one we care about most)."""

# Add the feature data columns to the SHAP interaction score dataframes
snp_shap_int.insert(4, "Feature1_Data", "SNP")
snp_shap_int.insert(5, "Feature2_Data", "SNP")
pav_shap_int.insert(4, "Feature1_Data", "PAV")
pav_shap_int.insert(5, "Feature2_Data", "PAV")
cnv_shap_int.insert(4, "Feature1_Data", "CNV")
cnv_shap_int.insert(5, "Feature2_Data", "CNV")
spc_shap_int.insert(4, "Feature1_Data", spc_shap_int.apply(
    lambda x: x["Feature1"].split(";")[1] if ";" in x["Feature1"] else "SNP",
    axis=1))
spc_shap_int.insert(5, "Feature2_Data", spc_shap_int.apply(
    lambda x: x["Feature2"].split(";")[1] if ";" in x["Feature2"] else "SNP",
    axis=1))
sp_shap_int.insert(4, "Feature1_Data", sp_shap_int.apply(
    lambda x: x["Feature1"].split(";")[1] if ";" in x["Feature1"] else "SNP",
    axis=1))
sp_shap_int.insert(5, "Feature2_Data", sp_shap_int.apply(
    lambda x: x["Feature2"].split(";")[1] if ";" in x["Feature2"] else "SNP",
    axis=1))
sc_shap_int.insert(4, "Feature1_Data", sc_shap_int.apply(
    lambda x: x["Feature1"].split(";")[1] if ";" in x["Feature1"] else "SNP",
    axis=1))
sc_shap_int.insert(5, "Feature2_Data", sc_shap_int.apply(
    lambda x: x["Feature2"].split(";")[1] if ";" in x["Feature2"] else "SNP",
    axis=1))
pc_shap_int.insert(4, "Feature1_Data", pc_shap_int.apply(
    lambda x: x["Feature1"].split(";")[1], axis=1))
pc_shap_int.insert(5, "Feature2_Data", pc_shap_int.apply(
    lambda x: x["Feature2"].split(";")[1], axis=1))


def fix_orf_names(orf_pd_series, axis=0):
    '''Prepare ORF identifiers for mapping to genes:'''
    orf_pd_series = orf_pd_series.apply(lambda x: re.sub("^X", "", x))
    orf_pd_series = orf_pd_series.apply(lambda x: re.sub(";PAV$", "", x))
    orf_pd_series = orf_pd_series.apply(lambda x: re.sub(";CNV$", "", x))
    # orf_pd_series = orf_pd_series.apply(lambda x: re.sub("\.0$", "", x)) # remove the .0 from the end of the column names
    orf_pd_series = orf_pd_series.apply(lambda x: re.sub(
        "\.", "-", x))  # replace . with - in the column names
    return orf_pd_series


pav_shap_int.Feature1 = fix_orf_names(pav_shap_int.Feature1)
pav_shap_int.Feature2 = fix_orf_names(pav_shap_int.Feature2)
cnv_shap_int.Feature1 = fix_orf_names(cnv_shap_int.Feature1)
cnv_shap_int.Feature2 = fix_orf_names(cnv_shap_int.Feature2)
spc_shap_int.Feature1 = fix_orf_names(spc_shap_int.Feature1)
spc_shap_int.Feature2 = fix_orf_names(spc_shap_int.Feature2)
sp_shap_int.Feature1 = fix_orf_names(sp_shap_int.Feature1)
sp_shap_int.Feature2 = fix_orf_names(sp_shap_int.Feature2)
sc_shap_int.Feature1 = fix_orf_names(sc_shap_int.Feature1)
sc_shap_int.Feature2 = fix_orf_names(sc_shap_int.Feature2)
pc_shap_int.Feature1 = fix_orf_names(pc_shap_int.Feature1)
pc_shap_int.Feature2 = fix_orf_names(pc_shap_int.Feature2)


def feature2gene(feature):
    # Map features to genes
    try:
        return map_orfs_dict["gene"][feature]
    except KeyError:
        # try:
        return map_snps_dict["gene"][feature]


# SNP and ORF gene maps
map_snps = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")
map_snps_dict = map_snps[["snp", "gene"]].set_index(
    "snp").to_dict(orient="dict")
map_orfs_dict = map_orfs[["orf", "gene"]].set_index(
    "orf").to_dict(orient="dict")


snp_shap_int.insert(0, "Gene2", snp_shap_int.apply(
    lambda x: map_snps_dict["gene"][x["Feature2"]], axis=1))
snp_shap_int.insert(0, "Gene1", snp_shap_int.apply(
    lambda x: map_snps_dict["gene"][x["Feature1"]], axis=1))
pav_shap_int.insert(0, "Gene2", pav_shap_int.apply(
    lambda x: map_orfs_dict["gene"][x["Feature2"]], axis=1))
pav_shap_int.insert(0, "Gene1", pav_shap_int.apply(
    lambda x: map_orfs_dict["gene"][x["Feature1"]], axis=1))
cnv_shap_int.insert(0, "Gene2", cnv_shap_int.apply(
    lambda x: map_orfs_dict["gene"][x["Feature2"]], axis=1))
cnv_shap_int.insert(0, "Gene1", cnv_shap_int.apply(
    lambda x: map_orfs_dict["gene"][x["Feature1"]], axis=1))
spc_shap_int.insert(0, "Gene2", spc_shap_int.apply(
    lambda x: feature2gene(x["Feature2"]), axis=1))
spc_shap_int.insert(0, "Gene1", spc_shap_int.apply(
    lambda x: feature2gene(x["Feature1"]), axis=1))
sp_shap_int.insert(0, "Gene2", sp_shap_int.apply(
    lambda x: feature2gene(x["Feature2"]), axis=1))
sp_shap_int.insert(0, "Gene1", sp_shap_int.apply(
    lambda x: feature2gene(x["Feature1"]), axis=1))
sc_shap_int.insert(0, "Gene2", sc_shap_int.apply(
    lambda x: feature2gene(x["Feature2"]), axis=1))
sc_shap_int.insert(0, "Gene1", sc_shap_int.apply(
    lambda x: feature2gene(x["Feature1"]), axis=1))
pc_shap_int.insert(0, "Gene2", pc_shap_int.apply(
    lambda x: feature2gene(x["Feature2"]), axis=1))
pc_shap_int.insert(0, "Gene1", pc_shap_int.apply(
    lambda x: feature2gene(x["Feature1"]), axis=1))

# check for missing gene identifiers
for df in [snp_shap_int, pav_shap_int, cnv_shap_int, spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int]:
    # nothing was missing, no problems with the gene names.
    print(df.loc[spc_shap_int.isna().any(axis=1), :].shape)


def get_unique_gp(df, Gene1="Gene1", Gene2="Gene2"):
    '''
    Get the unique gene pairs from the dataframe
    '''
    df_gp = df.apply(lambda x: set(
        [x[Gene1], x[Gene2]]), axis=1).values  # gene pairs
    df_gp = {frozenset(sorted(set))
             for set in df_gp}  # get unique interactions
    return df_gp


# get the unique gene pairs
snp_shap_int_gp = get_unique_gp(snp_shap_int, "Gene1", "Gene2")
pav_shap_int_gp = get_unique_gp(pav_shap_int, "Gene1", "Gene2")
cnv_shap_int_gp = get_unique_gp(cnv_shap_int, "Gene1", "Gene2")
len(snp_shap_int_gp)  # 40780
len(pav_shap_int_gp)  # 455 (do not consider this model bc of poor performance)
len(cnv_shap_int_gp)  # 14430
len(snp_shap_int_gp.union(pav_shap_int_gp).union(
    cnv_shap_int_gp))  # 47,006 unique gene-gene interactions

spc_shap_int_gp = get_unique_gp(spc_shap_int, "Gene1", "Gene2")
sp_shap_int_gp = get_unique_gp(sp_shap_int, "Gene1", "Gene2")
sc_shap_int_gp = get_unique_gp(sc_shap_int, "Gene1", "Gene2")
pc_shap_int_gp = get_unique_gp(pc_shap_int, "Gene1", "Gene2")
len(spc_shap_int_gp)  # 64156
len(sp_shap_int_gp)  # 37573
len(sc_shap_int_gp)  # 36490
len(pc_shap_int_gp)  # 14280
len(spc_shap_int_gp.union(sp_shap_int_gp).union(
    sc_shap_int_gp).union(pc_shap_int_gp))  # 68420
len(snp_shap_int_gp.union(cnv_shap_int_gp).union(spc_shap_int_gp).
    union(sp_shap_int_gp).union(sc_shap_int_gp).union(pc_shap_int_gp))  # 69486 total unique gene-gene interactions (PAV model excl.)

tableS17_data = pd.concat([snp_shap_int, pav_shap_int, cnv_shap_int,
                           spc_shap_int, sp_shap_int, sc_shap_int, pc_shap_int], axis=0,
                          ignore_index=True)

len(get_unique_gp(tableS17_data, "Gene1", "Gene2"))  # 69495 (w/PAV model incl.)

# Save all interactions to a file
tableS17_data.isna().sum()  # no missing values!
tableS17_data.to_csv(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_benchmark_RF_models_SHAP_interactions.txt",
    sep="\t", index=False)

# the PAV model did not perform well, so we will not consider it
tableS17_data = dt.fread(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_benchmark_RF_models_SHAP_interactions.txt").to_pandas()
ts17_no_pav = tableS17_data[tableS17_data.Model != "PAV"].copy(deep=True)


def calc_z_score(interaction_values):
    '''Calculate the z-score from a pandas series of SHAP interaction values
    The z-scores have a mean of 0 and std of 1, thus a normalized SHAP interaction
    score indicates how extreme (how unusually large or small) the interaction
    score is for a feature pair compared to all other feature pairs in the model.
    (The number of standard deviations away from the mean.)'''
    return ((interaction_values - interaction_values.mean()) / interaction_values.std())


# apply normalization to the interaction scores within each model
# this will help with comparisons across models.
for model in ts17_no_pav.Model.unique():
    model_mask = ts17_no_pav.Model == model
    ts17_no_pav.loc[model_mask, "Z_Interaction"] = calc_z_score(
        ts17_no_pav.loc[model_mask, "Interaction"])


# check for duplicates in ts17_no_pav
ts17_no_pav["Pair"] = ts17_no_pav.apply(
    lambda x: tuple(sorted([x.Feature1, x.Feature2])), axis=1)
ts17_no_pav.sort_values("Z_Interaction", ascending=False, inplace=True)
# some pairs appear 2 to 6 times (but encodings are a combo of PAV and/or CNV)
ts17_no_pav.groupby("Pair").count().value_counts()

# variant-variant interactions; will be greater than the gene-gene interactions
ts17_no_pav.groupby(["Feature1_Data", "Feature2_Data"]).count()
# excluding the PAV model, 69486 unique gene-gene interactions
len(get_unique_gp(ts17_no_pav, "Gene1", "Gene2"))

# Need to collapse the variant-variant interaction types to 6 types:
ts17_no_pav["GI_Type"] = ts17_no_pav.Feature1_Data +\
    "-" + ts17_no_pav.Feature2_Data
ts17_no_pav.GI_Type.replace('CNV-PAV', 'PAV-CNV', inplace=True)
ts17_no_pav.GI_Type.replace('PAV-SNP', 'SNP-PAV', inplace=True)
ts17_no_pav.GI_Type.replace('CNV-SNP', 'SNP-CNV', inplace=True)
ts17_no_pav.GI_Type.nunique()  # 6 :)

# How many unique variant-variant interactions are there in total?
# some PAV-CNV or CNV-PAV interactions are the "same", except which feature has
# a PAV or CNV encoding is different in each case.
ts17_no_pav[['GI_Type', 'Pair']].drop_duplicates(keep='first').shape  # 210758

# How many unique gene-gene interactions are there per variant-variant type?
snp_snp = ts17_no_pav[ts17_no_pav.GI_Type == "SNP-SNP"]
len(get_unique_gp(snp_snp, "Gene1", "Gene2"))  # 65859 SNP-SNP
pav_pav = ts17_no_pav[ts17_no_pav.GI_Type == "PAV-PAV"]
len(get_unique_gp(pav_pav, "Gene1", "Gene2"))  # 24376 PAV-PAV
cnv_cnv = ts17_no_pav[ts17_no_pav.GI_Type == "CNV-CNV"]
len(get_unique_gp(cnv_cnv, "Gene1", "Gene2"))  # 18148 CNV-CNV
snp_pav = ts17_no_pav[ts17_no_pav.GI_Type == "SNP-PAV"]
len(get_unique_gp(snp_pav, "Gene1", "Gene2"))  # 48704 SNP-PAV
snp_cnv = ts17_no_pav[ts17_no_pav.GI_Type == "SNP-CNV"]
len(get_unique_gp(snp_cnv, "Gene1", "Gene2"))  # 27204 SNP-CNV
pav_cnv = ts17_no_pav[ts17_no_pav.GI_Type == "PAV-CNV"]
len(get_unique_gp(pav_cnv, "Gene1", "Gene2"))  # 9698 PAV-CNV

# How many unique variant-variant interactions by type?
snp_snp.Pair.nunique()  # 65859
pav_pav.Pair.nunique()  # 24376
cnv_cnv.Pair.nunique()  # 18148
snp_pav.Pair.nunique()  # 62387
snp_cnv.Pair.nunique()  # 30290
pav_cnv.Pair.nunique()  # 9698    total: 210,758

### Figure S12. UpSet diagram of variant-variant shap interactions------------------#
sets = {"SNP-SNP": get_unique_gp(snp_snp, "Gene1", "Gene2"),
        "PAV-PAV": get_unique_gp(pav_pav, "Gene1", "Gene2"),
        "CNV-CNV": get_unique_gp(cnv_cnv, "Gene1", "Gene2"),
        "SNP-PAV": get_unique_gp(snp_pav, "Gene1", "Gene2"),
        "SNP-CNV": get_unique_gp(snp_cnv, "Gene1", "Gene2"),
        "PAV-CNV": get_unique_gp(pav_cnv, "Gene1", "Gene2")}


def pair2sets(sets_dict):
    # Invert the dictionary: map each gene pair to all sets it belongs to.
    # A set refers to all the variant-variant type categories.
    pair_to_sets = defaultdict(set)
    for set_name, pairs in sets_dict.items():  # var-var: set of gene-gene pairs
        for pair in pairs:
            pair_to_sets[pair].add(set_name)
    # create the var-var memberships list
    memberships = list(pair_to_sets.values())
    # create a count for each membership (all 1s if unique)
    data = from_memberships(memberships, data=[1]*len(memberships))
    return data


UpSet(pair2sets(sets), subset_size="count", show_percentages="{:.1%}").plot()
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_6/shap_interaction/Figure_S12_num_unique_GI_in_shap_interactions_upset_plot.pdf")
plt.close()
del sets
### ---------------------------------------------------------------------------#

'''Table S. Which variant-variant types tend to have higher interaction scores?
Caveat: these are comparisons of the absolute value of the summed SHAP
        interaction values across all isolates. Within each variant-variant type,
        these interaction values can come from different models. A normalization
        I may be able to apply is to take the mean(abs(shap interaction value))
        of a feature pair across all isolates in the model.
        '''
variant_pairs = {"snp_snp": snp_snp, "pav_pav": pav_pav, "cnv_cnv": cnv_cnv,
                 "snp_pav": snp_pav, "snp_cnv": snp_cnv, "pav_cnv": pav_cnv}
res = {"Comparison": [], "Statistic": [], "p-value": []}
for name1, data1 in variant_pairs.items():
    for name2, data2 in variant_pairs.items():
        if name1 != name2:
            stat, p = mannwhitneyu(data1["Interaction"].abs(),
                                   data2["Interaction"].abs(), alternative="greater")
            print(f"{name1} vs {name2}: U={stat:.2f}, p={p:.4e}")
            res["Comparison"].append(f"{name1} > {name2}")
            res["Statistic"].append(f"{stat:.2f}")
            res["p-value"].append(f"{p:.4e}")
            del stat, p

pd.DataFrame(res).to_csv(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S_benomyl_benchmark_RF_models_SHAP_interactions_mwu_results.txt",
    sep="\t", index=False)

# smallest float values is: (useful to know when p-value is 0)
np.finfo(float).tiny  # 2.2250738585072014e-308
### ---------------------------------------------------------------------------#

################################################################################
# 4. Enrichment of known genetic interactions in the SHAP interactions
################################################################################
# BioGRID database genetic interactions
biogrid = dt.fread("Data/BioGRID/yeast_gi_biogrid.txt").to_pandas()
biogrid = biogrid.iloc[:, [5, 6, 7, 8, 11, 13, 14, 17, 36]]
biogrid.columns = ["Systematic Name Interactor A", "Systematic Name Interactor B",
                   "Standard Name Interactor A", "Standard name Interactor B",
                   "Evidence", "Author", "PMID", "Throughput", "Organism"]
biogrid = biogrid.loc[biogrid.Organism ==
                      "Saccharomyces cerevisiae (S288c)", :]
biogrid = biogrid.loc[biogrid.Evidence.str.strip().isin([
    "Synthetic Growth Defect", "Synthetic Lethality", "Synthetic Rescue",
    "Negative Genetic", "Positive Genetic"]), :]  # remove overexpression gene pairs

biogrid_gp = get_unique_gp(
    biogrid, "Systematic Name Interactor A", "Systematic Name Interactor B")
len(biogrid_gp)  # 438546

# Costanzo et al. 2021 benomyl genetic interactions
costanzo = pd.read_excel(
    "Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",
    engine="openpyxl", sheet_name="Genome-scale_Benomyl")

# filter the genetic interaction network using the strict criteria, which
# represents benomyl interactions with high confidence
costanzo_ben_strict = costanzo.loc[(costanzo.mean_condition_epsilon.abs() > 0.12) &
                                   (costanzo.condition_p_value < 0.05), :]
costanzo_ben_strict.insert(0, "array_gene", costanzo_ben_strict.apply(
    lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ben_strict.insert(0, "query_gene", costanzo_ben_strict.apply(
    lambda x: x.query_orf.split("_")[0], axis=1))

costanzo_ctrl_strict = costanzo.loc[(costanzo.mean_reference_epsilon.abs() > 0.12) &
                                    (costanzo.reference_p_value < 0.05), :]  # control condition interactions with high confidence
costanzo_ctrl_strict.insert(0, "array_gene", costanzo_ctrl_strict.apply(
    lambda x: x.array_orf.split("_")[0], axis=1))
costanzo_ctrl_strict.insert(0, "query_gene", costanzo_ctrl_strict.apply(
    lambda x: x.query_orf.split("_")[0], axis=1))

# get the unique gene pairs
costanzo_ben_strict_gp = get_unique_gp(
    costanzo_ben_strict, "query_gene", "array_gene")
costanzo_ctrl_strict_gp = get_unique_gp(
    costanzo_ctrl_strict, "query_gene", "array_gene")
len(costanzo_ben_strict_gp)  # 3472 unique interactions
len(costanzo_ctrl_strict_gp)  # 3417
len(biogrid_gp.union(costanzo_ctrl_strict_gp).union(
    costanzo_ben_strict_gp))  # 441,520 unique GIs

### Figure S11 Experimentally verified GI -------------------------------------#
sets = {"Costanzo Control": costanzo_ctrl_strict_gp,
        "Costanzo Benomyl": costanzo_ben_strict_gp,
        "BioGRID": biogrid_gp}
venn(sets)  # total: 441520 unique GIs
plt.title("Experimentally verified genetic interactions")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_6/shap_interaction/Figure_S11_experimentally_verified_GIs_strict_venn_diagram.pdf")
plt.close()
del sets
### ---------------------------------------------------------------------------#

# How many Costanzo & BioGRID interactions were identified by SHAP?
len(get_unique_gp(ts17_no_pav, "Gene1", "Gene2").intersection(biogrid_gp))  # 6354
len(get_unique_gp(ts17_no_pav, "Gene1", "Gene2").intersection(
    costanzo_ben_strict_gp))  # 1
len(get_unique_gp(ts17_no_pav, "Gene1", "Gene2").intersection(
    costanzo_ctrl_strict_gp))  # 59
len(get_unique_gp(ts17_no_pav, "Gene1", "Gene2").intersection(
    biogrid_gp)).intersection(costanzo_ctrl_strict_gp).intersection(
    costanzo_ben_strict_gp)  # 6358


verified_GIs_strict = {"ben": {}, "ctrl": {}, "biogrid": {}}
for name, data in variant_pairs.items():  # from line 605
    # strict threshold
    data_gp = get_unique_gp(data, "Gene1", "Gene2")
    print(f"Number of {name} interactions in benomyl GI network:\
        {len(data_gp.intersection(costanzo_ben_strict_gp))} and the control GI network:\
        {len(data_gp.intersection(costanzo_ctrl_strict_gp))} using the strict threshold.")
    verified_GIs_strict["ben"][name] = data_gp.intersection(
        costanzo_ben_strict_gp)
    verified_GIs_strict["ctrl"][name] = data_gp.intersection(
        costanzo_ctrl_strict_gp)
    print(f"Number of {name} interactions in BioGRID:\
        {len(data_gp.intersection(biogrid_gp))}")
    verified_GIs_strict["biogrid"][name] = data_gp.intersection(biogrid_gp)
    #
    del data_gp

'''
Number of snp_snp interactions in benomyl GI network: 1
               and the control GI network: 55
Number of snp_snp interactions in BioGRID: 6092
Number of pav_pav interactions in benomyl GI network: 1
               and the control GI network: 22
Number of pav_pav interactions in BioGRID: 2419
Number of cnv_cnv interactions in benomyl GI network: 0
               and the control GI network: 8
Number of cnv_cnv interactions in BioGRID: 1897
Number of snp_pav interactions in benomyl GI network: 1
               and the control GI network: 47
Number of snp_pav interactions in BioGRID: 4611
Number of snp_cnv interactions in benomyl GI network: 0
               and the control GI network: 14
Number of snp_cnv interactions in BioGRID: 2327
Number of pav_cnv interactions in benomyl GI network: 1
               and the control GI network: 3
Number of pav_cnv interactions in BioGRID: 798
'''

# Add columns to Table S17
ts17_no_pav.insert(11, "In_BioGRID", ts17_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in biogrid_gp else 0, axis=1))
ts17_no_pav.insert(12, "In_Costanzo_Benomyl_Strict", ts17_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ben_strict_gp else 0, axis=1))
ts17_no_pav.insert(13, "In_Costanzo_Ctrl_Strict", ts17_no_pav.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in costanzo_ctrl_strict_gp else 0, axis=1))
'''sanity check:
tmp = ts17_no_pav.apply(lambda x: 1 if frozenset(
    (x.Gene1, x.Gene2)) in biogrid_gp else 0, axis=1)
tmp = pd.concat([ts17_no_pav[['Gene1', 'Gene2', 'GI_Type']],
                tmp], axis=1, ignore_index=True)
len(get_unique_gp(tmp[tmp[3]==1], 0, 1)) # 6354 (see line 691)
len(get_unique_gp(tmp.loc[(tmp[3]==1) & (tmp[2]=="SNP-SNP")], 0, 1)) # 6092
len(get_unique_gp(tmp.loc[(tmp[3]==1) & (tmp[2]=="PAV-PAV")], 0, 1)) # 2419
'''  # looking good!

ts17_no_pav.to_excel(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_benchmark_RF_models_SHAP_interactions_expanded_pav_model_removed.xlsx",
    index=False)

# How many unique gene-gene interactions were identified different genetic variant pair types?
ts17_no_pav = pd.read_excel(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_benchmark_RF_models_SHAP_interactions_expanded_pav_model_removed.xlsx")

unique_GIs = get_unique_gp(ts17_no_pav, "Gene1", "Gene2")
len(unique_GIs)  # 69486 unique gene-gene interactions (PAV model excl.)
gi_type_pairs = ts17_no_pav.groupby("GI_Type").apply(
    lambda x: get_unique_gp(x, "Gene1", "Gene2")).to_dict()
unique_GIs_variant_pair_types = {}
for gi_pair in unique_GIs:
    unique_GIs_variant_pair_types[gi_pair] = []  # list of variant pair types
    for var_pair_type, gi_pairs in gi_type_pairs.items():
        if gi_pair in gi_pairs:
            if var_pair_type not in unique_GIs_variant_pair_types[gi_pair]:
                unique_GIs_variant_pair_types[gi_pair].append(var_pair_type)


len(unique_GIs_variant_pair_types)  # 69486
unique_GIs_variant_pair_types = pd.DataFrame.from_dict(
    unique_GIs_variant_pair_types, orient="index")
unique_GIs_variant_pair_types.insert(
    6, "counts", unique_GIs_variant_pair_types.apply(
        lambda x: x.count(), axis=1))
unique_GIs_variant_pair_types.loc[
    unique_GIs_variant_pair_types.counts > 1].shape
# ^59610 gene-gene interactions were identified by 2 or more variant pair types
unique_GIs_variant_pair_types.index.name = "gene_pair"
unique_GIs_variant_pair_types.reset_index(inplace=True)
unique_GIs_variant_pair_types = pd.concat([
    pd.DataFrame(unique_GIs_variant_pair_types.gene_pair.to_list()),
    unique_GIs_variant_pair_types.iloc[:, 1:]], axis=1, ignore_index=True)
unique_GIs_variant_pair_types.columns = [
    "Gene1", "Gene2", "var_pair_1", "var_pair_2", "var_pair_3",
    "var_pair_4", "var_pair_5", "var_pair_6", "NumVarPairs"]

# How many of the experimentally validated GIs were identified by multiple
# variant-variant interaction types?
ts17_no_pav["Validated"] = ts17_no_pav[[
    "In_BioGRID", "In_Costanzo_Benomyl_Strict", "In_Costanzo_Ctrl_Strict"]].any(axis=1).astype(int)

val_GIs = get_unique_gp(ts17_no_pav[ts17_no_pav.Validated == 1])  # 6358
val_GI_type_pairs = ts17_no_pav[ts17_no_pav.Validated == 1].groupby("GI_Type").apply(
    lambda x: get_unique_gp(x, "Gene1", "Gene2")).to_dict()

{key: len(value) for key, value in val_GI_type_pairs.items()}
# {'CNV-CNV': 1898, 'PAV-CNV': 798, 'PAV-PAV': 2420, 'SNP-CNV': 2327, 'SNP-PAV': 4615, 'SNP-SNP': 6095}

val_GIs_variant_pair_types = {}
for gi_pair in val_GIs:
    val_GIs_variant_pair_types[gi_pair] = []  # list of variant pair types
    for var_pair_type, gi_pairs in val_GI_type_pairs.items():
        if gi_pair in gi_pairs:
            if var_pair_type not in val_GIs_variant_pair_types[gi_pair]:
                val_GIs_variant_pair_types[gi_pair].append(var_pair_type)


len(val_GIs_variant_pair_types)  # 6358
val_GIs_variant_pair_types = pd.DataFrame.from_dict(
    val_GIs_variant_pair_types, orient="index")
val_GIs_variant_pair_types.insert(
    6, "counts", val_GIs_variant_pair_types.apply(
        lambda x: x.count(), axis=1))
val_GIs_variant_pair_types.loc[
    val_GIs_variant_pair_types.counts > 1].shape  # 5645 genes
val_GIs_variant_pair_types.index.name = "gene_pair"
val_GIs_variant_pair_types.reset_index(inplace=True)
val_GIs_variant_pair_types = pd.concat([
    pd.DataFrame(val_GIs_variant_pair_types.gene_pair.to_list()),
    val_GIs_variant_pair_types.iloc[:, 1:]], axis=1, ignore_index=True)
val_GIs_variant_pair_types.columns = [
    "Gene1", "Gene2", "var_pair_1", "var_pair_2", "var_pair_3",
    "var_pair_4", "var_pair_5", "var_pair_6", "NumVarPairs"]

# Add the column to unique_GIs_variant_pair_types to indicate if the GI was validated
unique_GIs_variant_pair_types["Validated"] = unique_GIs_variant_pair_types.apply(
    lambda x: 1 if frozenset((x.Gene1, x.Gene2)) in val_GIs else 0, axis=1)
unique_GIs_variant_pair_types.Validated.sum()  # 6358
unique_GIs_variant_pair_types.to_csv(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S18_validated_GIs_identified_by_multiple_variant_pair_types.txt",
    index=False, sep="\t", header=True)

# Of the gene pairs that are not experimentally validated, are more GIs identified
# by multiple variant-variant interaction types compared to the validated GIs?
tmp = unique_GIs_variant_pair_types.loc[
    ~unique_GIs_variant_pair_types.set_index(
        ["Gene1", "Gene2"]).index.isin(
            val_GIs_variant_pair_types.set_index(["Gene1", "Gene2"]).index)]
# 63128 gene-gene interactions were not experimentally validated

tmp.NumVarPairs.value_counts()
# NumVarPairs
# 3    18781
# 2    18475
# 4    11119
# 1     9163
# 5     4636
# 6      954
# Name: count, dtype: int64
val_GIs_variant_pair_types.NumVarPairs.value_counts()
# NumVarPairs
# 3    2112
# 2    1821
# 4    1177
# 1     713
# 5     456
# 6      79
# Name: count, dtype: int64

# Determine enrichment of the combined validated GIs among SHAP interactions
ts17_no_pav["Validated"] = ts17_no_pav[[
    "In_BioGRID", "In_Costanzo_Benomyl_Strict", "In_Costanzo_Ctrl_Strict"]].any(axis=1).astype(int)
all_validated = biogrid_gp.union(
    costanzo_ctrl_strict_gp).union(costanzo_ben_strict_gp)

# get all possible feature pairs
# ignore if the PAV/CNV feature pair encodings are the same? Then, i'm looking at
# snp-snp, snp-orf, and orf-orf pairs.
dpath = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_rf"
benomyl_features = pd.read_csv(
    f"{dpath}/Features_one_variant_per_gene_benomyl_500ugml_snp_pav_cnv.txt", sep=",", nrows=1)
benomyl_features = pd.DataFrame(benomyl_features.columns[1:])
benomyl_features[0] = benomyl_features[0].str.replace("^X", "", regex=True)
benomyl_features[0] = benomyl_features[0].str.replace("\.", "-", regex=True)
benomyl_features[0] = benomyl_features[0].str.replace(";PAV$", "", regex=True)
benomyl_features[0] = benomyl_features[0].str.replace(";CNV$", "", regex=True)
benomyl_features = benomyl_features.drop_duplicates(keep="first")
all_pairs = pd.DataFrame(list(combinations(
    benomyl_features[0].values, 2)), columns=["Feature1", "Feature2"])
all_pairs.insert(2, "Gene1", all_pairs.apply(
    lambda x: feature2gene(x["Feature1"]), axis=1))
all_pairs.insert(3, "Gene2", all_pairs.apply(
    lambda x: feature2gene(x["Feature2"]), axis=1))
all_pairs.shape  # 258840 unique possible snp-snp, snp-orf, orf-orf pairs
# all_pairs.isna().sum() # all good!

# Drop duplicate ts17_no_pav feature pairs (ignore the variant-variant type).
# The reason is so that the contingency table counts sum up as they should.
# keep the first one (has the highest interaction score -- see caveat on line 637)
ts17_no_pav.shape  # (310356, 16)
ts17_no_pav = ts17_no_pav.sort_values('Interaction').drop_duplicates(
    subset=["Pair"], keep="first")
ts17_no_pav.shape  # (181758, 16)
ts17_no_pav.groupby("Pair").count().value_counts()  # now they only appear once
# before, they were appearing 2 to 6 times because variant-variant types differed

### See: ----------------------------------------------------------------------#
# Which verified GIs were identified by multiple variant-variant interactions?
overlap_counts = defaultdict(lambda: defaultdict(set))

for sample, categories in verified_GIs_strict.items():
    for category, pair_set in categories.items():
        for pair in pair_set:
            overlap_counts[pair][sample].add(category)

rows = []
for pair, sample_dict in overlap_counts.items():
    gene1, gene2 = sorted(pair)
    for sample, categories in sample_dict.items():
        if len(categories) > 1:
            rows.append({
                "gene1": gene1, "gene2": gene2,
                "sample": sample, "n_categories": len(categories),
                "categories": "; ".join(sorted(categories))})
            # print(f"{pair} appears in {len(categories)} categories for {sample}: {sorted(categories)}")

df = pd.DataFrame(rows)
df = df.sort_values(by=['sample', 'n_categories'], ascending=[True, False])
df.groupby("sample").head(20)
'''these are not necessarily the top [or the best] because many gene pairs were
identified by 4 variant pair types, for example.'''
df.shape  # (5694, 5)

# was found by all three datasets (biogrid, ben, and ctrl)
df.loc[(df.gene1 == 'YGR078C') & (df.gene2 == "YNL153C")]
### ----------------------------------------------------------------------------#


def enrichment_direction(k, n, C, G):
    """determine direction of enrichment
    if >= 1: + (overrepresented)
    if < 1: - (underrepresented)
    k: number of pairs above shap rank threshold and are experimentally validated GIs
    n: total number of genes in shap
    C: total number of genes (above shap rank + background) and are known GIs
    G: total number of genes (above shap rank + background)"""
    return ((k/C)/(n/G))


# Rank the SHAP interaction scores and calculate enrichment by percentiles
ts17_no_pav["Rank"] = np.nan
for model in ts17_no_pav["Model"].unique().tolist():
    ts17_no_pav.loc[ts17_no_pav.Model == model, "Rank"] = ts17_no_pav.loc[
        ts17_no_pav.Model == model, "Interaction"].rank(method="average", ascending=False)

# Add the GI_Type values into the Pair values tuples
ts17_no_pav["Pair_with_info"] = ts17_no_pav.apply(lambda x: tuple(sorted(
    [x.GI_Type, x.Feature1, x.Feature2])), axis=1)

# Calculate enrichment for each variant-variant type and percentile
results = []

for percentile in [.01, .05, .1, .15, .2, .25]:
    for model in ts17_no_pav["Model"].unique().tolist():
        # quantile returns the values between 0 <= percentile, so <= .01 is the top 1%
        percentile_rank_val = ts17_no_pav.loc[
            ts17_no_pav.Model == model, "Rank"].quantile(percentile)
        #
        # target subset of SHAP interactions
        top_shap_int = ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Rank <= percentile_rank_val), :]
        # background subset of SHAP interactions
        bg_shap_int = ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Rank > percentile_rank_val), :]
        assert (top_shap_int.shape[0] + bg_shap_int.shape[0]) == ts17_no_pav.loc[
            ts17_no_pav.Model == model, :].shape[0]
        #
        # shap interactions are validated and are above the rank threshold
        a = top_shap_int.loc[top_shap_int.Validated == 1, "Validated"].sum()
        # shap interactions are not validated and are above the rank threshold
        b = top_shap_int.loc[
            top_shap_int.Validated == 0, "Pair"].nunique()
        # shap interactions are validated and are not above the rank threshold
        c = bg_shap_int.loc[bg_shap_int.Validated == 1, "Validated"].sum()
        # shap interactions are not validated and are not above the rank threshold
        d = bg_shap_int.loc[
            bg_shap_int.Validated == 0, "Pair"].nunique()
        #
        assert (a+c) == ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Validated == 1), "Validated"].sum()
        assert (a+b) == ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Rank <= percentile_rank_val), "Pair"].nunique()
        # These other assertions are wrong, because some duplicate
        # (feature, feature, var-type) end up in both top and bg sets, so the
        # numbers on the left and right hand side are very similar but not equal.
        assert (c+d) == ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Rank > percentile_rank_val), "Pair"].nunique()
        assert (b+d) == ts17_no_pav.loc[(ts17_no_pav.Model == model) & (
            ts17_no_pav.Validated == 0), "Pair"].nunique()
        assert (a+b+c+d) == ts17_no_pav.loc[
            ts17_no_pav.Model == model, "Pair"].nunique()
        #
        odds, pval = fisher_exact([[a, b], [c, d]])
        if enrichment_direction(k=a, n=a+b, C=a+c, G=a+b+c+d) >= 1:
            direction = "+"
        else:
            direction = "-"
        # store the results
        # results[var_pair] = [a, b, c, d, odds, pval, direction]
        results.append([model, percentile, percentile_rank_val,
                        a, b, c, d, odds, pval, direction])
        print(
            f"odds ratio={odds:.2f}, p-value={pval:.4e}, enrichment direction={direction}")
        #
        # direction2 = (a/(a+c)) / (b/(b+d))  # chat gpt recommended formula, the results are very similar to the above
        # print(f"enrichment direction (method 2)={direction2:.2f}")
        del a, b, c, d, odds, pval, direction
        del top_shap_int, bg_shap_int, percentile_rank_val


out = pd.DataFrame(results, columns=["model", "percentile", "percentile_rank_val",
                                     "a", "b", "c", "d", "Odds_Ratio", "p-value",
                                     "Enrichment_Direction"])
q_values = false_discovery_control(
    out["p-value"], method="bh")  # Adjust the p-values
out["q-value"] = q_values
out.sort_values(by="q-value").to_excel(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S19_SHAP_interaction_validated_GI_enrichment.xlsx")

### ---------------------------------------------------------------------------#

################################################################################
# 5. SHAP interaction plots (Figures 6B-F)
################################################################################
# Which gene pairs have the highest interaction scores?
# Load the shap interaction data
ts17_no_pav = pd.read_excel(
    "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_benchmark_RF_models_SHAP_interactions_expanded_pav_model_removed.xlsx")
ts17_no_pav = ts17_no_pav.sort_values('Interaction', ascending=False)

# Which gene pairs have the highest number of feature interactions?
ts17_no_pav.groupby(["Gene1", "Gene2"]).count().sort_values(
    "Model", ascending=False)

# highest shap interaction scores by model type
top5_snp_model = ts17_no_pav.loc[ts17_no_pav.Model == "SNP"].iloc[:5, :]
top5_snp_model["rank"] = np.arange(1, 6)
top5_cnv_model = ts17_no_pav.loc[ts17_no_pav.Model == "CNV"].iloc[:5, :]
top5_cnv_model["rank"] = np.arange(1, 6)
top5_snp_pav_model = ts17_no_pav.loc[
    ts17_no_pav.Model == "SNP + PAV"].iloc[:5, :]
top5_snp_pav_model["rank"] = np.arange(1, 6)
top5_snp_cnv_model = ts17_no_pav.loc[
    ts17_no_pav.Model == "SNP + CNV"].iloc[:5, :]
top5_snp_cnv_model["rank"] = np.arange(1, 6)
top5_pav_cnv_model = ts17_no_pav.loc[
    ts17_no_pav.Model == "PAV + CNV"].iloc[:5, :]
top5_pav_cnv_model["rank"] = np.arange(1, 6)
top5_snp_pav_cnv_model = ts17_no_pav.loc[
    ts17_no_pav.Model == "SNP + PAV + CNV"].iloc[:5, :]
top5_snp_pav_cnv_model["rank"] = np.arange(1, 6)

# Also plot the interaction for the benomyl network gene pair ('YGR078C', 'YNL153C')
# first, rank the features for each model
ts17_no_pav["rank"] = np.nan
for model in ts17_no_pav["Model"].unique().tolist():
    ts17_no_pav.loc[ts17_no_pav.Model == model, "rank"] = ts17_no_pav.loc[
        ts17_no_pav.Model == model, "Interaction"].rank(method="average", ascending=False)


verified_ben30_gi = pd.concat([
    ts17_no_pav.loc[(ts17_no_pav.Gene1 == 'YGR078C') & (
        ts17_no_pav.Gene2 == "YNL153C")],
    ts17_no_pav.loc[(ts17_no_pav.Gene2 == 'YGR078C') & (
        ts17_no_pav.Gene1 == "YNL153C")]], axis=0)


### Draw shap interaction plots -----------------------------------------------#
# Load the feature tables
snp = dt.fread("Data/Peter_2018/geno_corrected.csv").to_pandas()
snp = snp.set_index("ID")
pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
pav.columns = fix_orf_names(pd.Series(pav.columns))
cnv = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
cnv.columns = fix_orf_names(pd.Series(cnv.columns))

# Load the SHAP values
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_rf"
shap_snp = pd.read_csv(
    f"{d}/snp/SHAP_values_sorted_snp_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)
shap_cnv = pd.read_csv(
    f"{d}/cnv/SHAP_values_sorted_cnv_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)
shap_cnv.columns = fix_orf_names(pd.Series(shap_cnv.columns))
shap_pav_cnv = pd.read_csv(
    f"{d}/pav_cnv/SHAP_values_sorted_integrated_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)
shap_snp_cnv = pd.read_csv(
    f"{d}/snp_cnv/SHAP_values_sorted_integrated_snp_cnv_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)
shap_snp_pav = pd.read_csv(
    f"{d}/snp_pav/SHAP_values_sorted_integrated_snp_pav_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)
shap_snp_pav_cnv = pd.read_csv(
    f"{d}/snp_pav_cnv/SHAP_values_sorted_integrated_snp_pav_cnv_one_variant_per_gene_benomyl_500ugml_training.txt", sep="\t", index_col=0)


def interaction_plot(feature1, feature2, shap1, shap2, row, ax, cbar_label="geno", reverse=False):
    """Plot the interaction between two features and their SHAP values. Median
    SHAP values are calculated across feature1 and feature2 categories.
    X-axis: feature1 values
    Y-axis: SHAP values for feature1
    Hue: feature2 values
    """
    if cbar_label == "geno":
        hue = feature2
        if reverse:
            title = f"{row.Gene1} {row.Feature1_Data} genotypes"
        else:
            title = f"{row.Gene2} {row.Feature2_Data} genotypes"
    elif cbar_label == "SHAP":
        hue = shap2
        if reverse:
            title = f"SHAP {row.Gene1}; {row.Model} model"
        else:
            title = f"SHAP {row.Gene2}; {row.Model} model"
    #
    sns.stripplot(x=feature1, y=shap1, hue=hue, palette="RdYlBu_r", alpha=.4,
                  edgecolor="black", linewidth=0.05, legend=True, jitter=0.2,
                  dodge=True, size=3, ax=ax)
    sns.lineplot(x=feature1, y=shap1, errorbar=('pi', 50), color="black",
                 linewidth=1, marker="o", markersize=5, ax=ax,
                 label=f"Median SHAP across {row.Feature1} {row.Feature1_Data}s")
    # interactions with a pav variant only have one hue, and thus a pointplot cannot be drawn.
    if hue.nunique() > 1:
        sns.pointplot(x=feature1, y=shap1, hue=hue,
                      palette="RdYlBu_r", order=np.sort(feature1.unique()),
                      estimator=np.median, dodge=.4, linestyle="none",
                      errorbar=("pi", 50), marker="o", markersize=5,
                      markeredgewidth=.5, capsize=.1, ax=ax)
        ax.legend(title=title, loc="upper right",
                  bbox_to_anchor=(1.25, 1), fontsize=5)
    if reverse:
        ax.set_xlabel(
            f"{row.Feature2} ({row.Gene2}) {row.Feature2_Data}s",
            fontsize=5)
        ax.set_ylabel(
            f"SHAP {row.Gene2}; {row.Model} model",
            fontsize=5)
    else:
        ax.set_xlabel(
            f"{row.Feature1} ({row.Gene1}) {row.Feature1_Data}s",
            fontsize=5)
        ax.set_ylabel(
            f"SHAP {row.Gene1}; {row.Model} model",
            fontsize=5)
    return ax


# Loop through the top 20 gene pairs and plot their feature values
for top_df in [verified_ben30_gi, top5_snp_model, top5_cnv_model, top5_snp_pav_model,
               top5_snp_cnv_model, top5_pav_cnv_model, top5_snp_pav_cnv_model]:
    for i, row in top_df.iterrows():
        # Get the genotype values
        if row.Feature1_Data == "CNV":
            feature1 = cnv.loc[:, row.Feature1]
            RF_feature1_name = f"{row.Feature1};CNV"
        elif row.Feature1_Data == "PAV":
            feature1 = pav.loc[:, row.Feature1]
            RF_feature1_name = f"{row.Feature1};PAV"
        elif row.Feature1_Data == "SNP":
            feature1 = snp.loc[:, row.Feature1]
            RF_feature1_name = row.Feature1
        #
        if row.Feature2_Data == "CNV":
            feature2 = cnv.loc[:, row.Feature2]
            RF_feature2_name = f"{row.Feature2};CNV"
        elif row.Feature2_Data == "PAV":
            feature2 = pav.loc[:, row.Feature2]
            RF_feature2_name = f"{row.Feature2};PAV"
        elif row.Feature2_Data == "SNP":
            feature2 = snp.loc[:, row.Feature2]
            RF_feature2_name = row.Feature2
        #
        # Get the shap values (only the 625 training isolates)
        if row.Model == "SNP":
            shap1 = shap_snp.loc[:, row.Feature1]
            shap2 = shap_snp.loc[:, row.Feature2]
        elif row.Model == "CNV":
            shap1 = shap_cnv.loc[:, row.Feature1]
            shap2 = shap_cnv.loc[:, row.Feature2]
        elif row.Model == "PAV + CNV":
            shap1 = shap_pav_cnv.loc[:, RF_feature1_name]
            shap2 = shap_pav_cnv.loc[:, RF_feature2_name]
        elif row.Model == "SNP + CNV":
            shap1 = shap_snp_cnv.loc[:, RF_feature1_name]
            shap2 = shap_snp_cnv.loc[:, RF_feature2_name]
        elif row.Model == "SNP + PAV":
            shap1 = shap_snp_pav.loc[:, RF_feature1_name]
            shap2 = shap_snp_pav.loc[:, RF_feature2_name]
        elif row.Model == "SNP + PAV + CNV":
            shap1 = shap_snp_pav_cnv.loc[:, RF_feature1_name]
            shap2 = shap_snp_pav_cnv.loc[:, RF_feature2_name]
        #
        if os.path.exists(f"Scripts/Data_Vis/Section_6/shap_interaction/Figure_6B-F_rank_{row['rank']}_model_{row.Model}_{row.GI_Type}_{row.Gene1}_{row.Gene2}_interaction.pdf"):
            continue
        #
        # Gene1 SHAP vs Gene1 feature values; Hue = Gene2 feature values
        fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(8, 8))
        ax[0, 0] = interaction_plot(
            feature1, feature2, shap1, shap2, row, ax[0, 0], cbar_label='geno')
        #
        # Gene2 SHAP vs Gene2 feature values; Hue = Gene1 feature values
        ax[0, 1] = interaction_plot(feature1=feature2, feature2=feature1,
                                    shap1=shap2, shap2=shap1, row=row, ax=ax[0, 1],
                                    cbar_label='geno', reverse=True)
        #
        # Gene1 SHAP vs Gene1 feature values; Hue = Gene2 SHAP values
        sns.stripplot(x=feature1, y=shap1, hue=shap2, palette="RdBu_r", alpha=0.5,
                      edgecolor="black", linewidth=0.1, legend=True,
                      jitter=.2, dodge=True, size=3, ax=ax[1, 0])
        ax[1, 0].set_xlabel(
            f"{row.Feature1} ({row.Gene1}) {row.Feature1_Data}s",
            fontsize=5)
        ax[1, 0].set_ylabel(
            f"SHAP {row.Gene1}; {row.Model} model", fontsize=5)
        ax[1, 0].get_legend().set_title(
            f"SHAP {row.Gene2}; {row.Model} model")
        #
        # Gene2 SHAP vs Gene2 feature values; Hue = Gene1 SHAP values
        sns.stripplot(x=feature2, y=shap2, hue=shap1, palette="RdBu_r", alpha=0.5,
                      edgecolor="black", linewidth=0.1, legend=True,
                      jitter=.2, dodge=True, size=3, ax=ax[1, 1])
        ax[1, 1].set_xlabel(
            f"{row.Feature2} ({row.Gene2}) {row.Feature2_Data}s",
            fontsize=5)
        ax[1, 1].set_ylabel(
            f"SHAP {row.Gene2}; {row.Model} model",
            fontsize=5)
        ax[1, 1].get_legend().set_title(
            f"SHAP {row.Gene1}; {row.Model} model")
        #
        for axis in ax.flat:
            axis.set_box_aspect(1)
        #
        plt.tight_layout()
        plt.savefig(
            f"Scripts/Data_Vis/Section_6/shap_interaction/Figure_6B-F_rank_{row['rank']}_model_{row.Model}_{row.GI_Type}_{row.Gene1}_{row.Gene2}_interaction.pdf")
        plt.close('all')


### ----------------------------------------------------------------------------#
