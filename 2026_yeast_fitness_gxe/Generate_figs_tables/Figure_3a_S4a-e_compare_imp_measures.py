#!/usr/bin/env python3
############################################################################
# This script generates:
# 1. Figure 3A (SHAP variant comparison heatmap for optimized RF models)
# 2. Figures S4A-E (A-B: Gini vs SHAP heatmap for complete RF models; C-E: Gini
#    or SHAP variant comparison heatmaps for complete and/or optimized RF models)
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Only plot for these environments
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]

# Figure S4A-B Gini vs SHAP Heatmaps
ab_data = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Table_S5_gini_vs_shap_rank_per_corr.tsv", sep="\t")
ab_data["annotations"] = ab_data.apply(
    lambda x: f"{x.rho:.2f}\n({x.NumShared})", axis=1)
complete_rho = ab_data.loc[(ab_data["Model Type"] == "complete") & (
    ab_data.Env.isin(target_envs)), :].pivot_table(
        index="Data Type", columns="Env", values="rho")
complete_ann = ab_data.loc[(ab_data["Model Type"] == "complete") & (
    ab_data.Env.isin(target_envs)), :].pivot(
        index="Data Type", columns="Env", values="annotations")
optimized_rho = ab_data.loc[(ab_data["Model Type"] == "optimized") & (
    ab_data.Env.isin(target_envs)), :].pivot_table(
        index="Data Type", columns="Env", values="rho")
optimized_ann = ab_data.loc[(ab_data["Model Type"] == "optimized") & (
    ab_data.Env.isin(target_envs)), :].pivot(
        index="Data Type", columns="Env", values="annotations")

fig, ax = plt.subplots(1, 2, figsize=(10, 5))
sns.heatmap(complete_rho, cmap="RdBu_r", annot=complete_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[0], cbar=True, square=True, center=0.5)
sns.heatmap(optimized_rho, cmap="RdBu_r", annot=optimized_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[1], cbar=True, square=True, center=0.5)
ax[0].set_title("Gini vs SHAP Rank Percentiles Baseline Models")
ax[1].set_title("Gini vs SHAP Rank Percentiles Optimized Models")
cbar = ax[0].collections[0].colorbar
cbar.set_label("rho")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_3/Figure_S4ab_gini_vs_shap_rank_per_pd_spearman.pdf")
plt.close()

# Figure 3A, S4C-E variant vs variant heatmaps
c_data = pd.read_csv(
    "Scripts/Data_Vis/Section_3/Table_S6_rank_per_corr_btwn_data_types_with_n_shared_genes.tsv", sep="\t")
c_data["annotations"] = c_data.apply(
    lambda x: f"{x.rho:.2f}\n({x.NumSharedGenes})", axis=1)
gini_opt_rho = c_data.loc[(c_data["Importance Type"].str.contains("gini")) & (
    c_data["Model Type"] == "optimized") & (c_data.Env.isin(target_envs)), :].pivot_table(
    index="Comparison", columns="Env", values="rho")
gini_opt_ann = c_data.loc[(c_data["Importance Type"].str.contains("gini")) & (
    c_data["Model Type"] == "optimized") & (c_data.Env.isin(target_envs)), :].pivot(
    index="Comparison", columns="Env", values="annotations")
shap_opt_rho = c_data.loc[(c_data["Importance Type"].str.contains("SHAP")) & (
    c_data["Model Type"] == "optimized") & (c_data.Env.isin(target_envs)), :].pivot_table(
    index="Comparison", columns="Env", values="rho")
shap_opt_ann = c_data.loc[(c_data["Importance Type"].str.contains("SHAP")) & (
    c_data["Model Type"] == "optimized") & (c_data.Env.isin(target_envs)), :].pivot(
    index="Comparison", columns="Env", values="annotations")
gini_comp_rho = c_data.loc[(c_data["Importance Type"].str.contains("gini")) & (
    c_data["Model Type"] == "complete") & (c_data.Env.isin(target_envs)), :].pivot_table(
    index="Comparison", columns="Env", values="rho")
gini_comp_ann = c_data.loc[(c_data["Importance Type"].str.contains("gini")) & (
    c_data["Model Type"] == "complete") & (c_data.Env.isin(target_envs)), :].pivot(
    index="Comparison", columns="Env", values="annotations")
shap_comp_rho = c_data.loc[(c_data["Importance Type"].str.contains("SHAP")) & (
    c_data["Model Type"] == "complete") & (c_data.Env.isin(target_envs)), :].pivot_table(
    index="Comparison", columns="Env", values="rho")
shap_comp_ann = c_data.loc[(c_data["Importance Type"].str.contains("SHAP")) & (
    c_data["Model Type"] == "complete") & (c_data.Env.isin(target_envs)), :].pivot(
    index="Comparison", columns="Env", values="annotations")

fig, ax = plt.subplots(2, 2, figsize=(10, 10))
sns.heatmap(gini_opt_rho, cmap="RdBu_r", annot=gini_opt_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[0][0], cbar=True, square=True, center=0)
sns.heatmap(shap_opt_rho, cmap="RdBu_r", annot=shap_opt_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[0][1], cbar=True, square=True, center=0)
sns.heatmap(gini_comp_rho, cmap="RdBu_r", annot=gini_comp_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[1][0], cbar=True, square=True, center=0)
sns.heatmap(shap_comp_rho, cmap="RdBu_r", annot=shap_comp_ann, fmt="",
            annot_kws={"size": 6}, ax=ax[1][1], cbar=True, square=True, center=0)
ax[0][0].set_title("Gini, Variant comparison, Optimized Models")
ax[0][1].set_title("SHAP, Variant comparison, Optimized Models")
ax[1][0].set_title("Gini, Variant comparison, Complete Models")
ax[1][1].set_title("SHAP, Variant comparison, Complete Models")
cbar = ax[0][0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[0][1].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[1][0].collections[0].colorbar
cbar.set_label("rho")
cbar = ax[1][1].collections[0].colorbar
cbar.set_label("rho")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_3/Figure_3a_S4c-e_rf_rank_per_variant_comparison_spearman.pdf")
plt.close()
