#!/usr/bin/env python3
"""This script was used for the following tasks:
- Cluster SHAP values and plot the cluster dendrogram, heatmap and clustered
  fitness label values. (Figures 5, S6, S7, S8, S9, S10; Supplementary files 15 & 16)
- Conduct a linear regression to compare fitness to median absolute SHAP values.
  (Table S14)
"""

__author__ = "Kenia Segura Abá"

import os
import pandas as pd
import datatable as dt
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from scipy.cluster import hierarchy
from scipy.stats import linregress, mannwhitneyu, pearsonr


os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/")


def subset_top20(f, target_envs, dtype=""):
    """ 
    Subset the 20 most important features based on median absolute SHAP values
    """
    env = [e for e in f.split("/")[-1].split("_") if e in target_envs][0]
    if env in target_envs:
        # SHAP values dataframe
        df = pd.read_csv(f"{f}", sep="\t", index_col=0).T
        # select the top 20 features based on median absolute shap value across isolates
        top_feat = df.abs().median(axis=1).sort_values(
            ascending=False)[:20].index
        df_top = df.loc[top_feat, :]
        if dtype == "snp":
            snp_id_map = pd.read_csv(
                "Data/Peter_2018/0_raw_data/mapping_of_geno.csv_feature_names_to_actual_feature_names_09102025.csv",
                index_col=0)
            df_top.index = df_top.index.map(snp_id_map["actual_feature"])
        return df_top


def sort_clusters(fitness, dend):  # df_top, iso_col
    """
    Sort the dendrogram clusters by median fitness value
    """
    # Get the cluster median fitness
    cluster_fitness = fitness.groupby(
        dend["leaves_color_list"]).median().sort_values()
    # Sort the clusters by median fitness
    fitness_sorted = {}
    map_df = pd.DataFrame({"ID": fitness.index.tolist(),
                          "cluster": dend["leaves_color_list"]})
    for c in cluster_fitness.index:
        fitness_sorted[c] = fitness[map_df.set_index("ID")["cluster"] == c]
    fitness_tmp = pd.DataFrame()
    sorted_colnames = []
    for c in fitness_sorted.keys():
        fitness_tmp = pd.concat([fitness_tmp, fitness_sorted[c]])
        # insert an extra 2 columns between clusters, for better visualization
        sorted_colnames.extend(fitness_sorted[c].index.tolist())
        # add two NaN columns for spacing
        sorted_colnames.extend([np.nan, np.nan])
    return (fitness_tmp, sorted_colnames)


def cluster_shap_combined(snp_files, pav_files, cnv_files, target_envs, pheno, iso_colors, dend_thrsh, cbar_lim, save_dir, by="snp"):
    '''Cluster SHAP values based on the SNP model and re-order the PAV and CNV models'''
    correlations = {}  # store the shap vs fitness linear regression results
    for env in target_envs:
        print("Clustering SHAP values for", env)
        #
        # Get file and read in SHAP values; get the top 20 features
        f_snp = [f"{f}" for f in snp_files if [
            f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_pav = [f"{f}" for f in pav_files if [
            f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        f_cnv = [f"{f}" for f in cnv_files if [
            f2 for f2 in f.split("_") if f2 in target_envs][0] == env][0]
        top_snp = subset_top20(f_snp, target_envs, dtype="snp")
        top_pav = subset_top20(f_pav, target_envs)
        top_cnv = subset_top20(f_cnv, target_envs)
        if by == "snp":
            # Cluster isolates by SNP SHAP values
            Z_iso = hierarchy.linkage(
                top_snp.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(
                Z_iso, top_snp.T, metric="euclidean")
        if by == "pav":
            # Cluster isolates by PAV SHAP values
            Z_iso = hierarchy.linkage(
                top_pav.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(
                Z_iso, top_pav.T, metric="euclidean")
        if by == "cnv":
            # Cluster isolates by CNV SHAP values
            Z_iso = hierarchy.linkage(
                top_cnv.T, method="ward", metric="euclidean")
            Z_iso_ordered = hierarchy.optimal_leaf_ordering(
                Z_iso, top_cnv.T, metric="euclidean")
        #
        optimal_leaves_list = hierarchy.leaves_list(Z_iso_ordered)
        dend = hierarchy.dendrogram(
            Z_iso_ordered, no_plot=True, color_threshold=dend_thrsh[env])
        #
        # Cluster features in the SNP, PAV, and CNV SHAP value dataframe
        Z_feat_snp = hierarchy.linkage(
            top_snp, method="ward", metric="euclidean")
        Z_feat_snp_ordered = hierarchy.optimal_leaf_ordering(
            Z_feat_snp, top_snp, metric="euclidean")
        dend_feat_snp = hierarchy.dendrogram(Z_feat_snp_ordered, no_plot=True)
        #
        Z_feat_pav = hierarchy.linkage(
            top_pav, method="ward", metric="euclidean")
        Z_feat_pav_ordered = hierarchy.optimal_leaf_ordering(
            Z_feat_pav, top_pav, metric="euclidean")
        dend_feat_pav = hierarchy.dendrogram(Z_feat_pav_ordered, no_plot=True)
        #
        Z_feat_cnv = hierarchy.linkage(
            top_cnv, method="ward", metric="euclidean")
        Z_feat_cnv_ordered = hierarchy.optimal_leaf_ordering(
            Z_feat_cnv, top_cnv, metric="euclidean")
        dend_feat_cnv = hierarchy.dendrogram(Z_feat_cnv_ordered, no_plot=True)
        #
        # Assign each dendrogram cluster a color
        cmap = mpl.colormaps.get_cmap("hsv").resampled(
            len(np.unique(dend["leaves_color_list"])))  # sample hsv color palette
        color_list = {np.unique(dend["leaves_color_list"])[i]: rgb2hex(
            cmap(i)[:3]) for i in range(cmap.N)}  # assign colors to clades
        #
        # Plot the unsorted dendrogram with colored clusters
        # set the color palette for the dendrogram
        hierarchy.set_link_color_palette(list(color_list.values()))
        fig, ax = plt.subplots(5, 1, figsize=(20, 20))
        if by == "snp":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                                            color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                                            labels=top_snp.iloc[:, optimal_leaves_list].columns)
        if by == "pav":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                                            color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                                            labels=top_pav.iloc[:, optimal_leaves_list].columns)
        if by == "cnv":
            dend_iso = hierarchy.dendrogram(Z_iso_ordered, ax=ax[0],
                                            color_threshold=dend_thrsh[env], leaf_font_size=2.5,
                                            labels=top_cnv.iloc[:, optimal_leaves_list].columns)
        print("BEFORE", color_list)
        new_color_list = {np.unique(dend_iso["leaves_color_list"])[
            i]: rgb2hex(cmap(i)[:3]) for i in range(cmap.N)}
        print("AFTER", new_color_list)
        color_list = color_list.update(new_color_list)
        iso_colors["Clade_Color"] = iso_colors["Clade_Color"].map(
            new_color_list)  # assign colors to isolate
        #
        # Re-order the SNP SHAP values dataframe and fitness values
        top_snp = top_snp.iloc[dend_feat_snp["leaves"], dend_iso["leaves"]]
        top_pav = top_pav.iloc[dend_feat_pav["leaves"], dend_iso["leaves"]]
        top_cnv = top_cnv.iloc[dend_feat_cnv["leaves"], dend_iso["leaves"]]
        fitness = pheno.loc[top_snp.columns, env]
        #
        # Sort the clusters by median fitness value and reorder all shap dataframes
        fitness_sorted, sorted_colnames = sort_clusters(fitness, dend_iso)
        # get the indices of the non-NaN columns
        no_na_cols = [col for col in sorted_colnames if not pd.isna(col)]
        # reorder the SNP dataframe
        top_snp_sorted = top_snp.loc[:, no_na_cols]
        top_pav_sorted = top_pav.loc[:, no_na_cols]
        top_cnv_sorted = top_cnv.loc[:, no_na_cols]
        na_cols = [i for i in range(len(sorted_colnames)) if pd.isna(
            sorted_colnames[i])]  # get the indices of the extra columns
        na_cols.pop(-1)
        # remove the last two NaN columns, which are not needed for the heatmap
        na_cols.pop(-1)
        #
        # Add the masked values to the remaining 2 dataframes to visualize clusters in heatmap
        for i in na_cols:
            top_snp_sorted.insert(i, f"NaN{i}", pd.Series(
                np.nan, index=top_snp_sorted.index))
            top_pav_sorted.insert(i, f"NaN{i}", pd.Series(
                np.nan, index=top_pav_sorted.index))
            top_cnv_sorted.insert(i, f"NaN{i}", pd.Series(
                np.nan, index=top_cnv_sorted.index))
            #
        top_snp_sorted.clip(
            lower=cbar_lim["snp"][env][0], upper=cbar_lim["snp"][env][1], inplace=True)
        sns.heatmap(top_snp_sorted, cmap="RdBu_r", center=0, ax=ax[1],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location": "right"})
        #
        top_pav_sorted.clip(
            lower=cbar_lim["pav"][env][0], upper=cbar_lim["pav"][env][1], inplace=True)
        sns.heatmap(top_pav_sorted, cmap="RdBu_r", center=0, ax=ax[2],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location": "right"})
        #
        top_cnv_sorted.clip(
            lower=cbar_lim["cnv"][env][0], upper=cbar_lim["cnv"][env][1], inplace=True)
        sns.heatmap(top_cnv_sorted, cmap="RdBu_r", center=0, ax=ax[3],
                    yticklabels=True, xticklabels=False,
                    cbar_kws={"orientation": "vertical", "location": "right"})
        #
        # Plot the violin plot of fitness values by dendrogram cluster
        fitness_toplot = pd.DataFrame({
            "Fitness": fitness,
            "Cluster_Color": dend_iso["leaves_color_list"]}).reset_index()
        cluster_medians = fitness_toplot.groupby("Cluster_Color").\
            median("Fitness").sort_values(
                by="Fitness", ascending=True)  # median fitness
        fitness_toplot["Cluster_Color"] = pd.Categorical(
            fitness_toplot["Cluster_Color"], categories=cluster_medians.index,
            ordered=True)  # convert cluster color into an ordered categorical variable
        fitness_toplot = fitness_toplot.sort_values(
            "Cluster_Color")  # sort by median cluster fitness
        sns.violinplot(data=fitness_toplot, x="Cluster_Color", y="Fitness",
                       hue="Cluster_Color", palette=color_list, inner="quart", ax=ax[4])
        plt.savefig(
            f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_sorted_v6.pdf")
        plt.close()
        #
        # Save the sorted SHAP values and fitness values to files
        top_snp_sorted.to_csv(
            f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_snp_sorted_v6.tsv", sep="\t")
        top_pav_sorted.to_csv(
            f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_pav_sorted_v6.tsv", sep="\t")
        top_cnv_sorted.to_csv(
            f"{save_dir}/SHAP_clustered_by_{by}_{env}_top20_cnv_sorted_v6.tsv", sep="\t")
        fitness_toplot.to_csv(
            f"{save_dir}/SHAP_clustered_by_{by}_{env}_fitness_sorted_v6.tsv", sep="\t")
        #
        # Determine the relationship between cluster fitness and absolute shap values
        if by == "snp":
            top = top_snp_sorted.copy(deep=True)
        elif by == "pav":
            top = top_pav_sorted.copy(deep=True)
        elif by == "cnv":
            top = top_cnv_sorted.copy(deep=True)
        #
        # relationship between each feature and fitness
        for i in range(top.shape[0]):
            # 625 SHAP values for feature i
            shap_row = top.iloc[i, :].loc[no_na_cols]
            print(shap_row, fitness_sorted)
            m, b, r, p, se = linregress(
                shap_row, fitness_sorted.values.flatten())
            correlations[f'sorted_by_{by}:{env}:{top.index[i]}'] = \
                {'m': m, 'b': b, 'r': r, 'p-value': p, 'se': se}
            del m, b, r, p, se
        # relationship between median shap values and fitness
        print(np.median(np.abs(top.loc[:, no_na_cols]), axis=0))
        m, b, r, p, se = linregress(np.median(np.abs(top.loc[:, no_na_cols]), axis=0),
                                    fitness_sorted.values.flatten())
        correlations[f'sorted_by_{by}:{env}:median_SHAP'] = \
            {'m': m, 'b': b, 'r': r, 'p-value': p, 'se': se}
        shap_v_fitness_df = pd.DataFrame.from_dict(
            correlations, orient='index')
        del top_snp, top_pav, top_cnv, fitness
    return shap_v_fitness_df


def save_results(mwu_res, shap_v_fitness_res, mwu_save, shap_v_fitness_save):
    # Convert Mann-Whitney U test results to a dataframe
    mwu_df = pd.DataFrame.from_dict({(i, j, k, l, m): mwu_res[i][j][k][l][m]
                                     for i in mwu_res.keys()
                                     for j in mwu_res[i].keys()
                                     for k in mwu_res[i][j].keys()
                                     for l in mwu_res[i][j][k].keys()
                                     for m in mwu_res[i][j][k][l].keys()},
                                    orient='index')
    mwu_df = mwu_df.droplevel(0)
    mwu_df.sort_index(inplace=True)
    mwu_df.reset_index(level=[3], inplace=True)
    mwu_df.insert(0, "Cluster1", [i[0] for i in mwu_df.level_3])
    mwu_df.insert(1, "Cluster2", [i[1] for i in mwu_df.level_3])
    mwu_df.drop(columns="level_3").to_csv(mwu_save)

    # Median absolute SHAP values vs median fitness of clusters linear regression
    shap_v_fitness_df = pd.DataFrame.from_dict({(i, j, k): shap_v_fitness_res[i][j][k]
                                                for i in shap_v_fitness_res.keys()
                                                for j in shap_v_fitness_res[i].keys()
                                                for k in shap_v_fitness_res[i][j].keys()},
                                               orient='index')
    shap_v_fitness_df = shap_v_fitness_df.droplevel(0)
    shap_v_fitness_df.sort_index(inplace=True)
    shap_v_fitness_df.to_csv(shap_v_fitness_save)


if __name__ == "__main__":
    # Clade & fitness information
    clades = pd.read_excel("Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls",
                           sheet_name="Table S1", skiprows=3, nrows=1011)  # isolate clades
    pheno = pd.read_csv("Data/Peter_2018/pheno.csv",
                        index_col=0)  # isolate fitness

    # Map clades to isolates and create a colormap
    clades = clades[["Standardized name", "Clades"]]  # subset relevant columns
    clades = clades.loc[clades["Standardized name"].isin(
        pheno.index)]  # diploid isolates
    clades.set_index("Standardized name", inplace=True)
    # replace NaN with Unknown
    clades.loc[clades.Clades.isna(), "Clades"] = "Unknown"

    cmap = mpl.colormaps.get_cmap("hsv").resampled(
        clades.Clades.nunique())  # sample hsv color palette
    color_list = {clades.Clades.unique()[i]: rgb2hex(
        cmap(i)[:3]) for i in range(cmap.N)}  # assign colors to clades
    iso_colors = clades['Clades'].map(color_list)  # assign colors to isolates
    iso_colors = pd.concat([clades["Clades"], iso_colors], axis=1)
    iso_colors.columns = ["Clades", "Clade_Color"]

    # map of SNPs to genes
    map_snps = pd.read_csv(
        "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
    # map of ORFs to genes
    map_orfs = pd.read_csv(
        "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")

    ############## Read in SHAP value files for optimized models ###############
    snp_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/SNP/fs"
    files = os.listdir(snp_path)
    snp_files = [
        f"{snp_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    pav_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/PAV/fs"
    files = os.listdir(pav_path)
    pav_files = [
        f"{pav_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]
    cnv_path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/CNV/fs"
    files = os.listdir(cnv_path)
    cnv_files = [
        f"{cnv_path}/{f}" for f in files if f.startswith("SHAP_values_sorted_Y")]

    # Cluster SHAP values
    target_envs = ["YPDCAFEIN40", "YPDCAFEIN50",
                   "YPDSODIUMMETAARSENITE", "YPDBENOMYL500", "YPDCUSO410MM"]
    snp_files = [f for f in snp_files if len(
        [f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    pav_files = [f for f in pav_files if len(
        [f2 for f2 in f.split("_") if f2 in target_envs]) > 0]
    cnv_files = [f for f in cnv_files if len(
        [f2 for f2 in f.split("_") if f2 in target_envs]) > 0]

    dend_thrsh_snp = {"YPDCAFEIN40": 0.014, "YPDCAFEIN50": 0.035,
                      "YPDSODIUMMETAARSENITE": 0.05, "YPDBENOMYL500": 0.012, "YPDCUSO410MM": 0.09}
    dend_thrsh_pav = {"YPDCAFEIN40": 0.08, "YPDCAFEIN50": 0.055,
                      "YPDSODIUMMETAARSENITE": 0.08, "YPDBENOMYL500": 0.09, "YPDCUSO410MM": 0.2}
    dend_thrsh_cnv = {"YPDCAFEIN40": 0.1, "YPDCAFEIN50": 0.07,
                      "YPDSODIUMMETAARSENITE": 0.25, "YPDBENOMYL500": 0.08, "YPDCUSO410MM": 0.8}
    cbar_lim = {"snp": {"YPDCAFEIN40": [-0.002, 0.002], "YPDCAFEIN50": [-0.005, 0.005], "YPDSODIUMMETAARSENITE": [-0.02, 0.02], "YPDBENOMYL500": [-0.001, 0.001], "YPDCUSO410MM": [-0.03, 0.03]},  # SHAP value limits for colorbar
                "pav": {"YPDCAFEIN40": [-0.04, 0.04], "YPDCAFEIN50": [-0.01, 0.01], "YPDSODIUMMETAARSENITE": [-0.04, 0.04], "YPDBENOMYL500": [-0.01, 0.01], "YPDCUSO410MM": [-0.04, 0.04]},
                "cnv": {"YPDCAFEIN40": [-0.04, 0.04], "YPDCAFEIN50": [-0.02, 0.02], "YPDSODIUMMETAARSENITE": [-0.06, 0.06], "YPDBENOMYL500": [-0.01, 0.01], "YPDCUSO410MM": [-0.5, 0.5]}}

    shap_v_fit_snp = cluster_shap_combined(snp_files, pav_files, cnv_files,
                                           target_envs, pheno, iso_colors, dend_thrsh_snp, cbar_lim=cbar_lim,
                                           save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="snp")
    shap_v_fit_pav = cluster_shap_combined(snp_files, pav_files, cnv_files,
                                           target_envs, pheno, iso_colors, dend_thrsh_pav, cbar_lim=cbar_lim,
                                           save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="pav")
    shap_v_fit_cnv = cluster_shap_combined(snp_files, pav_files, cnv_files,
                                           target_envs, pheno, iso_colors, dend_thrsh_cnv, cbar_lim=cbar_lim,
                                           save_dir="Scripts/Data_Vis/Section_6/shap_clusters/fs_top20", by="cnv")

    # Save the SHAP values vs fitness linear regression results
    shap_v_fitness = pd.concat(
        [shap_v_fit_snp, shap_v_fit_pav, shap_v_fit_cnv], axis=0)
    shap_v_fitness.index = shap_v_fitness.index.str.split(":", expand=True)
    shap_v_fitness.index = shap_v_fitness.index.set_names(
        ["Sorting Method", "Environment", "Feature"])
    shap_v_fitness.to_csv(
        "Scripts/Data_Vis/Section_6/shap_clusters/fs_top20/SHAP_cluster_values_vs_fitness_linreg_RF_FS_models_v5.csv")
