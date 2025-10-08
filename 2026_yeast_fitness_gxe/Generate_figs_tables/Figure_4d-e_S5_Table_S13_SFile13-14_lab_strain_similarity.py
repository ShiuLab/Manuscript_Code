#!/usr/bin/env python3
"""
Figure 4: SHAP values of top genes/benchmark genes VS genetic similarity to
S288C.

This script performs the following tasks:
1. Generate genetic distance matrices from SNP (0,1,2 encoding) and PAV data
   (Supplementary files 13 & 14)
2. K-means clustering of the genetic distance matrices (Figs. 4D, S5)
3. Compare the shap values between the cluster containing the lab strains and
   the most distinct cluster to it. (Fig. 4E; Table S13)
4. Determine how many of the optimized ORFs are absent in S288C or how many
   important SNP features are intergenic
"""

import os
import re
import pandas as pd
import numpy as np
import datatable as dt
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import statsmodels.api as sm
from scipy.spatial.distance import pdist, squareform, euclidean
from glob import glob
from functools import partial
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.stats import mannwhitneyu

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

################################################################################
# 1. Generate the genetic distance matrices from SNP and PAV data
################################################################################
# Read in the SNP, PAV, and test set data
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
geno = dt.fread("Data/Peter_2018/0_raw_data/geno_012.csv").to_pandas()
geno.set_index("ID", inplace=True)

# Add S288C SNP genotypes as a row of 0s (homozygous for the reference allele)
geno = geno.T
geno["S288C"] = 0
geno = geno.T
# geno.to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/geno_012_with_S288C.csv")
# geno = dt.fread("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/geno_with_S288C.csv").to_pandas()
# geno.set_index("ID", inplace=True)

# Calculate the SNP-based euclidean distance (S13 File)
# remove the test set before calculating genetic distances
geno_train = geno.loc[~geno.index.isin(test[0]), :]
eu_dist_snp = pdist(geno_train.astype(int).values, metric="euclidean")
eu_dist_snp = pd.DataFrame(squareform(eu_dist_snp),
                           columns=geno_train.index, index=geno_train.index)
# eu_dist_snp.to_csv("Scripts/Data_Vis/Section_4/genetic_distance_snp012_euclidean_to_S288C.csv")
eu_dist_snp = pd.read_csv(
    "Scripts/Data_Vis/Section_4/genetic_distance_snp012_euclidean_to_S288C.csv", index_col=0)

# Read in ORF presence/absence matrix
pav = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv", max_nrows=1,
               header=False)  # Peter et al., 2018 ORFs
pav = pav[:, 1:].to_pandas().T  # exclude 'ID' column
pav["orf"] = pav.apply(lambda x: re.sub("^X", "", x[0]), axis=1)  # fix orf IDs
pav["orf"] = pav.apply(lambda x: re.sub("\.", "-", x["orf"]), axis=1)
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")
pav = pav.merge(map_orfs[map_orfs.organism == "Saccharomyces cerevisiae S288C"],
                left_on="orf", right_on="orf", how="left")  # ORF to gene map is based on BLAST results
pav["in_ref_blast"] = pav.apply(
    lambda x: 1 if not pd.isna(x["gene"]) else 0, axis=1)
in_ref_orfs = pav.loc[pav["in_ref_blast"] == 1, "gene"].to_list()
len(in_ref_orfs)  # 5517
len(np.intersect1d(in_ref_orfs, pav.gene.to_list()))  # 5490

# Calculate the PAV-based euclidean distance (S14 File)
# read in PAV genotype data
pav_df = dt.fread("Data/Peter_2018/ORFs_pres_abs.csv").to_pandas()
pav_df.set_index("ID", inplace=True)
pav_df = pav_df.T.merge(pav.set_index(
    0)["in_ref_blast"], left_index=True, right_index=True)  # add the S288C PAV genotypes
pav_df.rename(columns={"in_ref_blast": "S288C"}, inplace=True)
pav_df.to_csv("Data/Peter_2018/ORFs_pres_abs_with_S288C.csv")

# remove the test set before calculating genetic distances
pav_train = pav_df.loc[:, ~pav_df.columns.isin(test[0])]
pav_train = pav_train.astype(int)
eu_dist_pav = pdist(pav_train.T.values, metric="euclidean")
eu_dist_pav = pd.DataFrame(squareform(eu_dist_pav),
                           columns=pav_train.columns, index=pav_train.columns)
# eu_dist_pav.to_csv("Scripts/Data_Vis/Section_4/genetic_distance_pav_euclidean_to_S288C.csv")

################################################################################
# 2. K-means clustering of genetic distance matrices
################################################################################
# Genetic distance matrices
eu_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_4/genetic_distance_snp012_euclidean_to_S288C.csv",
                          index_col=0)
eu_dist_pav = pd.read_csv("Scripts/Data_Vis/Section_4/genetic_distance_pav_euclidean_to_S288C.csv",
                          index_col=0)

# Fitness data
pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)

# to get training instances
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)

snp_train = eu_dist_snp.loc[~eu_dist_snp.index.isin(
    test[0]), ~eu_dist_snp.index.isin(test[0])]
pav_train = eu_dist_pav.loc[~eu_dist_pav.index.isin(
    test[0]), ~eu_dist_pav.index.isin(test[0])]

inertia_snp = []
inertia_pav = []
for k in range(2, 11):
    kmeans_snp = KMeans(n_clusters=k, random_state=42).fit(snp_train)
    kmeans_pav = KMeans(n_clusters=k, random_state=42).fit(pav_train)
    inertia_snp.append(kmeans_snp.inertia_)
    inertia_pav.append(kmeans_pav.inertia_)

# Plot the elbow plot
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
ax[0].plot([k for k in range(2, 11)], inertia_snp, marker="o")
ax[0].set_title("Elbow Plot for Kmeans clustering of SNP genetic distance")
ax[0].set_xlabel("Number of clusters (k)")
ax[0].set_ylabel("Inertia")
ax[1].plot([k for k in range(2, 11)], inertia_pav, marker="o")
ax[1].set_title("Elbow Plot for Kmeans clustering of PAV genetic distance")
ax[1].set_xlabel("Number of clusters (k)")
ax[1].set_ylabel("Inertia")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Kmeans_snp_or_pav_genetic_distance_elbow_plot_v3.pdf",
            bbox_inches="tight", dpi=300)
plt.close()

# Fit K-means with the optimal number of clusters
kmeans_snp = KMeans(n_clusters=6, random_state=42).fit(snp_train)
kmeans_pav = KMeans(n_clusters=4, random_state=42).fit(pav_train)

# Visualize the clusters using PCA
pca = PCA(n_components=2)
pca_snp = pca.fit(snp_train)
pca_snp_df = pca_snp.transform(snp_train)
# variance explained by each component
vexp_pca_snp = pca_snp.explained_variance_ratio_

pca_p = PCA(n_components=2)
pca_pav = pca_p.fit(pav_train)
pca_pav_df = pca_p.transform(pav_train)
vexp_pca_pav = pca_pav.explained_variance_ratio_

# SNP plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
scatter1 = axes[0].scatter(pca_snp_df[:, 0], pca_snp_df[:, 1], c=kmeans_snp.labels_,
                           label=kmeans_snp.labels_, cmap='tab10', s=60, alpha=0.7)
axes[0].set_title('K-means Clusters (SNP Features)')
axes[0].set_xlabel(f'PC1 ({vexp_pca_snp[0]:.2%} variance explained)')
axes[0].set_ylabel(f'PC2 ({vexp_pca_snp[1]:.2%} variance explained)')
axes[0].scatter(pca_snp_df[snp_train.index.get_loc("S288C"), 0],
                pca_snp_df[snp_train.index.get_loc("S288C"), 1],
                color='red', s=100, edgecolor='black', label='S288C')
legend1 = axes[0].legend(*scatter1.legend_elements(),
                         loc="upper right", title="Clusters", fontsize=6)
axes[0].add_artist(legend1)

# PAV plot
scatter2 = axes[1].scatter(pca_pav_df[:, 0], pca_pav_df[:, 1], c=kmeans_pav.labels_,
                           label=kmeans_pav.labels_, cmap='tab10', s=60, alpha=0.7)
axes[1].set_title('K-means Clusters (PAV Features)')
axes[1].set_xlabel(f'PC1 ({vexp_pca_pav[0]:.2%} variance explained)')
axes[1].set_ylabel(f'PC2 ({vexp_pca_pav[1]:.2%} variance explained)')
axes[1].scatter(pca_pav_df[pav_train.index.get_loc("in_S288C"), 0],
                pca_pav_df[pav_train.index.get_loc("in_S288C"), 1],
                color='red', s=100, edgecolor='black', label='S288C')
legend2 = axes[1].legend(*scatter2.legend_elements(),
                         loc="upper right", title="Clusters", fontsize=6)
axes[1].add_artist(legend2)
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_4b_S5_PCA_snp_or_pav_genetic_distance.pdf",
            bbox_inches="tight", dpi=300)
plt.close()


def cluster_distance(centroids):
    # Calculate distance between the cluster centers
    distances = np.zeros((len(centroids), len(centroids)))
    for i in range(len(centroids)):
        for j in range(i + 1, len(centroids)):
            distances[i, j] = euclidean(centroids[i], centroids[j])
            distances[j, i] = distances[i, j]
    return pd.DataFrame(distances)


centroids_snp = kmeans_snp.cluster_centers_  # Cluster centers
centroids_pav = kmeans_pav.cluster_centers_
dist_snp = cluster_distance(centroids_snp)
dist_pav = cluster_distance(centroids_pav)

# Identify the clusters most distinct to the cluster in which S288C and W303 are in
cluster_assignments_snp = pd.DataFrame(kmeans_snp.labels_, index=kmeans_snp.feature_names_in_,
                                       columns=["Cluster"])  # S288C is in cluster 4
cluster_assignments_pav = pd.DataFrame(kmeans_pav.labels_, index=kmeans_pav.feature_names_in_,
                                       columns=["Cluster"])  # S288C is in cluster 3

# Re-assign cluster numbers so that S288C is in cluster 0
cluster_assignments_snp["New_Cluster"] = cluster_assignments_snp["Cluster"].replace(
    {0: 2, 1: 3, 2: 4, 3: 5, 4: 0, 5: 1})
cluster_assignments_pav["New_Cluster"] = cluster_assignments_pav["Cluster"].replace(
    {0: 1, 1: 2, 2: 3, 3: 0})

# cluster 5 (originally 3) is the most distinct to the S288C cluster 0 (originally cluster 4)
[4, dist_snp.loc[4, :].idxmax()]
s288c_distinct_cluster_snp = 5
# cluster 3 (originally 2) is the most distinct to the S288C cluster 0 (originally cluster 3)
[3, dist_pav.loc[3, :].idxmax()]
s288c_distinct_cluster_pav = 3

################################################################################
# 3. Compare the fitness distributions and shap values between the cluster
#    containing the lab strains and the most distinct cluster to it. (Figs. 4E,
#    S5B; Table S13)
################################################################################
# Benchmark genes validated in S288C
ben_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_meta = dt.fread(
    "Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt", sep="\t").to_pandas()
cu_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

ben_s288c = ben_meta.loc[ben_meta["Strain Background"]
                         == "S288C", "Gene Systematic Name"].unique()
caf_s288c = caf_meta.loc[caf_meta["Strain Background"]
                         == "S288C", "Gene Systematic Name"].unique()
cu_s288c = cu_meta.loc[cu_meta["Strain Background"]
                       == "S288C", "Gene Systematic Name"].unique()
sma_s288c = sma_meta.loc[sma_meta["Strain Background"]
                         == "S288C", "Gene Systematic Name"].unique()

# Benchmark genes present in SNP and ORF data
map_snps = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")
map_snps.set_index("snp", inplace=True)
map_orfs.set_index("orf", inplace=True)

ben_snp = map_snps.loc[map_snps["Benomyl"] == 1, "gene"].unique()
ben_orf = map_orfs.loc[map_orfs["Benomyl"] == 1, "gene"].unique()
caf_snp = map_snps.loc[map_snps["Caffeine"] == 1, "gene"].unique()
caf_orf = map_orfs.loc[map_orfs["Caffeine"] == 1, "gene"].unique()
cu_snp = map_snps.loc[map_snps["CuSO4"] == 1, "gene"].unique()
cu_orf = map_orfs.loc[map_orfs["CuSO4"] == 1, "gene"].unique()
sma_snp = map_snps.loc[map_snps["Sodium_meta-arsenite"] == 1, "gene"].unique()
sma_orf = map_orfs.loc[map_orfs["Sodium_meta-arsenite"] == 1, "gene"].unique()


def fix_orf_ids(df):
    """ Fix the ORF IDs in the dataframe by removing the prefix 'X' and
    replacing '.' with '-' """
    df.columns = df.apply(lambda x: re.sub(
        "^X", "", x.name), axis=0)  # fix orf IDs
    df.columns = df.apply(lambda x: re.sub("\.", "-", x.name), axis=0)
    return df


snp_clust_map = {0: 2, 1: 3, 2: 4, 3: 5, 4: 0, 5: 1}  # {original: new}
pav_clust_map = {0: 1, 1: 2, 2: 3, 3: 0}
invert_snp_clust_map = {v: k for k, v in snp_clust_map.items()}
invert_pav_clust_map = {v: k for k, v in pav_clust_map.items()}

d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"
mwu_res = {"snp": {}, "pav": {}}
benchmark_genes = ["benomyl", "caffeine",
                   "cuso4", "sodium_meta-arsenite", "caffeine"]
for i, bench_list in enumerate([ben_s288c, caf_s288c, cu_s288c, sma_s288c, caf_s288c]):
    for mtype in ["fs", "baseline"]:  # optimized & complete RF models
        print(len(bench_list))
        if i == 0:
            env = "YPDBENOMYL500"
        elif i == 1:
            env = "YPDCAFEIN40"
        elif i == 2:
            env = "YPDCUSO410MM"
        elif i == 3:
            env = "YPDSODIUMMETAARSENITE"
        elif i == 4:
            env = "YPDCAFEIN50"
        #
        # SNPs, read in the feature name map
        snp_feat_map_shap = dt.fread(glob(
            f"{d}/SNP/{mtype}/SHAP_values_sorted_average_{env}_snp_rf_{mtype}_*_with_actual_feature_names_*")[0]).to_pandas()
        # set column names
        snp_feat_map_shap.columns = snp_feat_map_shap.iloc[0, :]
        snp_feat_map_shap = snp_feat_map_shap.iloc[1:, :]  # remove first row
        snp_feat_map_shap.set_index(snp_feat_map_shap.columns[1], inplace=True)
        snp_feat_map_shap = snp_feat_map_shap.drop(
            columns=0).to_dict()
        # read in shap values from the optimized RF models
        snp_shap = dt.fread(glob(
            f"{d}/SNP/{mtype}/SHAP_values_sorted_{env}_snp_rf_{mtype}_*")[0]).to_pandas()
        snp_shap.set_index("ID", inplace=True)
        snp_shap.columns = snp_shap.columns.map(
            snp_feat_map_shap["actual_feature"])
        #
        # drop the feature with zero shap values across all isolates
        print('snp_shap before dropping zero columns', snp_shap.shape)
        snp_shap = snp_shap.loc[:, ~snp_shap.eq(0).all(axis=0)]
        print('snp_shap after dropping zero columns', snp_shap.shape)
        #
        # subset the benchmark genes validated in s288c for an environment
        # note: all snps per gene are considered
        print('Number of S288C benchmark genes in list:', len(bench_list))
        s288c_snp_shap = snp_shap.loc[
            :, snp_shap.columns.isin(map_snps.loc[map_snps.gene.isin(bench_list)].index)].abs()  # absolute shap values
        print('Number of S288C benchmark genes in SNP shap features:', len(
            s288c_snp_shap.columns.map(map_snps.gene).unique()))
        #
        if s288c_snp_shap.shape[1] != 0:
            # plot the distributions of the benchmark gene shap values in a violin plot
            s288c_snp_shap.insert(
                0, 'Cluster', cluster_assignments_snp.loc[s288c_snp_shap.index, "New_Cluster"])
            s288c_snp_shap_box = s288c_snp_shap.melt(
                id_vars=["Cluster"], var_name="ID", value_name="SHAP")
            # x1000 and log the values for better visualization
            s288c_snp_shap_box["SHAP_log"] = np.log10(
                s288c_snp_shap_box["SHAP"] * 1000)
            sns.violinplot(data=s288c_snp_shap_box, x="Cluster",
                           y="SHAP_log", hue="Cluster", cut=0, palette='tab10')
            plt.title(f"{env} SNP Benchmark Gene SHAP Values")
            plt.ylabel("log10(SHAP Value x 1000)")
            plt.xlabel("Cluster")
            plt.savefig(
                f"Scripts/Data_Vis/Section_4/Figure_4e_{env}_snp_rf_{mtype}_bench_gene_shap_violin.pdf")
            plt.close()
            #
            # subset the cluster 0 shap values (S288C is in cluster '0')
            clust0_snp_shap = s288c_snp_shap.loc[
                s288c_snp_shap.Cluster == 0].drop(columns=["Cluster"])
            print('clust0_snp_shap', clust0_snp_shap.shape)
            # calculate the median absolute shap value per gene per isolate
            clust0_snp_shap_med = clust0_snp_shap.T.groupby(
                map_snps.loc[clust0_snp_shap.columns].gene).median()
            print('clust0_snp_shap_med', clust0_snp_shap_med.shape)
            for clust in [1, 2, 3, 4, 5]:
                print(f"Comparing cluster 0 to cluster {clust} for {env} SNPs")
                # subset the clust median absolute shap values
                clust_snp_shap = s288c_snp_shap.loc[
                    s288c_snp_shap.Cluster == clust].drop(columns=["Cluster"])
                print(f'clust{clust}_snp_shap', clust_snp_shap.shape)
                # calculate the median absolute shap value per gene per isolate
                clust_snp_shap_med = clust_snp_shap.T.groupby(
                    map_snps.loc[clust0_snp_shap.columns].gene).median()
                print(f'clust{clust}_snp_shap_med', clust_snp_shap_med.shape)
                # Are clust0 shap values significantly greater than the other cluster's shap?
                try:
                    x = clust0_snp_shap_med.values.flatten()
                    y = clust_snp_shap_med.values.flatten()
                    u, pval = mannwhitneyu(x, y, alternative="greater")
                    mwu_res["snp"][f"{env}.bench_{mtype}_med_shap_per_gene.s288c_clust0_vs_clust{clust}"] = \
                        [u, pval, clust0_snp_shap_med.shape[0], dist_snp.loc[4, invert_snp_clust_map[clust]]] +\
                        clust0_snp_shap_med.median().describe().tolist() +\
                        clust_snp_shap_med.median().describe().tolist()
                    # count, mean, std, min, 25%, 50%, 75%, max
                except ValueError:
                    # there were no benchmark genes in the optimized RF model features
                    mwu_res["snp"][f"{env}.bench_{mtype}_med_shap_per_gene.s288c_clust0_vs_clust{clust}"] = \
                        [np.nan, np.nan, 0, dist_snp.loc[4, invert_snp_clust_map[clust]]] +\
                        clust0_pav_shap_med.median().describe().tolist() +\
                        clust_pav_shap_med.median().describe()
        #
        # PAVs
        # read in shap values from the optimized RF models
        pav_shap = dt.fread(glob(
            f"{d}/PAV/{mtype}/SHAP_values_sorted_{env}_pav_rf_{mtype}_*")[0]).to_pandas()
        pav_shap.set_index("ID", inplace=True)
        pav_shap = fix_orf_ids(pav_shap)
        #
        # drop the feature with zero shap values across all isolates
        print('pav_shap before dropping zero columns', pav_shap.shape)
        pav_shap = pav_shap.loc[:, ~pav_shap.eq(0).all(axis=0)]
        print('pav_shap after dropping zero columns', pav_shap.shape)
        #
        # subset the benchmark genes validated in s288c for an environment
        # note: all orfs per gene are considered
        s288c_pav_shap = pav_shap.loc[
            :, pav_shap.columns.isin(map_orfs.loc[map_orfs.gene.isin(bench_list)].index)].abs()
        print('Number of S288C benchmark genes in PAV shap features:', len(
            s288c_pav_shap.columns.map(map_orfs.gene).unique()))
        #
        if s288c_pav_shap.shape[1] != 0:
            # plot the distributions of the benchmark gene shap values in a boxplot
            s288c_pav_shap.insert(
                0, 'Cluster', cluster_assignments_pav.loc[s288c_pav_shap.index, "New_Cluster"])
            s288c_pav_shap_box = s288c_pav_shap.melt(
                id_vars=["Cluster"], var_name="ID", value_name="SHAP")
            s288c_pav_shap_box["SHAP_log"] = np.log10(
                s288c_pav_shap_box["SHAP"] * 1000)
            sns.violinplot(data=s288c_pav_shap_box, x="Cluster",
                           y="SHAP_log", hue="Cluster", cut=0, palette='tab10')
            plt.title(f"{env} PAV Benchmark Gene SHAP Values")
            plt.ylabel("log10(SHAP Value x 1000)")
            plt.xlabel("Cluster")
            plt.savefig(
                f"Scripts/Data_Vis/Section_4/Figure_4e_{env}_pav_rf_{mtype}_bench_gene_shap_violin.pdf")
            plt.close()
            # subset the cluster 0 shap values (S288C is in cluster '0')
            clust0_pav_shap = s288c_pav_shap.loc[
                s288c_pav_shap.Cluster == 0].drop(columns=["Cluster"])
            # calculate the median shap value per gene per isolate
            clust0_pav_shap_med = clust0_pav_shap.T.groupby(
                map_orfs.loc[clust0_pav_shap.columns].gene).median()
            #
            for clust in [1, 2, 3]:
                print(
                    f"Comparing cluster 0 to cluster {clust} for {env} PAVs")
                # subset the cluster 2 shap values
                clust_pav_shap = s288c_pav_shap.loc[
                    s288c_pav_shap.Cluster == clust].drop(columns=["Cluster"])
                print(f'clust{clust}_pav_shap', clust_pav_shap.shape)
                # calculate the median shap value per gene per isolate
                clust_pav_shap_med = clust_pav_shap.T.groupby(
                    map_orfs.loc[clust_pav_shap.columns].gene).median()
                print(f'clust{clust}_pav_shap_med', clust_pav_shap_med.shape)
                #
                # Are clust0 shap values significantly greater than the other cluster's shap?
                try:
                    x = clust0_pav_shap_med.values.flatten()
                    y = clust_pav_shap_med.values.flatten()
                    u, pval = mannwhitneyu(x, y, alternative="greater")
                    mwu_res["pav"][f"{env}.bench_{mtype}_med_shap_per_gene.s288c_clust0_vs_clust{clust}"] =\
                        [u, pval, clust0_pav_shap_med.shape[0], dist_pav.loc[3, invert_pav_clust_map[clust]]] +\
                        clust0_pav_shap_med.median().describe().tolist() +\
                        clust_pav_shap_med.median().describe().tolist()
                except ValueError:
                    mwu_res["pav"][f"{env}.bench_{mtype}_med_shap_per_gene.s288c_clust0_vs_clust{clust}"] =\
                        [np.nan, np.nan, 0, dist_pav.loc[3, invert_pav_clust_map[clust]]] +\
                        clust0_pav_shap_med.median().describe().tolist() +\
                        clust_pav_shap_med.median().describe().tolist()
        #
        del snp_shap, s288c_snp_shap, clust0_snp_shap, clust0_snp_shap_med
        del clust_snp_shap, clust_snp_shap_med
        del pav_shap, s288c_pav_shap
        try:
            del clust0_pav_shap, clust0_pav_shap_med, clust_pav_shap, clust_pav_shap_med
        except NameError:
            continue
    del bench_list


mwu_res_df = pd.json_normalize(mwu_res, sep='.').T
mwu_res_df.index = pd.MultiIndex.from_tuples(
    mwu_res_df.index.str.split(".").map(tuple))
mwu_res_df.index.names = ["Variant Type", "Environment", "Description",
                          "Comparison"]
mwu_res_df = mwu_res_df.apply(lambda x: pd.Series(x[0]), axis=1)
mwu_res_df.columns = ["U", "pval", "num_benchmark_features",
                      "distance_between_clusters", "s288c_clust_shap_count",
                      "s288c_clust_shap_mean", "s288c_clust_shap_std",
                      "s288c_clust_shap_min", "s288c_clust_shap_25%",
                      "s288c_clust_shap_50%", "s288c_clust_shap_75%",
                      "s288c_clust_shap_max", "other_clust_shap_count",
                      "other_clust_shap_mean", "other_clust_shap_std",
                      "other_clust_shap_min", "other_clust_shap_25%",
                      "other_clust_shap_50%", "other_clust_shap_75%",
                      "other_clust_shap_max"]

mwu_res_df["difference_median_shap_btwn_clusters"] = mwu_res_df["s288c_clust_shap_50%"] - \
    mwu_res_df["other_clust_shap_50%"]
mwu_res_df.reset_index().to_excel(
    "Scripts/Data_Vis/Section_4/Table_S13_bench_gene_shap_btwn_clusters_mwu_stats.xlsx",
    index=False)

################################################################################
# 4. Determine how many of the optimized features are ORFs absent in S288C or
#    intergenic SNPs.
################################################################################
# Feature to gene maps
map_snps = map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                                  sep="\t", header=None, names=["snp", "chr", "pos", "gene"], index_col=0)
map_orfs = pd.read_csv("Data/Peter_2018/final_map_orf_to_gene_CORRECTED.tsv",
                       sep="\t", index_col=0)

# Optimized feature sets for all single environment models
snp = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_snp.tsv", sep="\t", index_col=0)
pav = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_pav.tsv", sep="\t", index_col=0)
cnv = pd.read_csv(
    "Scripts/Data_Vis/Section_3/RF_optimized_gini_cnv.tsv", sep="\t", index_col=0)

# presence/absence of ORFs in S288C based on BLAST results
pav_info = pd.read_csv(
    "Data/Peter_2018/ORFs_pres_abs_with_S288C.csv", index_col=0)
pav_info = fix_orf_ids(pav_info.T).T

# For each env, count how many of the optimized features are intergenic/unknown ORFs
res = {}
for env in snp.columns[2:]:
    env_snp = snp[[env]].dropna()
    env_pav = pav[[env]].dropna()
    env_cnv = cnv[[env]].dropna()
    res[env] = {"SNP": [sum(map_snps.loc[env_snp.index, "gene"] == "intergenic"),  # number of intergenic snps
                        len(env_snp)],  # total number of snps
                "PAV": [sum(pd.merge(pav_info["S288C"], env_pav, how="right",
                                     left_index=True, right_index=True).S288C == 0),  # number of unknown orfs
                        len(env_pav)],  # total number of orfs
                "CNV": [sum(pd.merge(pav_info["S288C"], env_cnv, how="right",  # number of unknown orfs
                                     left_index=True, right_index=True).S288C == 0),
                        len(env_cnv)]}  # total number of orfs
    del env_snp, env_pav, env_cnv


res = pd.concat(
    {env: pd.DataFrame(
        data, index=["absent or intergenic", "Total"]).T for env, data in res.items()},
    axis=1
)
res.stack().T.to_csv(
    "Scripts/Data_Vis/Section_4/Number_unknown_or_intergenic_features_optimized_RF_models.csv")


# For the 5 target environments, what percentage of the optimized features are
# intergenic/ORFs absent in S288C (based on the BLAST analysis)?
res = res.stack().T
target_envs = ["YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL500", "YPDCUSO410MM",
               "YPDSODIUMMETAARSENITE"]
res_targ = res.loc[target_envs, :]
res_targ_percent = []
for env in target_envs:
    for data_type in ["SNP", "PAV", "CNV"]:
        df = res_targ.loc[env, data_type]
        res_targ_percent.append(
            [data_type, env, df["absent or intergenic"] / df["Total"]])

pd.DataFrame(res_targ_percent).groupby(0).describe().T
#               CNV       PAV       SNP
#   count  5.000000  5.000000  5.000000
#   mean   0.564063  0.676562  0.387569
#   std    0.184415  0.048349  0.026600
#   min    0.250000  0.601562  0.345000
#   25%    0.562500  0.656250  0.381167
#   50%    0.617188  0.695312  0.395333
#   75%    0.687500  0.710938  0.402344
#   max    0.703125  0.718750  0.414000
