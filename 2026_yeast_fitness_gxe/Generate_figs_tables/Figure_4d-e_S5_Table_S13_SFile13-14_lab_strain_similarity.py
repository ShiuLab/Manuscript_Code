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
import gc
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
geno = geno.astype(int)

# Add S288C SNP genotypes as a row of 0s (homozygous for the reference allele)
geno = geno.T
geno["S288C"] = 0
geno = geno.T
# geno.to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/geno_012_with_S288C_v2.csv") # S12 File
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
# eu_dist_snp = pd.read_csv("Scripts/Data_Vis/Section_4/genetic_distance_snp012_euclidean_to_S288C.csv",
#                           index_col=0)
# eu_dist_pav = pd.read_csv("Scripts/Data_Vis/Section_4/genetic_distance_pav_euclidean_to_S288C.csv",
#                           index_col=0)
eu_dist_snp = pd.read_excel('Scripts/Data_Vis/S13_File.xlsx', engine='openpyxl')
eu_dist_pav = pd.read_excel('Scripts/Data_Vis/S14_File.xlsx', engine='openpyxl')
eu_dist_snp.set_index("ID", inplace=True)
eu_dist_pav.set_index(eu_dist_pav.columns[0], inplace=True)

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
    print(k)
    kmeans_snp = KMeans(n_clusters=k, random_state=42, n_init=10).fit(snp_train)
    kmeans_pav = KMeans(n_clusters=k, random_state=42, n_init=10).fit(pav_train)
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
kmeans_snp = KMeans(n_clusters=6, random_state=42, n_init=10).fit(snp_train)
kmeans_pav = KMeans(n_clusters=4, random_state=42, n_init=10).fit(pav_train)

# Visualize the clusters using PCA
pca = PCA(n_components=2)
pca_snp = pca.fit(snp_train)
pca_snp_df = pca_snp.transform(snp_train)
pca_snp_df = pd.DataFrame(pca_snp_df, index=snp_train.index, columns=["PC1", "PC2"])
# variance explained by each component
vexp_pca_snp = pca_snp.explained_variance_ratio_

pca_p = PCA(n_components=2)
pca_pav = pca_p.fit(pav_train)
pca_pav_df = pca_p.transform(pav_train)
pca_pav_df = pd.DataFrame(pca_pav_df, index=pav_train.index, columns=["PC1", "PC2"])
vexp_pca_pav = pca_pav.explained_variance_ratio_

# SNP plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
scatter1 = axes[0].scatter(pca_snp_df[["PC1"]], pca_snp_df[["PC2"]], c=kmeans_snp.labels_,
                           label=kmeans_snp.labels_, cmap='tab10', s=60, alpha=0.7)
axes[0].set_title('K-means Clusters (SNP Features)')
axes[0].set_xlabel(f'PC1 ({vexp_pca_snp[0]:.2%} variance explained)')
axes[0].set_ylabel(f'PC2 ({vexp_pca_snp[1]:.2%} variance explained)')
axes[0].scatter(pca_snp_df.loc["S288C", "PC1"],
                pca_snp_df.loc["S288C", "PC2"],
                color='red', s=100, edgecolor='black', label='S288C')
axes[0].scatter(pca_snp_df.loc["SACE_GAV", "PC1"],
                pca_snp_df.loc["SACE_GAV", "PC2"],
                color='blue', s=100, edgecolor='black', label='W303')
legend1 = axes[0].legend(*scatter1.legend_elements(),
                         loc="upper right", title="Clusters", fontsize=6)
axes[0].add_artist(legend1)

# PAV plot
scatter2 = axes[1].scatter(pca_pav_df[["PC1"]], pca_pav_df[["PC2"]], c=kmeans_pav.labels_,
                           label=kmeans_pav.labels_, cmap='tab10', s=60, alpha=0.7)
axes[1].set_title('K-means Clusters (PAV Features)')
axes[1].set_xlabel(f'PC1 ({vexp_pca_pav[0]:.2%} variance explained)')
axes[1].set_ylabel(f'PC2 ({vexp_pca_pav[1]:.2%} variance explained)')
axes[1].scatter(pca_pav_df.loc["S288C", "PC1"],
                pca_pav_df.loc["S288C", "PC2"],
                color='red', s=100, edgecolor='black', label='S288C')
axes[1].scatter(pca_pav_df.loc["SACE_GAV", "PC1"],
                pca_pav_df.loc["SACE_GAV", "PC2"],
                color='blue', s=100, edgecolor='black', label='W303')
legend2 = axes[1].legend(*scatter2.legend_elements(),
                         loc="upper right", title="Clusters", fontsize=6)
axes[1].add_artist(legend2)
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_4/Figure_4b_S5_PCA_snp_or_pav_genetic_distance_v2.pdf",
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

# Identify the clusters most distinct to the cluster in which S288C and W303 (SACE_GAV) are in
cluster_assignments_snp = pd.DataFrame(kmeans_snp.labels_, index=kmeans_snp.feature_names_in_,
                                       columns=["Cluster"])  # S288C is in cluster 0; W303 is in cluster 0
cluster_assignments_pav = pd.DataFrame(kmeans_pav.labels_, index=kmeans_pav.feature_names_in_,
                                       columns=["Cluster"])  # S288C is in cluster 1; W303 is in cluster 2

# Re-assign cluster numbers so that S288C is in cluster 0
cluster_assignments_snp["New_Cluster"] = cluster_assignments_snp["Cluster"].replace(
    {0: 0, 1: 3, 2: 5, 3: 4, 4: 2, 5: 1})
cluster_assignments_pav["New_Cluster"] = cluster_assignments_pav["Cluster"].replace(
    {0: 3, 1: 0, 2: 2, 3: 1})

# cluster 5 (originally 2) is the most distinct to the S288C cluster 0 (originally cluster 0)
dist_snp.loc[0, :].idxmax() # 2
s288c_distinct_cluster_snp = 5
w303_distinct_cluster_snp = 5
# cluster 3 (originally 0) is the most distinct to the S288C cluster 0 (originally cluster 1)
dist_pav.loc[1, :].idxmax() # 0
s288c_distinct_cluster_pav = 3
dist_pav.loc[2, :].idxmax() # 1 (w303 is originally cluster 2); 1 is now cluster 0
w303_distinct_cluster_pav = 0

# assign dist_snp and dist_pav new cluster numbers
dist_snp.index = dist_snp.index.map({0: "0", 1: "3", 2: "5", 3: "4", 4: "2", 5: "1"})
dist_snp.columns = dist_snp.columns.map({0: "0", 1: "3", 2: "5", 3: "4", 4: "2", 5: "1"})
dist_pav.index = dist_pav.index.map({0: "3", 1: "0", 2: "2", 3: "1"})
dist_pav.columns = dist_pav.columns.map({0: "3", 1: "0", 2: "2", 3: "1"})

# sanity check
dist_snp.loc['0', :].idxmax() # '5' (most distinct to S288C and W303 cluster)
dist_pav.loc['0', :].idxmax() # '3' (most distinct to S288C cluster)
dist_pav.loc['2', :].idxmax() # '0' (most distinct to W303 cluster)

################################################################################
# 3. Compare the fitness distributions and shap values between the cluster
#    containing the lab strains and the most distinct cluster to it. (Figs. 4E,
#    S5B; Table S13)
################################################################################
# Benchmark genes validated in S288C or W303
ben_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/benomyl_phenotype_annotations_sensitive_genes.txt", sep="\t")
caf_meta = dt.fread(
    "Data/SGD_Experiment_Genes/caffeine_phenotype_annotations_sensitive_genes.txt", sep="\t").to_pandas()
cu_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/copperII_sulfate_phenotype_annotations_sensitive_genes.txt", sep="\t")
sma_meta = pd.read_csv(
    "Data/SGD_Experiment_Genes/sodium_arsenite_phenotype_annotations_sensitive_genes.txt", sep="\t")

ben_s288c = ben_meta.loc[ben_meta["Strain Background"]=="S288C", "Gene Systematic Name"].unique()
caf_s288c = caf_meta.loc[caf_meta["Strain Background"]=="S288C", "Gene Systematic Name"].unique()
cu_s288c = cu_meta.loc[cu_meta["Strain Background"]=="S288C", "Gene Systematic Name"].unique()
sma_s288c = sma_meta.loc[sma_meta["Strain Background"]=="S288C", "Gene Systematic Name"].unique()

ben_w303 = ben_meta.loc[ben_meta["Strain Background"]=="W303", "Gene Systematic Name"].unique()
caf_w303 = caf_meta.loc[caf_meta["Strain Background"]=="W303", "Gene Systematic Name"].unique()
cu_w303 = cu_meta.loc[cu_meta["Strain Background"]=="W303", "Gene Systematic Name"].unique()
sma_w303 = sma_meta.loc[sma_meta["Strain Background"]=="W303", "Gene Systematic Name"].unique()

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


snp_clust_map = {0: 0, 1: 3, 2: 5, 3: 4, 4: 2, 5: 1}  # {original: new}
pav_clust_map = {0: 3, 1: 0, 2: 2, 3: 1}
invert_snp_clust_map = {v: k for k, v in snp_clust_map.items()}
invert_pav_clust_map = {v: k for k, v in pav_clust_map.items()}


def benchmark_shap_clust_comparisons(env, bench_list, mtype, vtype, lab_strain, method, median, mwu_res):
    """Compare the distribution of benchmark gene SHAP values between clusters
    using Mann-Whitney U test. Clusters of isolates were determined based on the
    SNP- or PAV-derived genetic distances (Euclidean).
    
    Args:
        env (str): Environment name
        bench_list (list): List of benchmark genes validated in S288C or W303 for the environment
        mtype (str): Model type ("fs" for optimized RF models, "baseline" for complete RF models)
        vtype (str): Variant type ("snp" or "pav")
        lab_strain (str): Lab strain name ("S288C" or "W303")
        method (str): whether to take the absolute value of the SHAP values ("abs"),
            to keep the original SHAP values ("raw"), to keep only the positive
            SHAP values ("pos"), or to keep only the negative SHAP values ("neg")
        median (bool): whether to use the median SHAP value per gene (True) or
            the local SHAP values (False) for the Mann-Whitney U test
        mwu_res (dict): dictionary to store the Mann-Whitney U test results
    
    Returns: mwu_res (dict): updated dictionary with the Mann-Whitney U test
        results for the environment, variant type, lab strain, and "method".
    """
    d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"
    if vtype == "snp":
        # SNPs, read in the feature name map
        feat_map_shap = dt.fread(glob(
            f"{d}/SNP/{mtype}/SHAP_values_sorted_average_{env}_snp_rf_{mtype}_*_with_actual_feature_names_*")[0]).to_pandas()
        # set column names
        feat_map_shap.columns = feat_map_shap.iloc[0, :]
        feat_map_shap = feat_map_shap.iloc[1:, :]  # remove first row
        feat_map_shap.set_index(feat_map_shap.columns[1], inplace=True)
        feat_map_shap = feat_map_shap.drop(columns=0).to_dict()
        # read in shap values from the optimized RF models
        shap = dt.fread(glob(
            f"{d}/SNP/{mtype}/SHAP_values_sorted_{env}_snp_rf_{mtype}_*")[0]).to_pandas()
        shap.set_index("ID", inplace=True)
        shap.columns = shap.columns.map(feat_map_shap["actual_feature"])
    if vtype == "pav":
        # read in shap values from the optimized RF models
        shap = dt.fread(glob(
            f"{d}/PAV/{mtype}/SHAP_values_sorted_{env}_pav_rf_{mtype}_*")[0]).to_pandas()
        shap.set_index("ID", inplace=True)
        shap = fix_orf_ids(shap)
    
    # drop the features with zero shap values across all isolates
    print("shape before dropping zero columns", shap.shape)
    shap = shap.loc[:, ~shap.eq(0).all(axis=0)]
    print("shape after dropping zero columns", shap.shape)
    
    # subset the benchmark genes validated in s288c or w303 for an environment
    print(f"Number of {lab_strain} benchmark genes in list:", len(bench_list))
    if vtype == "snp":
        labstrain_shap = shap.loc[
            :, shap.columns.isin(map_snps.loc[map_snps.gene.isin(bench_list)].index)]
        n_bench = len(labstrain_shap.columns.map(map_snps.gene).unique())
        print(f"Number of {lab_strain} benchmark genes in SNP shap features: {n_bench}")
    if vtype == "pav":
        labstrain_shap = shap.loc[
            :, shap.columns.isin(map_orfs.loc[map_orfs.gene.isin(bench_list)].index)]
        n_bench = len(labstrain_shap.columns.map(map_orfs.gene).unique())
        print(f"Number of {lab_strain} benchmark genes in PAV shap features: {n_bench}" )
    print(f"Number of features in the lab strain SHAP subset:", labstrain_shap.shape[1])
    # print(labstrain_shap.head())
    # method for dealing with SHAP values
    if method == "abs":
        labstrain_shap = labstrain_shap.abs()  # absolute shap values
    elif method == "pos":
        labstrain_shap = labstrain_shap.mask(labstrain_shap <= 0, None)  # keep only positive shap values
    elif method == "neg":
        labstrain_shap = labstrain_shap.mask(labstrain_shap >= 0, None)  # keep only negative shap values
    else:
        pass  # keep the original shap values
    
    # plot the SHAP value distributions and calculate Mann-Whitney U test statistics between clusters
    if labstrain_shap.shape[1] != 0:
        # insert the updated isolate cluster assignments ("New_Cluster")
        if vtype == "snp":
            labstrain_shap.insert(
                0, "Cluster", cluster_assignments_snp.loc[labstrain_shap.index, "New_Cluster"])
        if vtype == "pav":
            labstrain_shap.insert(
                0, "Cluster", cluster_assignments_pav.loc[labstrain_shap.index, "New_Cluster"])
        
        labstrain_shap_box = labstrain_shap.reset_index().melt(
            id_vars=["Cluster", "ID"], var_name="Feature", value_name="SHAP") # reshape to long format for plotting
        labstrain_shap_box.dropna(inplace=True) # drop masked values
        
        # x1000 and log the values for better visualization
        if method not in ["neg", "raw"]:
            labstrain_shap_box["SHAP_log"] = np.log10(labstrain_shap_box["SHAP"] * 1000)
            
            # plot the distributions of the benchmark gene shap values
            sns.violinplot(data=labstrain_shap_box, x="Cluster", y="SHAP_log",
                           hue="Cluster", cut=0, palette="tab10")
            plt.ylabel("log10(SHAP Value x 1000)")
        else:
            labstrain_shap_box["SHAP_log"] = np.log10(labstrain_shap_box["SHAP"].abs() * 1000)
            sns.violinplot(data=labstrain_shap_box, x="Cluster", y="SHAP",
                           hue="Cluster", cut=0, palette="tab10")
            plt.ylabel("SHAP Value")
        
        plt.title(f"{env} {vtype.upper()} {lab_strain} Benchmark Gene SHAP Values")
        plt.xlabel("Cluster")
        if median:
            plt.savefig(
                f"Scripts/Data_Vis/Section_4/Figure_4e_{env}_{vtype}_rf_{mtype}_{lab_strain}_bench_gene_median_{method}_shap_violin.pdf")
        else:
            plt.savefig(
                f"Scripts/Data_Vis/Section_4/Figure_S5_{env}_{vtype}_rf_{mtype}_{lab_strain}_bench_gene_local_{method}_shap_violin.pdf")
        plt.close()
        
        if lab_strain == "S288C":
            if median:
                clustlab_shap = labstrain_shap.loc[labstrain_shap.Cluster == 0].drop(columns=["Cluster"])
                n_iso_clustlab = clustlab_shap.shape[0]
                n_fea_clustlab = clustlab_shap.shape[1]
            else:
                # subset the cluster 0 shap values (S288C is in cluster "0" in dist_snp and dist_pav)
                clustlab_shap = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0].drop(columns=["Cluster", "ID", "Feature", "SHAP_log"])
                n_iso_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0]["ID"].nunique()
                n_fea_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0]["Feature"].nunique()
        if lab_strain == "W303":
            if vtype == "snp": # W303 is in cluster "0"
                if median:
                    clustlab_shap = labstrain_shap.loc[labstrain_shap.Cluster == 0].drop(columns=["Cluster"])
                    n_iso_clustlab = clustlab_shap.shape[0]
                    n_fea_clustlab = clustlab_shap.shape[1]
                else:
                    clustlab_shap = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0].drop(columns=["Cluster", "ID", "Feature", "SHAP_log"])
                    n_iso_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0]["ID"].nunique()
                    n_fea_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 0]["Feature"].nunique()
            if vtype == "pav": # W303 is in cluster "2"
                if median:
                    clustlab_shap = labstrain_shap.loc[labstrain_shap.Cluster == 2].drop(columns=["Cluster"])
                    n_iso_clustlab = clustlab_shap.shape[0]
                    n_fea_clustlab = clustlab_shap.shape[1]
                else:
                    clustlab_shap = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 2].drop(columns=["Cluster", "ID", "Feature", "SHAP_log"])
                    n_iso_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 2]["ID"].nunique()
                    n_fea_clustlab = labstrain_shap_box.loc[labstrain_shap_box.Cluster == 2]["Feature"].nunique()
        
        # calculate the median shap value per gene per isolate
        if median:
            if vtype == "snp":
                clustlab_shap_med = clustlab_shap.T.groupby(
                    map_snps.loc[clustlab_shap.columns].gene).median() # median absolute SHAP if method = "abs"
            if vtype == "pav":
                clustlab_shap_med = clustlab_shap.T.groupby(
                    map_orfs.loc[clustlab_shap.columns].gene).median()
        
        # loop through the non-laboratory strain clusters
        for clust in sorted(labstrain_shap.Cluster.unique().tolist()):
            if (clust == 0 and lab_strain == "S288C"):
                lab_clust = 0
                continue
            elif (clust == 0 and vtype == "pav" and lab_strain == "W303"):
                lab_clust = 2
            elif (clust == 2 and vtype == "pav" and lab_strain == "W303"):
                lab_clust = 2
                continue
            elif (clust == 0 and vtype == "snp" and lab_strain == "W303"):
                lab_clust = 0
                continue
            print(f"Comparing cluster {lab_clust} to cluster {clust} for {env} {vtype.upper()} {lab_strain}")
            
            # subset the second cluster's shap values
            if median:
                clust_shap = labstrain_shap.loc[labstrain_shap.Cluster == clust].drop(columns=["Cluster"])
                n_iso_clust = clust_shap.shape[0]
                n_fea_clust = clust_shap.shape[1]
            else:
                clust_shap = labstrain_shap_box.loc[labstrain_shap_box.Cluster == clust].drop(columns=["Cluster", "ID", "Feature", "SHAP_log"])
                n_iso_clust = labstrain_shap_box.loc[labstrain_shap_box.Cluster == clust]["ID"].nunique()
                n_fea_clust = labstrain_shap_box.loc[labstrain_shap_box.Cluster == clust]["Feature"].nunique()
            print("clust_shap", clust_shap.shape)
            
            # calculate the median absolute shap value per gene per isolate
            if median:
                if vtype == "snp":
                    clust_shap_med = clust_shap.T.groupby(
                        map_snps.loc[clust_shap.columns].gene).median() # rows: genes, cols: isolates
                if vtype == "pav":
                    clust_shap_med = clust_shap.T.groupby(
                    map_orfs.loc[clust_shap.columns].gene).median()
                print("clust_shap_med", clust_shap_med.shape)
            
            # Are clustlab shap values significantly greater than the other cluster's shap?
            try:
                if median:
                    # use the median absolute SHAP values (when method = "abs")
                    x = clustlab_shap_med.values.flatten()
                    y = clust_shap_med.values.flatten()
                    u, pval = mannwhitneyu(x, y, alternative="greater")
                else: # use the local SHAP values (not medians)
                    u, pval = mannwhitneyu(clustlab_shap, clust_shap, alternative="greater")
                if vtype == "snp":
                    if median:
                        mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                            [u, pval, n_bench, labstrain_shap.shape[1]-1, dist_snp.loc[str(lab_clust), str(clust)]] +\
                            [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust]
                    else:
                        mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                            [u[0], pval[0], n_bench, labstrain_shap.shape[1]-1, dist_snp.loc[str(lab_clust), str(clust)]] +\
                            [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust] #+\
                            # clustlab_shap.describe().values.flatten().tolist() +\
                            # clust_shap.describe().values.flatten().tolist()
                if vtype == "pav":
                    if median:
                        mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                            [u, pval, n_bench, labstrain_shap.shape[1]-1, dist_pav.loc[str(lab_clust), str(clust)]] +\
                            [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust]
                    else:
                        mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                            [u[0], pval[0], n_bench, labstrain_shap.shape[1]-1, dist_pav.loc[str(lab_clust), str(clust)]] +\
                            [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust] #+\
                        # clustlab_shap.describe().values.flatten().tolist() +\
                        # clust_shap.describe().values.flatten().tolist()
            except ValueError:
                if vtype == "snp":
                    mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                        [np.nan, np.nan, n_bench, labstrain_shap.shape[1]-1, dist_snp.loc[str(lab_clust), str(clust)]] +\
                        [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust]
                if vtype == "pav":
                    mwu_res[(vtype, env, mtype, lab_strain, method, f"clust{lab_clust}_vs_clust{clust}")] = \
                        [np.nan, np.nan, n_bench, labstrain_shap.shape[1]-1, dist_pav.loc[str(lab_clust), str(clust)]] +\
                        [n_iso_clustlab, n_iso_clust, n_fea_clustlab, n_fea_clust]
            del clust_shap #, clust_shap_med
        del labstrain_shap_box, clustlab_shap, #clustlab_shap_med
    del shap, labstrain_shap
    gc.collect()
    gc.collect()
    return mwu_res


# benchmark_shap_clust_comparisons(env, bench_list, mtype, vtype, lab_strain, method, mwu_res)
for median in [True]: #, False]:
    mwu_res = {}
    for mtype in ["fs", "baseline"]:
        for method in ["abs", "raw", "pos", "neg"]:
            for vtype in ["snp", "pav"]:
                for env, bench_list in zip(["YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"],
                            [ben_s288c, caf_s288c, caf_s288c, cu_s288c, sma_s288c]):
                    print("Configs:", env, mtype, vtype, "S288C", method)
                    mwu_res = benchmark_shap_clust_comparisons(env, bench_list, mtype, vtype, "S288C", method, median, mwu_res=mwu_res)
                
                for env, bench_list in zip(["YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"],
                            [ben_w303, caf_w303, caf_w303, cu_w303, sma_w303]):
                    print("Configs:", env, mtype, vtype, "W303", method)
                    mwu_res = benchmark_shap_clust_comparisons(env, bench_list, mtype, vtype, "W303", method, median, mwu_res=mwu_res)
    
    
    mwu_res_df = pd.DataFrame.from_dict(mwu_res, orient="index")
    mwu_res_df.index = pd.MultiIndex.from_tuples(mwu_res_df.index)
    mwu_res_df.index.names = ["Variant Type", "Environment", "RF Model", "Lab Strain",
                            "SHAP Value Type", "Comparison"]
    mwu_res_df.columns = ["U", "pval", "num_benchmark_genes", "num_benchmark_features",
                        "distance_between_clusters", "num_isolates_clustlab",
                        "num_isolates_other_clust", "num_features_clustlab",
                        "num_features_other_clust"] #, "clustlab_shap_count",
                        # "clustlab_shap_mean", "clustlab_shap_std",
                        # "clustlab_shap_min", "clustlab_shap_25%",
                        # "clustlab_shap_50%", "clustlab_shap_75%",
                        # "clustlab_shap_max", "other_clust_shap_count",
                        # "other_clust_shap_mean", "other_clust_shap_std",
                        # "other_clust_shap_min", "other_clust_shap_25%",
                        # "other_clust_shap_50%", "other_clust_shap_75%",
                        # "other_clust_shap_max"]
    if median:
        # mwu_res_df["difference_median_shap_btwn_clusters"] = mwu_res_df["s288c_clust_shap_50%"] - \
        #     mwu_res_df["other_clust_shap_50%"]
        mwu_res_df.reset_index().to_excel(
            "Scripts/Data_Vis/Section_4/Table_S13_bench_gene_median_shap_btwn_clusters_mwu_stats_with_w303.xlsx",
            index=False) # these stats were calculated with median {method} SHAP values per gene
    else:
        mwu_res_df.reset_index().to_excel(
            "Scripts/Data_Vis/Section_4/Table_S13_bench_gene_local_shap_btwn_clusters_mwu_stats_with_w303.xlsx",
            index=False) # these stats were calculated with NON-MEDIAN SHAP values per gene


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
