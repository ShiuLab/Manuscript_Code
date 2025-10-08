#!/usr/bin/env python3

################################################################################
# In this script, I do the following:
# 1. Combine the average Gini importances for all environments (Table S9 & Supplementary file 7)
# 2. Combine the average absolute SHAP values for all environments (Table S9 & Supplementary file 7)
# 3. Add GO terms, pathways, and benchmark gene info to Table S9 and Supplementary file 7
#    Also add the annotations to the feature to gene mapping files (Supplementary files 8 & 9)
# 4. Calculate spearman's rho correlations between Gini and SHAP values for each environment (Table S5)
################################################################################

import os
import re
import pandas as pd
import datatable as dt
import numpy as np
from scipy.stats import spearmanr
from tqdm import tqdm
from scipy.stats import ks_2samp

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

# read feature to gene map files; will be used throughout the script
map_snps = pd.read_csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                       sep="\t", header=None, names=["snp", "chr", "pos", "gene"])
map_orfs = pd.read_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed.tsv", sep="\t")
map_snps.merge(map_orfs, how="inner", on="gene").gene.nunique()  # 5370 genes
map_snps["gene_with_intergenic"] = map_snps.apply(
    lambda row: f"intergenic//{row['snp']}" if row["gene"] == "intergenic" else row["gene"], axis=1)

################################################################################
# 1. Combine the average Gini importances for all environments
################################################################################
# paths to feature importance score files of optimized and complete RF models
# SNP optimized and complete RF modeling results
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/fs"
snp_rf_res = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")
snp_optimized_files = [os.path.join(dir, f"{x}_imp_with_actual_feature_names_09102025")
                       for x in snp_rf_res['ID']]  # optimized RF models
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/baseline"
snp_complete_files = [os.path.join(
    dir, f"{x}_rf_baseline_imp_with_actual_feature_names_09102025") for x in mapping.keys()]  # complete RF models

# ORF pres/abs optimized and complete RF modeling results
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/fs"
pav_rf_res = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")
pav_fs_files = [os.path.join(dir, f"{x}_imp") for x in pav_rf_res['ID']]
cnv_rf_res = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")  # CNV FS results
cnv_fs_files = [os.path.join(dir, f"{x}_imp") for x in cnv_rf_res['ID']]
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/ORF_yeast_RF_results/baseline"
pav_complete_files = [os.path.join(
    dir, f"{x}_pav_baseline_imp") for x in mapping.keys()]
cnv_complete_files = [os.path.join(
    dir, f"{x}_cnv_baseline_imp") for x in mapping.keys()]


def combine_imp_indiv(imp_files, map=map_snps, dtype="snp", save="", mapping=mapping):
    # combine gini importance for all envs per data type
    for i, env in enumerate(mapping.keys()):
        print(env)
        # Read gini importance file
        file = [f for f in imp_files if env in f]
        print(len(file))  # should be 1
        imp = dt.fread(file[0]).to_pandas()
        imp.set_index(imp.iloc[:, 0], inplace=True)  # feature names as index
        imp = imp.loc[:, "mean_imp"]  # use mean gini importances
        imp.rename(env, inplace=True)
        imp = pd.DataFrame(imp)
        if dtype != "snp":
            imp.index = imp.apply(lambda x: re.sub(
                "^X", "", x.name), axis=1)  # rename index
            imp.index = imp.apply(lambda x: re.sub("\.", "-", x.name), axis=1)
        if i == 0:
            merged = imp.copy(deep=True)
        else:
            merged = pd.concat([merged, imp], axis=1,
                               ignore_index=False)  # add to dictionary
        del imp
    print(merged.shape)
    # map to genes
    if dtype == "snp":
        merged = map_snps[["snp", "gene", "gene_with_intergenic"]].merge(
            merged, how="right", left_on="snp", right_index=True)
    else:
        merged = map_orfs[["orf", "gene"]].merge(
            merged, how="right", left_on="orf", right_index=True)
    merged.to_csv(save, sep="\t", index=False)
    return merged


combine_imp_indiv(snp_optimized_files, map=map_snps, dtype="snp",
                  save="Scripts/Data_Vis/Section_3/RF_optimized_gini_snp.tsv")
combine_imp_indiv(pav_fs_files, map=map_orfs, dtype="pav",
                  save="Scripts/Data_Vis/Section_3/RF_optimized_gini_pav.tsv")
combine_imp_indiv(cnv_fs_files, map=map_orfs, dtype="cnv",
                  save="Scripts/Data_Vis/Section_3/RF_optimized_gini_cnv.tsv")
combine_imp_indiv(snp_complete_files, map=map_snps, dtype="snp",
                  save="Scripts/Data_Vis/Section_3/RF_complete_gini_snp.tsv")
combine_imp_indiv(pav_complete_files, map=map_orfs, dtype="pav",
                  save="Scripts/Data_Vis/Section_3/RF_complete_gini_pav.tsv")
combine_imp_indiv(cnv_complete_files, map=map_orfs, dtype="cnv",
                  save="Scripts/Data_Vis/Section_3/RF_complete_gini_cnv.tsv")

################################################################################
# 2. Combine the average absolute SHAP values for all environments
################################################################################
# paths to feature average absolute SHAP value files
dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP"
snp_shap_optimized_files = [f"{dir}/SNP/fs/{file}" for file in os.listdir(
    dir + "/SNP/fs") if file.endswith("_with_actual_feature_names_09102025")]
snp_shap_complete_files = [f"{dir}/SNP/baseline/{file}" for file in os.listdir(
    dir + "/SNP/baseline") if file.endswith("_with_actual_feature_names_09102025")]
pav_shap_optimized_files = [f"{dir}/PAV/fs/{file}" for file in os.listdir(
    dir + "/PAV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
pav_shap_complete_files = [f"{dir}/PAV/baseline/{file}" for file in os.listdir(
    dir + "/PAV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_optimized_files = [f"{dir}/CNV/fs/{file}" for file in os.listdir(
    dir + "/CNV/fs") if file.startswith("SHAP_values_sorted_average_Y")]
cnv_shap_complete_files = [f"{dir}/CNV/baseline/{file}" for file in os.listdir(
    dir + "/CNV/baseline") if file.startswith("SHAP_values_sorted_average_Y")]


def combine_shap_indiv(shap_files, map=map_snps, merged={}, dtype="snp", save="", mapping=mapping):
    # combine shap values for all envs per data type
    for i, env in enumerate(mapping.keys()):
        print(env)
        # Read SHAP file
        file = [f for f in shap_files if env in f]
        print(len(file))  # should be 1
        shap = dt.fread(file[0]).to_pandas()
        shap.set_index(shap.iloc[:, 0], inplace=True)
        if dtype == "snp":
            shap = shap.iloc[:, 2:]  # only the numeric columns
            shap.rename(columns={"C2": env}, inplace=True)
        if dtype != "snp":
            shap = shap.iloc[:, 1:]
            shap.rename(columns={"C1": env}, inplace=True)
            shap.index = shap.apply(lambda x: re.sub(
                "^X", "", x.name), axis=1)  # rename index
            shap.index = shap.apply(
                lambda x: re.sub("\.", "-", x.name), axis=1)
        if i == 0:
            merged = shap.copy(deep=True)
        else:
            merged = pd.concat([merged, shap], axis=1,
                               ignore_index=False)  # add to dictionary
        del shap
    print(merged.shape)
    # map to genes
    if dtype == "snp":
        merged = map_snps[["snp", "gene", "gene_with_intergenic"]].merge(
            merged, how="right", left_on="snp", right_index=True)
        merged.iloc[1:, :].to_csv(save, sep="\t", index=False)
    else:
        merged = map_orfs[["orf", "gene"]].merge(
            merged, how="right", left_on="orf", right_index=True)
        merged.to_csv(save, sep="\t", index=False)
    return merged


combine_shap_indiv(snp_shap_optimized_files, map=map_snps, dtype="snp",
                   save="Scripts/Data_Vis/Section_3/RF_optimized_shap_snp.tsv")
combine_shap_indiv(pav_shap_optimized_files, map=map_orfs, dtype="pav",
                   save="Scripts/Data_Vis/Section_3/RF_optimized_shap_pav.tsv")
combine_shap_indiv(cnv_shap_optimized_files, map=map_orfs, dtype="cnv",
                   save="Scripts/Data_Vis/Section_3/RF_optimized_shap_cnv.tsv")
combine_shap_indiv(snp_shap_complete_files, map=map_snps, dtype="snp",
                   save="Scripts/Data_Vis/Section_3/RF_complete_shap_snp.tsv")
combine_shap_indiv(pav_shap_complete_files, map=map_orfs, dtype="pav",
                   save="Scripts/Data_Vis/Section_3/RF_complete_shap_pav.tsv")
combine_shap_indiv(cnv_shap_complete_files, map=map_orfs, dtype="cnv",
                   save="Scripts/Data_Vis/Section_3/RF_complete_shap_cnv.tsv")

################################################################################
# 3. Add GO terms, pathways, and benchmark gene info to Table S9 and Supplementary file 7
################################################################################
# Read in the annotation files
go = pd.read_csv("Data/yeast_GO/sgd_GO_BP_no_info.txt", sep="\t")
pwy_orf = pd.read_csv(
    "Data/Peter_2018/ORFs_and_S288C_genes_pwy_all_CORRECTED_16_removed.csv")
pwy_snp = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_pathway_map.csv")
map_snps = pd.read_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark.tsv", sep="\t")
map_orfs = pd.read_csv(
    "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark.tsv", sep="\t")
curated_snp = pd.read_csv(
    "Data/SGD_Experiment_Genes/manually_curated_genes_snps.txt", sep="\t")
curated_orf = pd.read_csv(
    "Data/SGD_Experiment_Genes/manually_curated_genes_orfs.txt", sep="\t")

# Delimet all the annotations with "; " for each gene
go_agg = go.groupby("Gene")["GO.ID"].agg(lambda x: "; ".join(x)).reset_index()
pwy_orf = pwy_orf[["Accession.1", "Pathways.of.gene"]].astype(str)
pwy_orf_agg = pwy_orf.groupby("Accession.1")["Pathways.of.gene"].agg(
    lambda x: "; ".join(x)).reset_index()
pwy_snp = pwy_snp[["gene", "pathway"]].astype(str)
pwy_snp_agg = pwy_snp.groupby("gene")["pathway"].agg(
    lambda x: "; ".join(x)).reset_index()

# Order of the environments to save the columns as in the final tables
env_order = pd.read_csv(
    "Scripts/Data_Vis/Section_2/Figure_2_r2_test_v5_env_order.txt", sep="\t", header=None)
inverted_mapping = {v: k for k, v in mapping.items()}
env_order["Envs"] = env_order[0].map(inverted_mapping)

# Read in the combined Gini and SHAP files and add annotations
for data_type in ["snp", "pav", "cnv"]:
    for model_type in ["optimized", "complete"]:
        for imp_type in ["shap", "gini"]:
            df = dt.fread(
                f"Scripts/Data_Vis/Section_3/RF_{model_type}_{imp_type}_{data_type}.tsv").to_pandas()
            # add GO terms
            df = df.merge(go_agg, how="left", left_on="gene", right_on="Gene")
            # add pathways and benchmark gene info
            if data_type != "snp":
                df = df.merge(pwy_orf_agg, how="left",
                              left_on="gene", right_on="Accession.1")
                df = df.merge(map_orfs[["orf", "Benomyl", "Caffeine", "CuSO4",
                                        "Sodium_meta-arsenite"]], how="left",
                              left_on="orf", right_on="orf")
                # re-order the df columns
                df = df[["orf", "gene", "GO.ID", "Pathways.of.gene", "Benomyl",
                        "Caffeine", "CuSO4", "Sodium_meta-arsenite"] +
                        env_order["Envs"].tolist()]
                df.insert(8, "Curated", np.where(
                    df["gene"].isin(curated_orf), 1, 0))
            else:
                df = df.merge(pwy_snp_agg, how="left",
                              left_on="gene", right_on="gene")
                df = df.merge(map_snps[["snp", "Benomyl", "Caffeine", "CuSO4",
                                        "Sodium_meta-arsenite"]], how="left",
                              left_on="snp", right_on="snp")
                # re-order the df columns
                df = df[["snp", "gene", "GO.ID", "pathway", "Benomyl",
                         "Caffeine", "CuSO4", "Sodium_meta-arsenite"] +
                        env_order["Envs"].tolist()]
                df.insert(8, "Curated", np.where(
                    df["gene"].isin(curated_snp), 1, 0))
            #
            # Save the updated files
            df.to_csv(
                f"Scripts/Data_Vis/Section_3/RF_{model_type}_{imp_type}_{data_type}_with_annotations.tsv", sep="\t", index=False)
            del df


# Add annotations to the feature to gene mapping files (Supplementary data files 8 & 9)
map_snps.insert(8, "Curated", np.where(
    map_snps["gene"].isin(curated_snp), 1, 0))
map_orfs.insert(7, "Curated", np.where(
    map_orfs["gene"].isin(curated_orf), 1, 0))
map_snps["Curated"].sum()  # 0, so drop
map_snps.drop(columns=["Curated"], inplace=True)
map_orfs["Curated"].sum()  # 4, so drop
map_orfs.drop(columns=["Curated"], inplace=True)
map_snps = map_snps.merge(go_agg, how="left", left_on="gene", right_on="Gene")
map_snps = map_snps.merge(pwy_snp_agg, how="left",
                          left_on="gene", right_on="gene")
map_orfs = map_orfs.merge(go_agg, how="left", left_on="gene", right_on="Gene")
map_orfs = map_orfs.merge(pwy_orf_agg, how="left",
                          left_on="gene", right_on="Accession.1")
map_snps.to_csv(
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_expanded_benchmark_with_annotations.tsv", sep="\t", index=False)
map_orfs.to_csv(
    "Data/Peter_2018/final_map_orf_to_gene_CORRECTED_16_removed_expanded_benchmark_with_annotations.tsv", sep="\t", index=False)

map_snps.isna().sum()
# snp                         0
# chr                         0
# pos                         0
# gene                        0
# Benomyl                     0
# Caffeine                    0
# CuSO4                       0
# Sodium_meta-arsenite        0
# Gene                    49971
# GO.ID                   49971
# pathway                 45585
map_orfs.isna().sum()
# orf                       0
# gene                      0
# organism                  0
# Benomyl                   0
# Caffeine                  0
# CuSO4                     0
# Sodium_meta-arsenite      0
# Gene                    759
# GO.ID                   759
# Accession.1             351
# Pathways.of.gene        351

################################################################################
# 4. Calculate spearman's rho correlations between Gini and SHAP values for each environment.
################################################################################
res = [["Model Type", "Data Type", "Env", "rho", "pval", "NumShared"]]
for data_type in ["snp", "pav", "cnv"]:
    # Optimized RF models
    gini = dt.fread(
        f"Scripts/Data_Vis/Section_3/RF_optimized_gini_{data_type}.tsv").to_pandas()
    shap = dt.fread(
        f"Scripts/Data_Vis/Section_3/RF_optimized_shap_{data_type}.tsv").to_pandas()
    #
    # Complete RF models
    gini_comp = dt.fread(
        f"Scripts/Data_Vis/Section_3/RF_complete_gini_{data_type}.tsv").to_pandas()
    shap_comp = dt.fread(
        f"Scripts/Data_Vis/Section_3/RF_complete_shap_{data_type}.tsv").to_pandas()
    if data_type == "snp":
        gini.set_index("snp", inplace=True)
        shap.set_index("snp", inplace=True)
        gini_comp.set_index("snp", inplace=True)
        shap_comp.set_index("snp", inplace=True)
    else:
        gini.set_index("orf", inplace=True)
        shap.set_index("orf", inplace=True)
        gini_comp.set_index("orf", inplace=True)
        shap_comp.set_index("orf", inplace=True)
    #
    # Convert data to rank percentiles
    gini_env_rank_out = pd.DataFrame()  # just FS features
    shap_env_rank_out = pd.DataFrame()
    gini_comp_rank_out = pd.DataFrame()  # just baseline features
    shap_comp_rank_out = pd.DataFrame()
    for env in mapping.keys():
        # First rank the optimized RF model features
        gini_env = gini.loc[:, env]
        shap_env = shap.loc[:, env]
        # remove the extra features from the feature selection dataset
        gini_env.dropna(inplace=True)
        shap_env.dropna(inplace=True)
        # drop features with zero importance
        gini_env = gini_env[gini_env != 0]
        shap_env = shap_env[shap_env != 0]
        # calculate the rank percentiles (1= most important)
        gini_env_rank = gini_env.sort_values(ascending=False).rank(
            axis=0, method="average", numeric_only=True, pct=True)
        shap_env_rank = shap_env.sort_values(ascending=False).rank(
            axis=0, method="average", numeric_only=True, pct=True)
        gini_env_rank_out = pd.concat(
            [gini_env_rank_out, gini_env_rank], axis=1, ignore_index=False)
        shap_env_rank_out = pd.concat(
            [shap_env_rank_out, shap_env_rank], axis=1, ignore_index=False)
        #
        # Now rank the complete RF model features
        gini_comp_env = gini_comp.loc[:, env].dropna()
        shap_comp_env = shap_comp.loc[:, env].dropna()
        gini_comp_env = gini_comp_env[gini_comp_env != 0]
        shap_comp_env = shap_comp_env[shap_comp_env != 0]
        gini_comp_env_rank = gini_comp_env.sort_values(ascending=False).rank(
            axis=0, method="average", numeric_only=True, pct=True)
        shap_comp_env_rank = shap_comp_env.sort_values(ascending=False).rank(
            axis=0, method="average", numeric_only=True, pct=True)
        gini_comp_rank_out = pd.concat(
            [gini_comp_rank_out, gini_comp_env_rank], axis=1, ignore_index=False)
        shap_comp_rank_out = pd.concat(
            [shap_comp_rank_out, shap_comp_env_rank], axis=1, ignore_index=False)
        #
        # Calculate spearman's rho
        # First, for the optimized RF model features
        df = pd.concat([gini_env_rank_out.loc[:, env],
                        shap_env_rank_out.loc[:, env]], ignore_index=False, axis=1).dropna()
        rho = df.corr(method=lambda x, y: spearmanr(
            x, y, alternative="two-sided").statistic)
        pval = df.corr(method=lambda x, y: spearmanr(
            x, y, alternative="two-sided").pvalue)
        res.append(["optimized", data_type, env, rho.iloc[0, 1],
                   pval.iloc[0, 1], len(df)])
        # Second, for the complete RF model features
        df = pd.concat([gini_comp_rank_out.loc[:, env],
                        shap_comp_rank_out.loc[:, env]], ignore_index=False, axis=1).dropna()
        rho = df.corr(method=lambda x, y: spearmanr(
            x, y, alternative="two-sided").statistic)
        pval = df.corr(method=lambda x, y: spearmanr(
            x, y, alternative="two-sided").pvalue)
        res.append(["complete", data_type, env,
                   rho.iloc[0, 1], pval.iloc[0, 1], len(df)])


res = pd.DataFrame(res)
res.columns = res.iloc[0, :]
res = res.iloc[1:, :]
res.sort_values(by='rho', ascending=False, inplace=True)
res.to_csv("Scripts/Data_Vis/Section_3/Table_S5_gini_vs_shap_rank_per_corr.tsv",
           sep="\t", index=False)

# Is model performance correlated with the correlation between Gini and SHAP values?
rho = pd.read_csv(
    'Scripts/Data_Vis/Section_3/Table_S5_gini_vs_shap_rank_per_corr.tsv', sep='\t')
snp = pd.read_csv(
    'Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt', sep='\t', index_col='Y')
pav = pd.read_csv(
    'Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt', sep='\t', index_col='Y')
cnv = pd.read_csv(
    'Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt', sep='\t', index_col='Y')

r2_v_rho_snp = rho[(rho['Model Type'] == 'optimized') & (rho['Data Type'] == 'snp')].
merge(snp['r2_test'], right_index=True, left_on="Env")
# np.float64(-0.032360053556767396)
r2_v_rho_snp[['rho', 'r2_test']].corr().iloc[0, 1]

r2_v_rho_pav = rho[(rho['Model Type'] == 'optimized') & (rho['Data Type'] == 'pav')].
merge(pav['r2_test'], right_index=True, left_on="Env")
# np.float64(0.3300096647369571)
r2_v_rho_pav[['rho', 'r2_test']].corr().iloc[0, 1]

r2_v_rho_cnv = rho[(rho['Model Type'] == 'optimized') & (rho['Data Type'] == 'cnv')].
merge(cnv['r2_test'], right_index=True, left_on="Env")
# np.float64(0.287642642641329)
r2_v_rho_cnv[['rho', 'r2_test']].corr().iloc[0, 1]

r2_v_rho_snp = rho[(rho['Model Type'] == 'complete') & (rho['Data Type'] == 'snp')].
merge(snp['r2_test'], right_index=True, left_on="Env")
# np.float64(-0.370472133387726)
r2_v_rho_snp[['rho', 'r2_test']].corr().iloc[0, 1]

r2_v_rho_pav = rho[(rho['Model Type'] == 'complete') & (rho['Data Type'] == 'pav')].
merge(pav['r2_test'], right_index=True, left_on="Env")
# np.float64(0.2079527477389604)
r2_v_rho_pav[['rho', 'r2_test']].corr().iloc[0, 1]

r2_v_rho_cnv = rho[(rho['Model Type'] == 'complete') & (rho['Data Type'] == 'cnv')].
merge(cnv['r2_test'], right_index=True, left_on="Env")
# np.float64(0.3529860550102224)
r2_v_rho_cnv[['rho', 'r2_test']].corr().iloc[0, 1]
