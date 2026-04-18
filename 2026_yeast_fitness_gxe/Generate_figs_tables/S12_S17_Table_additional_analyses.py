#! /usr/bin/env python3
"""
Additional analysis to supplement predicted single genes and shap interactions

0. How many SGD benchmark genes did I use in the original analysis?
	- What percentage of SGD benchmark genes were represented by the optimized YPD Benomyl 500 ug/ml models?

1. Do SGD benchmark genes have greater SHAP values than non-benchmark genes?

2. To what extent do the SHAP values of SGD benchmark genes agree with Costanzo
   et al., 2021 differential relative fitness (DRF) values?
	- How many benchmark genes are recovered in different percentiles by SHAP or by DRF?
	- How does the agreement between SHAP and DRF values of benchmark genes compare to non-benchmark genes?
	- Based on the findings, how suitable is the Costanzo DRF for validating our models?

3. Analysis of SHAP interaction values from the Figure 4A models (S17 Table).
	a. Combine SNP, PAV, and CNV "important non-benchmark gene features + 
	   benchmark genes" feature tables to train models on the integrated features.
	b. Calculate the SHAP interaction values for SNP, PAV, CNV, SNP + PAV, SNP +
	   CNV, PAV + CNV, and SNP + PAV + CNV models.
	c. Estimate the correlation of local SHAP interaction values and epsilon
	   values of the overlapping gene pairs.
	d. Prepare feature table to predict SHAP interaction values.

4. Examine gene CNV values (1, 2, 3, etc) and fitness (S12 Table)
- How do the effects of CNV features (SHAP values) differ across environments?
	- A feature's CNVs across isolates--> compare to (pearson)--> local SHAP vals
		- calculate for all 625 training isolates combined
		- calculate for each SHAP cluster of isolates
	- Compare differences in pearson correlations across environments
	- Correlate CNV values with fitness values across isolates (pearson)
		- How do the correlations change with median absolute feature importance?

salloc -N 1 -c 10 --mem=80GB --time=8:00:00
conda activate /mnt/home/seguraab/miniconda3/envs/shap
/mnt/home/seguraab/miniconda3/envs/shap/bin/python3
"""

import os, gc, re
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# ------------------------------------------------------------------------------
# 0. How many SGD benchmark genes did I use in the original analysis?
# ------------------------------------------------------------------------------
# read SNP and ORF gene maps
map_snps = pd.read_excel("Scripts/Data_Vis/S8_File.xlsx", engine="openpyxl")
map_orfs = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl")
map_snps.shape # (118382, 8)
map_orfs.shape # (5902, 7)

map_snps["gene"] = map_snps["gene"].str.split(",") # split up the genes into a list
map_snps = map_snps.explode(column="gene", ignore_index=True) # each gene gets its own row
map_snps.shape # (119616, 8)

map_snps["benchmark"] = map_snps[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1)
map_orfs["benchmark"] = map_orfs[["Benomyl", "Caffeine", "CuSO4", "Sodium_meta-arsenite"]].sum(axis=1)

map_snps.loc[map_snps["benchmark"] > 0, "gene"].nunique() # 1196 genes
map_orfs.loc[map_orfs["benchmark"] > 0, "gene"].nunique() # 1136 genes
len(set(map_snps.loc[map_snps["benchmark"] > 0, "gene"]).union(set(map_orfs.loc[map_orfs["benchmark"] > 0, "gene"]))) # 1213 genes
len(set(map_snps.loc[map_snps["Benomyl"] > 0, "gene"]).union(set(map_orfs.loc[map_orfs["Benomyl"] > 0, "gene"]))) # 376 SGD benomyl benchmark genes that mapped to SNPs and/or ORFs
len(set(map_snps.loc[map_snps["Caffeine"] > 0, "gene"]).union(set(map_orfs.loc[map_orfs["Caffeine"] > 0, "gene"]))) # 735 SGD caffeine benchmark genes
len(set(map_snps.loc[map_snps["CuSO4"] > 0, "gene"]).union(set(map_orfs.loc[map_orfs["CuSO4"] > 0, "gene"]))) # 158 SGD CuSO4 benchmark genes
len(set(map_snps.loc[map_snps["Sodium_meta-arsenite"] > 0, "gene"]).union(set(map_orfs.loc[map_orfs["Sodium_meta-arsenite"] > 0, "gene"]))) # 280 SGD Sodium_meta-arsenite benchmark genes

# What percentage of SGD benchmark genes were represented by the optimized YPD Benomyl 500 ug/ml SNP model?
opt_snp = pd.read_excel("Scripts/Data_Vis/S9_Table.xlsx", engine="openpyxl", sheet_name="shap_snp")
opt_pav = pd.read_excel("Scripts/Data_Vis/S9_Table.xlsx", engine="openpyxl", sheet_name="shap_pav")
opt_cnv = pd.read_excel("Scripts/Data_Vis/S9_Table.xlsx", engine="openpyxl", sheet_name="shap_cnv")

ypdbeno500_opt_snp = opt_snp.loc[~opt_snp["YPDBENOMYL500"].isna(),
	["SNP Feature", "Gene", "Benomyl", "YPDBENOMYL500"]] # has 6000 rows (all features are unique)
ypdbeno500_opt_pav = opt_pav.loc[~opt_pav["YPDBENOMYL500"].isna(),
	["ORF Feature", "Gene", "Benomyl", "YPDBENOMYL500"]] # has 128 rows (all features are unique)
ypdbeno500_opt_cnv = opt_cnv.loc[~opt_cnv["YPDBENOMYL500"].isna(),
	["ORF Feature", "Gene", "Benomyl", "YPDBENOMYL500"]] # has 128 rows (all features are unique)

ypdbeno500_opt_snp.groupby("Gene")["Benomyl"].max().value_counts()
# Benomyl
# 0    2002 <- includes "intergenic" group
# 1     128 <-- regardless if the average absolute SHAP value is 0, there are 128 benchmarks in the optimized model features
map_snps.loc[map_snps.snp.isin(ypdbeno500_opt_snp["SNP Feature"]), ["gene", "Benomyl"]].groupby("gene").max().value_counts()
# Benomyl
# 0          2030 <- the additional genes come from those delimeted by , and "intergenic" group
# 1           128 <- sanity check. 128/1213 = 10.55%

sub = ypdbeno500_opt_snp.loc[ypdbeno500_opt_snp["YPDBENOMYL500"] != 0,:]
sub.groupby("Gene")["Benomyl"].max().value_counts()
# Benomyl
# 0    1987 <- includes "intergenic" group
# 1     128


# ------------------------------------------------------------------------------
# 1. Do SGD benchmark genes have greater SHAP values than non-benchmark genes?
# ------------------------------------------------------------------------------
# salloc -N 1 -c 20 --mem=5GB --time=8:00:00
# conda activate /mnt/home/seguraab/miniconda3/envs/shap
# /mnt/home/seguraab/miniconda3/envs/shap/bin/python

"""_s288c_SHAP_comparisons_sgd_vs_imp_non_bench/
	- Intergenic SNPs were excluded. I did not expand the snp_map (snps with
	multiple gene assignments were excluded).
	- Optimized YPD Benomyl 500 ug/ml model was used to obtain the set of
	important non-benchmark features and a subset of the sgd benchmarks.
	- The SNP optimized model only has about 1/3 of the 376 benomyl sgd benchmark
	genes.
	- The PAV and CNV optimized models had 0 sgd benchmarks.
"""

import os, gc
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from scipy.stats import ks_2samp, mannwhitneyu

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Scripts/Data_Vis")

# SNP and ORF to gene map files
snp_map = pd.read_excel("S8_File.xlsx", engine="openpyxl")
orf_map = pd.read_excel("S9_File.xlsx", engine="openpyxl")

# SGD Benomyl benchmark genes -----------------------
sgd_benchmark = set(snp_map.loc[snp_map.Benomyl==1, "gene"]).union(
	set(orf_map.loc[orf_map.Benomyl==1, "gene"]))
len(sgd_benchmark) # 376 genes mapped to SNPs and/or ORFs
sgd_non_benchmark = set(snp_map.loc[snp_map.Benomyl==0, "gene"]).union(
	set(orf_map.loc[orf_map.Benomyl==0, "gene"]))
# sanity check
sgd_benchmark.isdisjoint(sgd_non_benchmark) # True (no shared elements)
# ---------------------------------------------------

# S288C Costanzo benchmark genes --------------------
# DRF: differential relative fitness values
costanzo_drf = pd.read_excel("../../Data/Costanzo_2021/2021_Costanzo_Data File S1_Conditions_Strains_Fitness.xlsx",
	engine="openpyxl", sheet_name="Diff. Mutant fitness_Conditions")
assert costanzo_drf["Strain ID"].nunique() == len(costanzo_drf) # unique mutant identifiers
costanzo_ben = costanzo_drf.set_index("Strain ID")[["Systematic Name", "Benomyl"]]
costanzo_ben.describe().to_dict()
# {'count': 4383.0, 'mean': -0.0430362051106548, 'std': 0.08260868299038501, 'min': -0.5325, '25%': -0.0815, '50%': -0.0385, '75%': 0.0, 'max': 0.3385}

# check if there are overlapping genes for mutants whose DRF < 0 and DRF >= 0
mask = costanzo_ben[costanzo_ben["Benomyl"] < 0]["Systematic Name"].isin(
	costanzo_ben[costanzo_ben["Benomyl"] >= 0]["Systematic Name"])
mask.sum() # 92 genes, calculate a median Benomyl DRF
costanzo_ben_median = costanzo_ben.groupby("Systematic Name")["Benomyl"].median()
costanzo_ben_median.dropna().shape # (4165,)

# define benchmark and non-benchmark sets
benchmark = costanzo_ben_median[costanzo_ben_median < 0] # mutants with negative differential relative fitness (DRF)
len(benchmark) # 3101 single mutants
non_benchmark = costanzo_ben_median[costanzo_ben_median >= 0] # mutants with non-negative DRF values
len(non_benchmark) # 1064 single mutants
assert set(benchmark.index).isdisjoint(set(non_benchmark.index)) # no overlapping genes between benchmark and non-benchmark sets
# ---------------------------------------------------

# SNP, PAV, and CNV SHAP files ----------------------
data_dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/"
snp_shap_files = glob(os.path.join(data_dir, "SNP/fs/SHAP_values_sorted_Y*.txt"))
pav_shap_files = glob(os.path.join(data_dir, "PAV/fs/SHAP_values_sorted_Y*.txt"))
cnv_shap_files = glob(os.path.join(data_dir, "CNV/fs/SHAP_values_sorted_Y*.txt"))
# snp_shap_files = glob(os.path.join(data_dir, "SNP/baseline/SHAP_values_sorted_Y*.txt"))
# pav_shap_files = glob(os.path.join(data_dir, "PAV/baseline/SHAP_values_sorted_Y*.txt"))
# cnv_shap_files = glob(os.path.join(data_dir, "CNV/baseline/SHAP_values_sorted_Y*.txt"))

snp_shap_file = [f for f in snp_shap_files if "YPDBENOMYL500" in f]
pav_shap_file = [f for f in pav_shap_files if "YPDBENOMYL500" in f]
cnv_shap_file = [f for f in cnv_shap_files if "YPDBENOMYL500" in f]

# SNP identifier mapping file
snp_id_map = pd.read_csv(
	"../../Data/Peter_2018/0_raw_data/mapping_of_geno.csv_feature_names_to_actual_feature_names_09102025.csv",
	index_col=0)

def fix_snp_ids(snp_shap_file):
	"""Fix the SNP feature names in the SHAP files to match the actual SNP
	identifiers"""
	snp_shap = dt.fread(snp_shap_file).to_pandas()
	snp_shap.set_index(snp_shap.columns[0], inplace=True)
	snp_shap.columns = snp_shap.columns.map(snp_id_map["actual_feature"])
	return snp_shap


def fix_orf_ids(orf_shap_file):
	"""Fix the ORF feature names in the SHAP files to match the actual ORF
	identifiers"""
	orf_shap = pd.read_csv(orf_shap_file, sep="\t", index_col=0)
	orf_shap.columns = orf_shap.columns.str.replace("^X", "", regex=True).str.replace("\.", "-", regex=True)
	return orf_shap

# read in SHAP files for the YPDBENOMYL500 environment
snp_shap = fix_snp_ids(snp_shap_file[0]) # 118,382 features (complete model); 6k features (optimized model)
pav_shap = fix_orf_ids(pav_shap_file[0]) # 7708 features (complete); 128 (optimized)
cnv_shap = fix_orf_ids(cnv_shap_file[0]) # 7708 features (complete; 128 (optimized)
# ---------------------------------------------------

# determine which SNPs and ORFs map to benchmark and non_benchmark sets
snp_shap_genes = pd.DataFrame(snp_shap.columns, columns=["snp"]).merge(
	snp_map[["snp", "gene"]], left_on="snp", right_on="snp", how="left")
pav_shap_genes = pd.DataFrame(pav_shap.columns, columns=["orf"]).merge(
	orf_map[["orf", "gene"]], left_on="orf", right_on="orf", how="left")
cnv_shap_genes = pd.DataFrame(cnv_shap.columns, columns=["orf"]).merge(
	orf_map[["orf", "gene"]], left_on="orf", right_on="orf", how="left")

benchmark_dict = benchmark.to_dict()
non_benchmark_dict = non_benchmark.to_dict()
def get_drf_value(row, benchmark=benchmark_dict):
	if row["gene"] in benchmark.keys():
		return benchmark[row["gene"]]
	else:
		return None

snp_shap_genes.insert(2, "Benomyl_DRF", snp_shap_genes.apply(lambda x: get_drf_value(x), axis=1))
pav_shap_genes.insert(2, "Benomyl_DRF", pav_shap_genes.apply(lambda x: get_drf_value(x), axis=1))
cnv_shap_genes.insert(2, "Benomyl_DRF", cnv_shap_genes.apply(lambda x: get_drf_value(x), axis=1))
snp_shap_genes.insert(3, "Non_benchmark", snp_shap_genes.apply(lambda x: get_drf_value(x, non_benchmark_dict), axis=1))
pav_shap_genes.insert(3, "Non_benchmark", pav_shap_genes.apply(lambda x: get_drf_value(x, non_benchmark_dict), axis=1))
cnv_shap_genes.insert(3, "Non_benchmark", cnv_shap_genes.apply(lambda x: get_drf_value(x, non_benchmark_dict), axis=1))
snp_shap_genes.insert(4, "SGD_benchmark", snp_shap_genes.apply(lambda x: 1 if x["gene"] in sgd_benchmark else 0, axis=1))
pav_shap_genes.insert(4, "SGD_benchmark", pav_shap_genes.apply(lambda x: 1 if x["gene"] in sgd_benchmark else 0, axis=1))
cnv_shap_genes.insert(4, "SGD_benchmark", cnv_shap_genes.apply(lambda x: 1 if x["gene"] in sgd_benchmark else 0, axis=1))

# how many SGD benchmark genes are in optimized SNP YPD Benomyl 500 ug/ml model?
snp_shap_genes.groupby("gene")["SGD_benchmark"].max().value_counts()
# SGD_benchmark
# 0    2002
# 1     128
snp_map.loc[snp_map.snp.isin(snp_shap_genes["snp"]), ["gene", "Benomyl"]].groupby("gene").max().value_counts()
# SGD_benchmark
# 0    2002
# 1     128
pav_shap_genes["SGD_benchmark"].value_counts()
# SGD_benchmark
# False    127
# True       1
cnv_shap_genes["SGD_benchmark"].value_counts()
# SGD_benchmark
# False    127
# True       1

print("Number of SNP features with negative DRF values:", snp_shap_genes.isna().value_counts())
'''s288c_SHAP_comparisons_sgd_vs_imp_non_bench/
snp    gene   Benomyl_DRF  Non_benchmark
False  False  True         True              3532    3532+654=4186
              False        True              1814
              True         False              654
'''
snp_shap_genes['snp'].duplicated().sum() # 62 snps (I exploded snp_map); 1234 snps (I exploded snp_map); 0 snps (I did not explode snp_map)

print("Number of PAV features with negative DRF values:", pav_shap_genes.isna().value_counts())
'''s288c_SHAP_comparisons_sgd_vs_imp_non_bench/
orf    gene   Benomyl_DRF  Non_benchmark
False  True   True         True             88
       False  False        True             19   benchmark
              True         True             18
                           False             3   non-benchmark
'''

print("Number of CNV features with negative DRF values:", cnv_shap_genes.isna().value_counts())
'''s288c_SHAP_comparisons_sgd_vs_imp_non_bench/
orf    gene   Benomyl_DRF  Non_benchmark
False  True   True         True               76
       False  True         True               34
              False        True               16  benchmark
              True         False               2
'''
print("Number of unique genes represented in SNP features:", snp_shap_genes.gene.nunique())
# 2130
print("Number of unique genes with negative DRF values:", snp_shap_genes.groupby("gene")["Benomyl_DRF"].median().isna().value_counts())
# Benomyl_DRF
# True     1071
# False    1059
print("Number of unique genes with non-negative DRF values:", snp_shap_genes.groupby("gene")["Non_benchmark"].median().isna().value_counts())
# Non_benchmark
# True     1752
# False     378
print("Number of unique genes represented in PAV features:", pav_shap_genes.gene.nunique())
# 40
print("Number of unique genes with negative DRF values:", pav_shap_genes.groupby("gene")["Benomyl_DRF"].median().isna().value_counts())
# Benomyl_DRF
# True     21
# False    19  <- genes that have negative DRF values
print("Number of unique genes with non-negative DRF values:", pav_shap_genes.groupby("gene")["Non_benchmark"].median().isna().value_counts())
# Non_benchmark
# True     37
# False     3  <- genes that have non-negative DRF values
print("Number of unique genes represented in CNV features:", cnv_shap_genes.gene.nunique())
# 52
print("Number of unique genes with negative DRF values:", cnv_shap_genes.groupby("gene")["Benomyl_DRF"].median().isna().value_counts())
# Benomyl_DRF
# True     36
# False    16
print("Number of unique genes with non-negative DRF values:", cnv_shap_genes.groupby("gene")["Non_benchmark"].median().isna().value_counts())
# Non_benchmark
# True     50
# False     2  <- genes that have non-negative DRF values

# functions to compare local SHAP value distributions of benchmark genes to that of the important non-benchmark genes
def plot_shap_violin(data1, data2, title, save_name, strain="W303"):
	"""Plot a violin plot comparing the distributions of two sets of SHAP values
	Args:
		data1 (np.ndarray) : SHAP values of genes that map to W303 benchmark genes
		data2 (np.ndarray): SHAP values of genes that do NOT map to W303 benchmark genes
	Note: data1 and data2 may have different lengths.
	"""
	df1 = pd.DataFrame({"Gene_type": f"{strain}_genes", "SHAP_value": data1})
	df2 = pd.DataFrame({"Gene_type": f"Non_{strain}_genes", "SHAP_value": data2})
	df_melted = pd.concat([df1, df2], axis=0)
	# plot distribution of values
	fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(6, 6))
	sns.violinplot(data=df_melted, x="Gene_type", y="SHAP_value", inner="quart", ax=ax[0, 0])
	sns.violinplot(data=df_melted, x="Gene_type", y="SHAP_value", inner="quart", ax=ax[0, 1])
	ax[0, 1].set_yscale("symlog", linthresh=1e-4)
	ax[0, 1].set_ylabel("SHAP value (symlog scale)", fontsize=7)
	fig.suptitle(title, fontsize=8)
	ax[0, 0].set_xlabel("")
	ax[0, 1].set_xlabel("")
	ax[0, 0].set_ylabel("SHAP value", fontsize=7)
	ax[0, 0].tick_params(axis="both", labelsize=6)
	ax[0, 1].tick_params(axis="both", labelsize=6)
	# plot CDF of values
	sns.ecdfplot(data=df1, x="SHAP_value", ax=ax[1, 0], label=f"{strain}_genes")
	sns.ecdfplot(data=df2, x="SHAP_value", ax=ax[1, 0], label=f"Non_{strain}_genes")
	ax[1, 0].set_xlabel("SHAP value", fontsize=7)
	ax[1, 0].set_ylabel("Cumulative density", fontsize=7)
	ax[1, 0].tick_params(axis="both", labelsize=6)
	ax[1, 0].legend(fontsize=6)
	plt.tight_layout()
	plt.savefig(save_name)
	plt.close("all")
	gc.collect()
	gc.collect()


def ks_2samp_shap(data1, data2, strain="w303", env="YPDBENOMYL500", vtype="SNP", tag=""):
	# Separate SHAP values by type or use absolute values
	if len(data1.shape) == 1: # global SHAP values
		raw_strain = data1.copy(deep=True)
		raw_non_strain = data2.copy(deep=True)
		abs_strain = data1.abs()
		abs_non_strain = data2.abs()
		pos_strain = data1.mask(data1 <= 0, None).dropna()
		pos_non_strain = data2.mask(data2 <= 0, None).dropna()
		neg_strain = data1.mask(data1 >= 0, None).dropna()
		neg_non_strain = data2.mask(data2 >= 0, None).dropna()
	else: # local SHAP values
		raw_strain = data1["shap_value"].copy(deep=True)
		raw_non_strain = data2["shap_value"].copy(deep=True)
		abs_strain = raw_strain.abs()
		abs_non_strain = raw_non_strain.abs()
		pos_strain = raw_strain.mask(raw_strain <= 0, None).dropna()
		pos_non_strain = raw_non_strain.mask(raw_non_strain <= 0, None).dropna()
		neg_strain = raw_strain.mask(raw_strain >= 0, None).dropna()
		neg_non_strain = raw_non_strain.mask(raw_non_strain >= 0, None).dropna()
	
	# store results
	out_ks = {}
	out_stats = {}
	if len(raw_strain) > 0:
		# raw SHAP values
		out_ks[("ks_test", "raw")] = list(ks_2samp(
			data1=raw_strain, data2=raw_non_strain, alternative="greater", method="asymp"))
		# plot
		plot_shap_violin(data1=raw_strain, data2=raw_non_strain, strain=strain.upper(),
			title=f"{env} {tag} {vtype} SHAP values", save_name=f"_{strain}_{env}_{tag}_{vtype}_raw_SHAP_violin.pdf")
		# summary statistics
		out_stats[("summary_stats", f"raw_{strain}")] = raw_strain.describe()
		out_stats[("summary_stats", f"raw_non_{strain}")] = raw_non_strain.describe()
		del raw_strain, raw_non_strain
	if len(abs_strain) > 0:
		# absolute SHAP values
		out_ks[("ks_test", "absolute")] = list(ks_2samp(
			data1=abs_strain, data2=abs_non_strain, alternative="greater", method="asymp"))
		# plot
		plot_shap_violin(data1=abs_strain, data2=abs_non_strain, strain=strain.upper(),
			title=f"{env} {tag} {vtype} absolute SHAP values", save_name=f"_{strain}_{env}_{tag}_{vtype}_abs_SHAP_violin.pdf")
		out_stats[("summary_stats", f"abs_{strain}")] = abs_strain.describe()
		out_stats[("summary_stats", f"abs_non_{strain}")] = abs_non_strain.describe()
		del abs_strain, abs_non_strain
	if len(pos_strain) > 0:
		# positive SHAP values only
		out_ks[("ks_test", "positive")] = list(ks_2samp(
			data1=pos_strain, data2=pos_non_strain, alternative="greater", method="asymp"))
		plot_shap_violin(data1=pos_strain, data2=pos_non_strain, strain=strain.upper(),
			title=f"{env} {tag} {vtype} Positive SHAP values", save_name=f"_{strain}_{env}_{tag}_{vtype}_pos_SHAP_violin.pdf")
		out_stats[("summary_stats", f"pos_{strain}")] = pos_strain.describe()
		out_stats[("summary_stats", f"pos_non_{strain}")] = pos_non_strain.describe()
		del pos_strain, pos_non_strain
	if len(neg_strain) > 0:
		# negative SHAP values only
		out_ks[("ks_test", "negative")] = list(ks_2samp(
			data1=neg_strain, data2=neg_non_strain, alternative="greater", method="asymp"))
		plot_shap_violin(data1=neg_strain, data2=neg_non_strain, strain=strain.upper(),
			title=f"{env} {tag} {vtype} Negative SHAP values", save_name=f"_{strain}_{env}_{tag}_{vtype}_neg_SHAP_violin.pdf")
		out_stats[("summary_stats", f"neg_{strain}")] = neg_strain.describe()
		out_stats[("summary_stats", f"neg_non_{strain}")] = neg_non_strain.describe()
		del neg_strain, neg_non_strain
	gc.collect()
	gc.collect()
	return out_ks, out_stats


def mwu_shap(data1, data2):
	"""Mann-Whitney U test.
	Null Ho: both distributions come from the same underlying distribution
	Alt ha: values in data1 tend to be greater than values in data2
	"""
	# Separate SHAP values by type or use absolute values
	if len(data1.shape) == 1: # global SHAP values
		raw_strain = data1.copy(deep=True)
		raw_non_strain = data2.copy(deep=True)
		abs_strain = data1.abs()
		abs_non_strain = data2.abs()
		pos_strain = data1.mask(data1 <= 0, None).dropna()
		pos_non_strain = data2.mask(data2 <= 0, None).dropna()
		neg_strain = data1.mask(data1 >= 0, None).dropna()
		neg_non_strain = data2.mask(data2 >= 0, None).dropna()
	else: # local SHAP values
		raw_strain = data1["shap_value"].copy(deep=True)
		raw_non_strain = data2["shap_value"].copy(deep=True)
		abs_strain = raw_strain.abs()
		abs_non_strain = raw_non_strain.abs()
		pos_strain = raw_strain.mask(raw_strain <= 0, None).dropna()
		pos_non_strain = raw_non_strain.mask(raw_non_strain <= 0, None).dropna()
		neg_strain = raw_strain.mask(raw_strain >= 0, None).dropna()
		neg_non_strain = raw_non_strain.mask(raw_non_strain >= 0, None).dropna()
	
	# store results
	out_mwu = {}
	if len(raw_strain) > 0:
		# raw SHAP values
		out_mwu[("mwu_test", "raw")] = list(mannwhitneyu(
			x=raw_strain, y=raw_non_strain, alternative="greater"))
		del raw_strain, raw_non_strain
	if len(abs_strain) > 0:
		# absolute SHAP values
		out_mwu[("mwu_test", "absolute")] = list(mannwhitneyu(
			x=abs_strain, y=abs_non_strain, alternative="greater"))
		del abs_strain, abs_non_strain
	if len(pos_strain) > 0:
		# positive SHAP values only
		out_mwu[("mwu_test", "positive")] = list(mannwhitneyu(
			x=pos_strain, y=pos_non_strain, alternative="greater"))
		del pos_strain, pos_non_strain
	if len(neg_strain) > 0:
		# negative SHAP values only
		out_mwu[("mwu_test", "negative")] = list(mannwhitneyu(
			x=neg_strain, y=neg_non_strain, alternative="greater"))
		del neg_strain, neg_non_strain
	gc.collect()
	gc.collect()
	return out_mwu


# S288C benchmark genes: compare local SHAP distributions ---------------
def get_shap_subsets(shap_df, map_ben_s288c, map_df, vtype="snp", non_bench_map=None):
	if vtype == "snp":
		# SNPs that map to S288C benchmark genes
		shap_s288c = shap_df.loc[:, shap_df.columns.isin(map_ben_s288c["snp"])]
		n_shap_s288c_genes = map_ben_s288c.loc[map_ben_s288c["snp"].isin(shap_s288c.columns), "gene"].nunique()
		n_shap_s288c_features = shap_s288c.shape[1]
		# SNPs that do not map to S288C benchmark genes
		if non_bench_map is None:
			shap_non = shap_df.loc[:, ~shap_df.columns.isin(map_ben_s288c["snp"])] # used to get results in s288c_SHAP_comparisons_non_exp_verified_included/
		else:
			shap_non = shap_df.loc[:, shap_df.columns.isin(non_bench_map["snp"])] # used to get results in s288c_SHAP_comparisons_non_exp_verified_excluded/
		# drop intergenic SNPs
		shap_non_map = map_df.loc[map_df["snp"].isin(shap_non.columns), ["snp", "gene"]].set_index("snp")
		shap_non = shap_non.loc[:, shap_non_map[shap_non_map.gene != "intergenic"].index.unique()]
		n_shap_non_genes = map_df.loc[map_df["snp"].isin(shap_non.columns), "gene"].nunique()
		n_shap_non_features = shap_non.shape[1]
	else:
		# ORFs that map to S288C benchmark genes
		shap_s288c = shap_df.loc[:, shap_df.columns.isin(map_ben_s288c["orf"])]
		n_shap_s288c_features = shap_s288c.shape[1]
		if shap_s288c.shape[1] > 0:
			n_shap_s288c_genes = map_ben_s288c.loc[map_ben_s288c["orf"].isin(shap_s288c.columns), "gene"].nunique()
		else:
			n_shap_s288c_genes = 0
		# ORFs that do not map to S288C benchmark genes
		shap_non = shap_df.loc[:, ~shap_df.columns.isin(map_ben_s288c["orf"])] # used to get results in s288c_SHAP_comparisons_non_exp_verified_included/
		n_shap_non_genes = map_df.loc[map_df["orf"].isin(shap_non.columns), "gene"].nunique()
		n_shap_non_features = shap_non.shape[1]
	print(f"{n_shap_s288c_genes} genes with {vtype} SHAP values that map to S288C benchmark genes")
	print(f"{n_shap_s288c_features} {vtype} features that map to S288C benchmark genes")
	print(f"{n_shap_non_genes} genes with {vtype} SHAP values that do NOT map to S288C benchmark genes")
	print(f"{n_shap_non_features} {vtype} features that do NOT map to S288C benchmark genes")
	return shap_s288c, shap_non, [n_shap_s288c_genes, n_shap_s288c_features, n_shap_non_genes, n_shap_non_features]


def get_global_shap(shap_df, map_ben_s288c, map_df, vtype="snp", is_non=False):
	"""Calculate the global SHAP value for each feature by summing the SHAP values across all isolates"""
	gl_shap = shap_df.abs().median(axis=0)
	gl_shap = gl_shap[gl_shap != 0] # remove un-important features
	print(f"{len(gl_shap)} {vtype} features with non-zero global SHAP values")
	n_shap_features = gl_shap.shape[0]
	if vtype == "snp": # re-count genes and features after removing un-important features
		if not is_non:
			n_shap_genes = map_ben_s288c.loc[map_ben_s288c["snp"].isin(gl_shap.index), "gene"].nunique()
		else:
			n_shap_genes = map_df.loc[map_df["snp"].isin(gl_shap.index), "gene"].nunique()
	else:
		if not is_non:
			n_shap_genes = map_ben_s288c.loc[map_ben_s288c["orf"].isin(gl_shap.index), "gene"].nunique()
		else:
			n_shap_genes = map_df.loc[map_df["orf"].isin(gl_shap.index), "gene"].nunique()
	return gl_shap, [n_shap_genes, n_shap_features]


def reshape_local_shap(s288c_shap_df, non_shap_df, map_df, strain="s288c", map_type="orf"):
	"""Re-arrange the local SHAP values into long format for plotting & performing ks-tests"""
	s288c_shap_long = s288c_shap_df.reset_index().melt(id_vars="ID", var_name="feature", value_name="shap_value")
	non_shap_long = non_shap_df.reset_index().melt(id_vars="ID", var_name="feature", value_name="shap_value")
	# drop un-important local SHAP values
	s288c_shap_long = s288c_shap_long[s288c_shap_long["shap_value"] != 0]
	non_shap_long = non_shap_long[non_shap_long["shap_value"] != 0]
	n_shap_s288c_genes = map_df.loc[map_df[map_type].isin(s288c_shap_long["feature"].unique()), "gene"].nunique()
	n_shap_non_genes = map_df.loc[map_df[map_type].isin(non_shap_long["feature"].unique()), "gene"].nunique()
	return s288c_shap_long, non_shap_long, [n_shap_s288c_genes, s288c_shap_long["feature"].nunique(), n_shap_non_genes, non_shap_long["feature"].nunique()]


# Compare the SHAP value distributions of SGD benchmarks to that of non-benchmarks
# define the snp and orf benchmark sets
snp_map_ben_s288c = snp_map.loc[snp_map["gene"].isin(sgd_benchmark), :]
orf_map_ben_s288c = orf_map.loc[orf_map["gene"].isin(sgd_benchmark), :]

# subset the SHAP values of the SGD benchmark genes and non-benchmark genes
snp_shap_s288c, snp_shap_non, snp_shap_counts = get_shap_subsets(snp_shap, snp_map_ben_s288c, snp_map, vtype="snp")
pav_shap_s288c, pav_shap_non, pav_shap_counts = get_shap_subsets(pav_shap, orf_map_ben_s288c, orf_map, vtype="orf") # not enough benchmarks
cnv_shap_s288c, cnv_shap_non, cnv_shap_counts = get_shap_subsets(cnv_shap, orf_map_ben_s288c, orf_map, vtype="orf") # pav and cnv only have 1 benchmark

# calculate the global shap  value (median absolute shap across instances per feature)
snp_gl_shap_s288c, snp_gl_shap_s288c_counts = get_global_shap(snp_shap_s288c, snp_map_ben_s288c, snp_map, vtype="snp", is_non=False)
snp_gl_shap_non, snp_gl_shap_non_counts = get_global_shap(snp_shap_non, snp_map_ben_s288c, snp_map, vtype="snp", is_non=True)

# re-arrange the local SHAP values into long format for plotting & performing ks-tests
snp_lc_shap_s288c, snp_lc_shap_non, snp_lc_shap_counts = reshape_local_shap(snp_shap_s288c, snp_shap_non, snp_map, strain="s288c", map_type="snp")

# perform the ks test and plot distributions
ks_snp_gl, stats_snp_gl = ks_2samp_shap(snp_gl_shap_s288c, snp_gl_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="SNP", tag=f"sgd_vs_imp_non_bench_global")
ks_snp_lc, stats_snp_lc = ks_2samp_shap(snp_lc_shap_s288c, snp_lc_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="SNP", tag=f"sgd_vs_imp_non_bench_local")

# perform the mann-whitney u test
mwu_snp_gl = mwu_shap(snp_gl_shap_s288c, snp_gl_shap_non)
mwu_snp_lc = mwu_shap(snp_lc_shap_s288c, snp_lc_shap_non)

# store results
ks_res = pd.concat([
	pd.DataFrame.from_dict(ks_snp_gl, orient="index", columns=["statistic", "pvalue"]),
	pd.DataFrame.from_dict(ks_snp_lc, orient="index", columns=["statistic", "pvalue"])], axis=0).reset_index()
ks_res.insert(0, ("Variant_Type", "SHAP_Type", "Comparison"), [("SNP", "global", "sgd_vs_imp_non_bench")]*3 + [("SNP", "local", "sgd_vs_imp_non_bench")]*4)
ks_res[["Variant_Type", "SHAP_Type", "Comparison"]] = pd.DataFrame(ks_res[("Variant_Type", "SHAP_Type", "Comparison")].tolist(), index=ks_res.index)
ks_res[["remove", "SHAP_Subset"]] = pd.DataFrame(ks_res["index"].tolist(), index=ks_res.index)
ks_res = ks_res[["Variant_Type", "SHAP_Type", "Comparison", "SHAP_Subset", "statistic", "pvalue"]]
ks_res.to_csv("_s288c_sgd_vs_imp_non_bench_ks_test_results.csv", index=False)

stats_res = {}
stats_res[("SNP", "global", "sgd_vs_imp_non_bench")] = pd.DataFrame.from_dict(stats_snp_gl, orient="index")
stats_res[("SNP", "local", "sgd_vs_imp_non_bench")] = pd.DataFrame.from_dict(stats_snp_lc, orient="index")
stats_res = pd.concat(stats_res)
stats_res.index.names = ["Variant_Type", "SHAP_Type", "Comparison", "remove", "SHAP_Subset"]
stats_res = stats_res.droplevel("remove")
stats_res.to_csv("_s288c_sgd_vs_imp_non_bench_shap_summary_stats.csv", index=True)

mwu_res = {}
mwu_res[("SNP", "global", "sgd_vs_imp_non_bench")] = pd.DataFrame.from_dict(mwu_snp_gl, orient="index", columns=["statistic", "pvalue"])
mwu_res[("SNP", "local", "sgd_vs_imp_non_bench")] = pd.DataFrame.from_dict(mwu_snp_lc, orient="index", columns=["statistic", "pvalue"])
mwu_res = pd.concat(mwu_res)
mwu_res.index.names = ["Variant_Type", "SHAP_Type", "Comparison", "remove"]
mwu_res.insert(0, "SHAP_Subset", ["raw", "absolute", "positive"]*2 + ["negative"])
mwu_res = mwu_res.droplevel("remove")
mwu_res.to_csv("_s288c_sgd_vs_imp_non_bench_mwu_test_results.csv", index=True)

counts_res = {}
counts_res[("SNP", "before_filtering", "sgd_vs_imp_non_bench")] = snp_shap_counts
counts_res[("SNP", "global", "sgd_vs_imp_non_bench")] = snp_gl_shap_s288c_counts + snp_gl_shap_non_counts
counts_res[("SNP", "local", "sgd_vs_imp_non_bench")] = snp_lc_shap_counts
counts_res = pd.DataFrame(counts_res).T
counts_res.columns = ["Num_S288C_Bench_Genes", "Num_S288C_Bench_Features", "Num_Non_S288C_Bench_Genes", "Num_Non_S288C_Bench_Features"]
counts_res.index.names = ["Variant_Type", "SHAP_Type", "Comparison"]
counts_res.to_csv("_s288c_sgd_vs_imp_non_bench_feature_gene_counts.csv", index=True)


# Compare the SHAP value distributions of Costanzo negative DRF benchmarks to that of non-benchmarks
ks_res = {}
mwu_res = {}
stats_res = {}
counts_res = {}
for q in [1, .0005, .005, .01, .025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9]:
	# Get the SNPs & ORFs that map to S288C benchmark genes
	# benchmark genes are those whose mutants have decreased fitness in benomyl compared to the reference condition
	# use the quantile to define mutants that are sufficiently negatively affected by benomyl
	print()
	if q == 1:
		print("Using all negative Benomyl DRF values to define benchmark genes")
		snp_map_ben_s288c = snp_map.loc[snp_map["gene"].isin(benchmark.index), :] # SNP Costanzo negative DRF benchmark genes
		orf_map_ben_s288c = orf_map.loc[orf_map["gene"].isin(benchmark.index), :] # ORF Costanzo negative DRF benchmark genes
	else:
		drf_thresh = benchmark.quantile(q, interpolation="linear")
		print(f"Benomyl DRF threshold for top {q*100}% of mutants: {drf_thresh}")
		benchmark_q = benchmark[benchmark <= drf_thresh] # the top q (%) of mutants most negatively affected by benomyl
		snp_map_ben_s288c = snp_map.loc[snp_map["gene"].isin(benchmark_q.index), :] # SNP Costanzo negative DRF benchmark genes
		orf_map_ben_s288c = orf_map.loc[orf_map["gene"].isin(benchmark_q.index), :] # ORF Costanzo negative DRF benchmark genes
	
	# subset the SHAP values of S288C benchmark genes and non-benchmark genes
	non_bench_map = snp_map.loc[snp_map["gene"].isin(non_benchmark.index), :]
	snp_shap_s288c, snp_shap_non, snp_shap_counts = get_shap_subsets(snp_shap, snp_map_ben_s288c, snp_map, vtype="snp", non_bench_map=non_bench_map)
	pav_shap_s288c, pav_shap_non, pav_shap_counts = get_shap_subsets(pav_shap, orf_map_ben_s288c, orf_map, vtype="pav")
	cnv_shap_s288c, cnv_shap_non, cnv_shap_counts = get_shap_subsets(cnv_shap, orf_map_ben_s288c, orf_map, vtype="cnv")
	
	# calculate the global SHAP value (median absolute shap across instances per feature); drop un-important features
	snp_gl_shap_s288c, snp_gl_shap_s288c_counts = get_global_shap(snp_shap_s288c, snp_map_ben_s288c, snp_map, vtype="snp", is_non=False)
	snp_gl_shap_non, snp_gl_shap_non_counts = get_global_shap(snp_shap_non, snp_map_ben_s288c, snp_map, vtype="snp", is_non=True)
	pav_gl_shap_s288c, pav_gl_shap_s288c_counts = get_global_shap(pav_shap_s288c, orf_map_ben_s288c, orf_map, vtype="pav", is_non=False)
	pav_gl_shap_non, pav_gl_shap_non_counts = get_global_shap(pav_shap_non, orf_map_ben_s288c, orf_map, vtype="pav", is_non=True)
	cnv_gl_shap_s288c, cnv_gl_shap_s288c_counts = get_global_shap(cnv_shap_s288c, orf_map_ben_s288c, orf_map, vtype="cnv", is_non=False)
	cnv_gl_shap_non, cnv_gl_shap_non_counts = get_global_shap(cnv_shap_non, orf_map_ben_s288c, orf_map, vtype="cnv", is_non=True)
	
	# re-arrange the local SHAP values into long format for plotting & performing ks-tests
	snp_lc_shap_s288c, snp_lc_shap_non, snp_lc_shap_counts = reshape_local_shap(snp_shap_s288c, snp_shap_non, snp_map, strain="s288c", map_type="snp")
	pav_lc_shap_s288c, pav_lc_shap_non, pav_lc_shap_counts = reshape_local_shap(pav_shap_s288c, pav_shap_non, orf_map, strain="s288c", map_type="orf")
	cnv_lc_shap_s288c, cnv_lc_shap_non, cnv_lc_shap_counts = reshape_local_shap(cnv_shap_s288c, cnv_shap_non, orf_map, strain="s288c", map_type="orf")
	
	# perform the ks test and plot distributions
	ks_snp_gl, stats_snp_gl = ks_2samp_shap(snp_gl_shap_s288c, snp_gl_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="SNP", tag=f"top_{q*100}pct_global")
	ks_pav_gl, stats_pav_gl = ks_2samp_shap(pav_gl_shap_s288c, pav_gl_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="PAV", tag=f"top_{q*100}pct_global")
	ks_cnv_gl, stats_cnv_gl = ks_2samp_shap(cnv_gl_shap_s288c, cnv_gl_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="CNV", tag=f"top_{q*100}pct_global")
	ks_snp_lc, stats_snp_lc = ks_2samp_shap(snp_lc_shap_s288c, snp_lc_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="SNP", tag=f"top_{q*100}pct_local")
	ks_pav_lc, stats_pav_lc = ks_2samp_shap(pav_lc_shap_s288c, pav_lc_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="PAV", tag=f"top_{q*100}pct_local")
	ks_cnv_lc, stats_cnv_lc = ks_2samp_shap(cnv_lc_shap_s288c, cnv_lc_shap_non, strain="s288c", env="YPDBENOMYL500", vtype="CNV", tag=f"top_{q*100}pct_local")
	
	# perform mann-whitney u test
	mwu_snp_gl = mwu_shap(snp_gl_shap_s288c, snp_gl_shap_non)
	mwu_pav_gl = mwu_shap(pav_gl_shap_s288c, pav_gl_shap_non)
	mwu_cnv_gl = mwu_shap(cnv_gl_shap_s288c, cnv_gl_shap_non)
	mwu_snp_lc = mwu_shap(snp_lc_shap_s288c, snp_lc_shap_non)
	mwu_pav_lc = mwu_shap(pav_lc_shap_s288c, pav_lc_shap_non)
	mwu_cnv_lc = mwu_shap(cnv_lc_shap_s288c, cnv_lc_shap_non)
	
	# store results
	ks_res[("SNP", "global", q)] = pd.DataFrame.from_dict(ks_snp_gl, orient="index", columns=["statistic", "pvalue"])
	stats_res[("SNP", "global", q)] = pd.DataFrame.from_dict(stats_snp_gl, orient="index")
	ks_res[("PAV", "global", q)] = pd.DataFrame.from_dict(ks_pav_gl, orient="index", columns=["statistic", "pvalue"])
	stats_res[("PAV", "global", q)] = pd.DataFrame.from_dict(stats_pav_gl, orient="index")
	ks_res[("CNV", "global", q)] = pd.DataFrame.from_dict(ks_cnv_gl, orient="index", columns=["statistic", "pvalue"])
	stats_res[("CNV", "global", q)] = pd.DataFrame.from_dict(stats_cnv_gl, orient="index")
	
	ks_res[("SNP", "local", q)] = pd.DataFrame.from_dict(ks_snp_lc, orient="index", columns=["statistic", "pvalue"])
	stats_res[("SNP", "local", q)] = pd.DataFrame.from_dict(stats_snp_lc, orient="index")
	ks_res[("PAV", "local", q)] = pd.DataFrame.from_dict(ks_pav_lc, orient="index", columns=["statistic", "pvalue"])
	stats_res[("PAV", "local", q)] = pd.DataFrame.from_dict(stats_pav_lc, orient="index")
	ks_res[("CNV", "local", q)] = pd.DataFrame.from_dict(ks_cnv_lc, orient="index", columns=["statistic", "pvalue"])
	stats_res[("CNV", "local", q)] = pd.DataFrame.from_dict(stats_cnv_lc, orient="index")
	
	mwu_res[("SNP", "global", q)] = pd.DataFrame.from_dict(mwu_snp_gl, orient="index", columns=["statistic", "pvalue"])
	mwu_res[("PAV", "global", q)] = pd.DataFrame.from_dict(mwu_pav_gl, orient="index", columns=["statistic", "pvalue"])
	mwu_res[("CNV", "global", q)] = pd.DataFrame.from_dict(mwu_cnv_gl, orient="index", columns=["statistic", "pvalue"])
	mwu_res[("SNP", "local", q)] = pd.DataFrame.from_dict(mwu_snp_lc, orient="index", columns=["statistic", "pvalue"])
	mwu_res[("PAV", "local", q)] = pd.DataFrame.from_dict(mwu_pav_lc, orient="index", columns=["statistic", "pvalue"])
	mwu_res[("CNV", "local", q)] = pd.DataFrame.from_dict(mwu_cnv_lc, orient="index", columns=["statistic", "pvalue"])
	
	counts_res[("SNP", "before_filtering", q)] = snp_shap_counts
	counts_res[("PAV", "before_filtering", q)] = pav_shap_counts
	counts_res[("CNV", "before_filtering", q)] = cnv_shap_counts
	counts_res[("SNP", "global", q)] = snp_gl_shap_s288c_counts + snp_gl_shap_non_counts
	counts_res[("PAV", "global", q)] = pav_gl_shap_s288c_counts + pav_gl_shap_non_counts
	counts_res[("CNV", "global", q)] = cnv_gl_shap_s288c_counts + cnv_gl_shap_non_counts
	counts_res[("SNP", "local", q)] = snp_lc_shap_counts
	counts_res[("PAV", "local", q)] = pav_lc_shap_counts
	counts_res[("CNV", "local", q)] = cnv_lc_shap_counts
	del snp_shap_s288c, snp_shap_non, snp_shap_counts, snp_gl_shap_s288c, snp_gl_shap_non, snp_gl_shap_s288c_counts, snp_gl_shap_non_counts, snp_lc_shap_s288c, snp_lc_shap_non, snp_lc_shap_counts, mwu_snp_gl, mwu_snp_lc #, ks_snp_lc, stats_snp_lc, ks_snp_gl, stats_snp_gl
	del pav_shap_s288c, pav_shap_non, pav_shap_counts, pav_gl_shap_s288c, pav_gl_shap_non, pav_gl_shap_s288c_counts, pav_gl_shap_non_counts, pav_lc_shap_s288c, pav_lc_shap_non, pav_lc_shap_counts, mwu_pav_gl, mwu_pav_lc #, ks_pav_lc, stats_pav_lc, ks_pav_gl, stats_pav_gl
	del cnv_shap_s288c, cnv_shap_non, cnv_shap_counts, cnv_gl_shap_s288c, cnv_gl_shap_non, cnv_gl_shap_s288c_counts, cnv_gl_shap_non_counts, cnv_lc_shap_s288c, cnv_lc_shap_non, cnv_lc_shap_counts, mwu_cnv_gl, mwu_cnv_lc #, ks_cnv_lc, stats_cnv_lc, ks_cnv_gl, stats_cnv_gl

# Save the KS-test results
df_list = []
num_empty = 0
for key in ks_res.keys():
	if key == ('CNV', 'local', 0.95):
		continue
	if not ks_res[key].empty:
		new_idx = [list(key) + list(x) for x in ks_res[key].reset_index()["index"]]
		ks_res[key].index = pd.MultiIndex.from_tuples(new_idx)
		df_list.append(ks_res[key])
	else:
		num_empty += 1

ks_res_df = pd.concat(df_list)
ks_res_df.index.names = ["Variant_Type", "SHAP_Type", "Benomyl_DRF_Quantile", "remove", "SHAP_Subset"]
ks_res_df.reset_index(inplace=True)
ks_res_df.drop(columns="remove", inplace=True)
ks_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_ks_test_results.csv", index=False)

assert len(df_list) + num_empty == len(ks_res.keys()) # pass!

# Save the Mann-Whitney U test results
df_list = []
num_empty = 0
for key in mwu_res.keys():
	if not mwu_res[key].empty:
		new_idx = [list(key) + list(x) for x in mwu_res[key].reset_index()["index"]]
		mwu_res[key].index = pd.MultiIndex.from_tuples(new_idx)
		df_list.append(mwu_res[key])
	else:
		num_empty += 1

mwu_res_df = pd.concat(df_list)
mwu_res_df.index.names = ["Variant_Type", "SHAP_Type", "Benomyl_DRF_Quantile", "remove", "SHAP_Subset"]
mwu_res_df.reset_index(inplace=True)
mwu_res_df.drop(columns="remove", inplace=True)
mwu_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_mwu_test_results.csv", index=False)

# apply multiple-testing correction
from statsmodels.stats.multitest import multipletests
mwu_res_df = mwu_res_df.pivot(index=["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type"], columns="SHAP_Subset", values="pvalue")
mwu_res_df.reset_index(inplace=True)

mwu_res_df[["raw_q", "absolute_q", "positive_q", "negative_q"]] = np.nan
for shap_type in mwu_res_df["SHAP_Type"].unique():
	for variant_type in mwu_res_df["Variant_Type"].unique():
		mask = (mwu_res_df["SHAP_Type"] == shap_type) & (mwu_res_df["Variant_Type"] == variant_type)
		sub = mwu_res_df.loc[mask, ["raw_q", "absolute_q", "positive_q", "negative_q"]]
		sub = mwu_res_df.loc[mask, ["raw", "absolute", "positive", "negative"]].apply(lambda x:multipletests(x, alpha=0.05, method='fdr_bh')[1], axis=0)
		mwu_res_df.loc[mask, ["raw_q", "absolute_q", "positive_q", "negative_q"]] = sub.values

mwu_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_mwu_test_results_with_qvals.csv", index=False)

# Save the distribution descriptive statistics
df_list = {k: v for k, v in stats_res.items() if not v.empty}
stats_res_df = pd.concat(df_list)
stats_res_df.index.names = ["Variant_Type", "SHAP_Type", "Benomyl_DRF_Quantile", "remove", "SHAP_Subset"]
stats_res_df.reset_index(inplace=True)
stats_res_df.drop(columns="remove", inplace=True)
stats_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_shap_summary_stats.csv", index=False)

counts_res_df = pd.DataFrame(counts_res).T
counts_res_df.index.names = ["Variant_Type", "SHAP_Type", "Benomyl_DRF_Quantile"]
counts_res_df.columns = ["Num_S288C_Bench_Genes", "Num_S288C_Bench_Features", "Num_Non_S288C_Bench_Genes", "Num_Non_S288C_Bench_Features"]
counts_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_feature_gene_counts.csv", index=True)

# Plot summary figure of ks_res_df and msu_res_df
ks_res_df = pd.read_csv("_s288c_costanzo_benomyl_benchmark_ks_test_results.csv")
ks_res_df = ks_res_df.pivot(index=["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type"], columns="SHAP_Subset", values="pvalue")
ks_res_df.reset_index(inplace=True)

# apply multiple testing correction to the p-values
ks_res_df[["raw_q", "absolute_q", "positive_q", "negative_q"]] = np.nan
for shap_type in ks_res_df["SHAP_Type"].unique():
	for variant_type in ks_res_df["Variant_Type"].unique():
		mask = (ks_res_df["SHAP_Type"] == shap_type) & (ks_res_df["Variant_Type"] == variant_type)
		sub = ks_res_df.loc[mask, ["raw_q", "absolute_q", "positive_q", "negative_q"]]
		sub = ks_res_df.loc[mask, ["raw", "absolute", "positive", "negative"]].apply(lambda x:multipletests(x, alpha=0.05, method='fdr_bh')[1], axis=0)
		ks_res_df.loc[mask, ["raw_q", "absolute_q", "positive_q", "negative_q"]] = sub.values

ks_res_df.to_csv("_s288c_costanzo_benomyl_benchmark_ks_test_results_with_qvals.csv", index=False)

mwu_res_df = pd.read_csv("_s288c_costanzo_benomyl_benchmark_mwu_test_results_with_qvals.csv")

# plot
def plot_ks_results(res_df, ylabel,save_name):
	fig, ax = plt.subplots(figsize=(11, 8), nrows=3, ncols=2, sharey=True, sharex=True)
	res_df[(res_df["SHAP_Type"] == "global") & (res_df["Variant_Type"] == "SNP")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[0, 0], title="SNP global SHAP values", linewidth=.8, ms=4, grid=True)
	res_df[(res_df["SHAP_Type"] == "local") & (res_df["Variant_Type"] == "SNP")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[0, 1], title="SNP local SHAP values", linewidth=.8, ms=4, grid=True)
	res_df[(res_df["SHAP_Type"] == "global") & (res_df["Variant_Type"] == "PAV")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[1, 0], title="PAV global SHAP values", linewidth=.8, ms=4, grid=True)
	res_df[(res_df["SHAP_Type"] == "local") & (res_df["Variant_Type"] == "PAV")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[1, 1], title="PAV local SHAP values", linewidth=.8, ms=4, grid=True)
	res_df[(res_df["SHAP_Type"] == "global") & (res_df["Variant_Type"] == "CNV")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[2, 0], title="CNV global SHAP values", linewidth=.8, ms=4, grid=True)
	res_df[(res_df["SHAP_Type"] == "local") & (res_df["Variant_Type"] == "CNV")].plot(
		x="Benomyl_DRF_Quantile", kind="line", marker="o", ax=ax[2, 1], title="CNV local SHAP values", linewidth=.8, ms=4, grid=True)
	# plot horizontal line at y=0.05
	for i in range(3):
		for j in range(2):
			ax[i, j].axhline(y=0.05, color="red", linestyle="--", linewidth=.8)
			ax[i, j].tick_params(axis="both", labelsize=7)
			ax[i, j].set_ylabel(ylabel)
	plt.savefig(save_name)
	plt.close("all")

plot_ks_results(ks_res_df[["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type", "absolute", "negative", "positive", "raw"]], "P-value", save_name='_s288c_costanzo_benomyl_benchmark_ks_test_results_pvals.pdf')
plot_ks_results(ks_res_df[["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type", "absolute_q", "negative_q", "positive_q", "raw_q"]], "Q-value", save_name='_s288c_costanzo_benomyl_benchmark_ks_test_results_qvals.pdf')
plot_ks_results(mwu_res_df[["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type", "absolute", "negative", "positive", "raw"]], "P-value", save_name='_s288c_costanzo_benomyl_benchmark_mwu_test_results_pvals.pdf')
plot_ks_results(mwu_res_df[["SHAP_Type", "Benomyl_DRF_Quantile", "Variant_Type", "absolute_q", "negative_q", "positive_q", "raw_q"]], "Q-value", save_name='_s288c_costanzo_benomyl_benchmark_mwu_test_results_qvals.pdf')


# ------------------------------------------------------------------------------
# 2. To what extent do the YPD Benomyl 500 ug/ml SHAP values of SGD benchmark
# genes agree with Costanzo et al., 2021 differential relative fitness values?
# ------------------------------------------------------------------------------
from glob import glob

# Load feature-to-gene mapping files
map_snps = pd.read_excel("Scripts/Data_Vis/S8_File.xlsx", engine="openpyxl")
map_orfs = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl")

# SNP identifier mapping file for fixing SNP feature IDs
snp_id_map = pd.read_csv(
	"Data/Peter_2018/0_raw_data/mapping_of_geno.csv_feature_names_to_actual_feature_names_09102025.csv",
	index_col=0)

# Functions to clean SNP & ORF feature IDs for comparing SHAP tables to the feature-to-gene mapping files
def fix_snp_ids(snp_shap_cols, snp_id_map):
	"""Fix the SNP feature names in the SHAP files to match the actual SNP
	identifiers"""
	return snp_shap_cols.map(snp_id_map["actual_feature"])


def fix_orf_ids(orf_shap_cols):
	"""Fix the ORF feature names in the SHAP files to match the actual ORF
	identifiers"""
	return orf_shap_cols.str.replace("^X", "", regex=True).str.replace("\.", "-", regex=True)


# Load in SHAP values
data_dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/"
snp_shap_f = glob(os.path.join(data_dir, "SNP/baseline/SHAP_values_sorted_YPDBENOMYL500*.txt"))
pav_shap_f = glob(os.path.join(data_dir, "PAV/baseline/SHAP_values_sorted_YPDBENOMYL500*.txt"))
cnv_shap_f = glob(os.path.join(data_dir, "CNV/baseline/SHAP_values_sorted_YPDBENOMYL500*.txt"))
snp_shap = dt.fread(snp_shap_f[0]).to_pandas().set_index("ID")
pav_shap = pd.read_csv(pav_shap_f[0], sep="\t", index_col=0)
cnv_shap = pd.read_csv(cnv_shap_f[0], sep="\t", index_col=0)
snp_shap.columns = fix_snp_ids(snp_shap.columns, snp_id_map)
pav_shap.columns = fix_orf_ids(pav_shap.columns)
cnv_shap.columns = fix_orf_ids(cnv_shap.columns)

# Load in Costanzo et al., DRF values
drf = pd.read_excel("Data/Costanzo_2021/2021_Costanzo_Data File S1_Conditions_Strains_Fitness.xlsx",
	engine="openpyxl", sheet_name="Diff. Mutant fitness_Conditions")
assert drf["Strain ID"].nunique() == len(drf) # unique mutant identifiers
drf_beno = drf.set_index("Strain ID")[["Systematic Name", "Benomyl"]]

# Take the absolute value of the DRF and summarize
drf_beno["Benomyl_abs"] = drf_beno["Benomyl"].abs()
drf_beno_summary = drf_beno.groupby("Systematic Name").describe()
drf_beno_summary = drf_beno_summary.dropna(thresh=7)
from scipy.stats import pearsonr
pearsonr(drf_beno_summary[('Benomyl', 'mean')], drf_beno_summary[('Benomyl', '50%')], alternative='two-sided')
# PearsonRResult(statistic=0.9991982949431055, pvalue=0.0)
pearsonr(drf_beno_summary[('Benomyl_abs', 'mean')], drf_beno_summary[('Benomyl_abs', '50%')], alternative='two-sided')
# PearsonRResult(statistic=0.999241692975126, pvalue=0.0)

# How many benchmark genes are recovered in different percentiles? ---
# SGD benomyl benchmark genes
sgd_benchmark = set(map_snps.loc[map_snps.Benomyl==1, "gene"]).union(set(map_orfs.loc[map_orfs.Benomyl==1, "gene"]))
len(sgd_benchmark) # 376 (370 are SNPs, 350 are ORFs)
len(set(drf_beno_summary.index).intersection(sgd_benchmark)) # 201 of the 376 SGD benomyl benchmark genes have DRF values
sgd_benchmarks_common = set(map_snps.loc[map_snps.Benomyl==1, "gene"]).intersection(
	set(map_orfs.loc[map_orfs.Benomyl==1, "gene"])).intersection(
	set(drf_beno_summary.index))
len(sgd_benchmarks_common) # 192 genes common to SNPs, ORFs, and Costanzo DRF


# calculate mean and median of absolute SHAP values
snp_shap_summary = snp_shap.abs().aggregate(["mean", "median", "max"], axis=0).T
pav_shap_summary = pav_shap.abs().aggregate(["mean", "median", "max"], axis=0).T
cnv_shap_summary = cnv_shap.abs().aggregate(["mean", "median", "max"], axis=0).T
# insert the gene information
snp_shap_summary = snp_shap_summary.merge(map_snps[["snp", "gene"]], left_index=True, right_on="snp", how="left").set_index("snp")
pav_shap_summary.insert(3, "gene", pav_shap_summary.index.map(lambda x: map_orfs.set_index("orf")["gene"].get(x, x)))
cnv_shap_summary.insert(3, "gene", cnv_shap_summary.index.map(lambda x: map_orfs.set_index("orf")["gene"].get(x, x)))

# calculate the number of benchmark genes above the SHAP and DRF thresholds
def calc_benchmark_threshold_counts(sgd_benchmarks_common):
	percentiles = [99, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50] # for SHAP & DRF
	shap_counts = []
	drf_counts = []
	drf_abs_counts = []
	for value in ["mean", "median", "max"]:
		# calculate the DRF threshold values
		if (value == "mean") or (value == "median"):
			if value == "median":
				value = "50%"
			drf_thresh_list = np.percentile(
				drf_beno_summary.loc[drf_beno_summary[('Benomyl', value)] != 0, ('Benomyl', value)], percentiles)
			drf_abs_thresh_list = np.percentile(
				drf_beno_summary.loc[drf_beno_summary[("Benomyl_abs", value)] != 0, ("Benomyl_abs", value)], percentiles)
		else:
			drf_thresh_list = np.percentile(
				drf_beno_summary.loc[drf_beno_summary[('Benomyl', "mean")] != 0, ('Benomyl', "mean")], percentiles)
			drf_abs_thresh_list = np.percentile(
				drf_beno_summary.loc[drf_beno_summary[("Benomyl_abs", "mean")] != 0, ("Benomyl_abs", "mean")], percentiles)
		
		# count the number of benchmark genes above the DRF thresholds
		for i in range(len(percentiles)):
			if (value == "mean") or (value == "50%"):
				drf_count = len(set(drf_beno_summary.loc[
					drf_beno_summary[('Benomyl', value)] >= drf_thresh_list[i], :].index).\
					intersection(sgd_benchmarks_common))
				drf_abs_count = len(set(drf_beno_summary.loc[
					drf_beno_summary[("Benomyl_abs", value)] >= drf_abs_thresh_list[i], :].index).\
					intersection(sgd_benchmarks_common))
				drf_counts.append([percentiles[i], value, "DRF", drf_thresh_list[i], drf_count])
				drf_abs_counts.append([percentiles[i], value, "absolute DRF", drf_abs_thresh_list[i], drf_abs_count])
			else:
				pass
		
		if value == "50%":
			value = "median"
		for vtype, shap_summary in zip(["snp", "pav", "cnv"], [snp_shap_summary, pav_shap_summary, cnv_shap_summary]):
			# calculate the SHAP threshold values
			shap_thresh_list = np.percentile(
				shap_summary.loc[shap_summary[value] != 0, value], percentiles)
			# count the number of benchmark genes above the SHAP and DRF thresholds
			for i in range(len(percentiles)):
				shap_count = len(set(shap_summary.loc[
					shap_summary[value] >= shap_thresh_list[i], "gene"].unique()).\
					intersection(sgd_benchmarks_common))
				shap_counts.append([percentiles[i], value, f"{vtype} absolute SHAP", shap_thresh_list[i], shap_count])
	
	out = pd.concat([pd.DataFrame(shap_counts), pd.DataFrame(drf_counts), pd.DataFrame(drf_abs_counts)], axis=0)
	out.columns=["percentile", "transformed value", "value category", "threshold", "benchmark_count"]
	out["value type"] = out["transformed value"] + " " + out["value category"]
	return out


def plot_benchmark_threshold_counts(df, save):
	df.plot.bar(figsize=(8, 4))
	g.set_xlabel("Percentile", fontsize=8)
	g.set_ylabel("Number of SGD benomyl benchmark genes", fontsize=8)
	g.tick_params(axis="both", labelsize=7)
	plt.savefig(save)
	plt.close("all")


out_sgd_common = calc_benchmark_threshold_counts(sgd_benchmarks_common)

# plot the mean DRF or abs DRF counts with the mean absolute SHAP value gene counts
plot_benchmark_threshold_counts(
	df=out_sgd_common[out_sgd_common["transformed value"]=="mean"].pivot(
		index="percentile", columns="value type", values="benchmark_count"),
	save="Scripts/Data_Vis/_shap_and_drf_mean_num_common_benchmark_genes_recovered.pdf")

# plot the median DRF or abs DRF counts with the median absolute SHAP value gene counts
plot_benchmark_threshold_counts(
	df=out_sgd_common[out_sgd_common["transformed value"].isin(["median", "50%"])].pivot(
		index="percentile", columns="value type", values="benchmark_count"),
	save="Scripts/Data_Vis/_shap_and_drf_median_num_common_benchmark_genes_recovered.pdf")

# plot the mean DRF or abs DRF counts with the max absolute SHAP value gene counts
subset = out_sgd_common[(out_sgd_common["transformed value"].isin(["max", "mean"])) & (
			  ~out_sgd_common["value type"].str.contains("mean snp absolute SHAP")) & (
			  ~out_sgd_common["value type"].str.contains("mean pav absolute SHAP")) & (
			  ~out_sgd_common["value type"].str.contains("mean cnv absolute SHAP"))]
plot_benchmark_threshold_counts(
	df=subset.pivot(index="percentile", columns="value type", values="benchmark_count"),
	save="Scripts/Data_Vis/_shap_max_and_drf_mean_num_common_benchmark_genes_recovered.pdf")

# plot the mean DRF or abs DRF counts with the median absolute SHAP value gene counts
subset = out_sgd_common[(out_sgd_common["transformed value"].isin(["median", "mean"])) & (
			  ~out_sgd_common["value type"].str.contains("mean snp absolute SHAP")) & (
			  ~out_sgd_common["value type"].str.contains("mean pav absolute SHAP")) & (
			  ~out_sgd_common["value type"].str.contains("mean cnv absolute SHAP"))]
plot_benchmark_threshold_counts(
	df=subset.pivot(index="percentile", columns="value type", values="benchmark_count"),
	save="Scripts/Data_Vis/_shap_median_and_drf_mean_num_common_benchmark_genes_recovered.pdf")



# Compare the SHAP values of SGD benchmark genes to the DRF values ---
# subset the SGD benchmark gene features from the non-benchmark gene features
snp_shap_beno = snp_shap.loc[:, snp_shap.columns.isin(map_snps.loc[map_snps.Benomyl==1, "snp"])]
snp_shap_nonb = snp_shap.loc[:, snp_shap.columns.isin(map_snps.loc[map_snps.Benomyl==0, "snp"])]
pav_shap_beno = pav_shap.loc[:, pav_shap.columns.isin(map_orfs.loc[map_orfs.Benomyl==1, "orf"])]
pav_shap_nonb = pav_shap.loc[:, pav_shap.columns.isin(map_orfs.loc[map_orfs.Benomyl==0, "orf"])]
cnv_shap_beno = cnv_shap.loc[:, cnv_shap.columns.isin(map_orfs.loc[map_orfs.Benomyl==1, "orf"])]
cnv_shap_nonb = cnv_shap.loc[:, cnv_shap.columns.isin(map_orfs.loc[map_orfs.Benomyl==0, "orf"])]


def shap_vs_drf(df_shap_beno, df_shap_nonb, drf_beno_av, map_df, vtype="snp", ylab="SNPs", drf_type="", outf=""):
	"""Calculate the pearson correlation between SHAP values and DRF values of:
		- SGD benomyl benchmark genes, or
	    - Non-benchmark genes
	   Mean (across features) of the mean absolute SHAP values (across isolates) were estimated.
	   Genes with zero DRF values were dropped.
	   Args:
		- df_shap_beno (pd.DataFrame): SHAP values of benomyl SGD benchmark gene features
	    - df_shap_nonb (pd.DataFrame): SHAP values of non-benchmark gene features
		- drf_beno_av (pd.DataFrame): mean DRF values per gene
	    - map_df (pd.DataFrame): map_snps or map_orfs
		- vtype (str): "snp" or "orf" (use "orf" for pav- and cnv-based SHAP values)
		- ylab (str): SNPs, PAVs, or CNVs, will be used in y-axis label
		- drf_type (str): DRF value type, will be used as x-axis label
		- outf (str): name to save plot as
	"""
	
	# get the benchmark and non-benchmark genes associated with the features
	beno_genes = map_df.loc[map_df[vtype].isin(df_shap_beno.columns), "gene"].unique()
	nonb_genes = map_df.loc[map_df[vtype].isin(df_shap_nonb.columns), "gene"].unique()
	# subset the DRF dataset
	drf_beno_genes = drf_beno_av.loc[drf_beno_av.index.isin(beno_genes)]
	drf_nonb_genes = drf_beno_av.loc[drf_beno_av.index.isin(nonb_genes)]
	# calculate the mean absolute SHAP value (across isolates)
	av_shap_beno = df_shap_beno.abs().mean(axis=0)
	av_shap_nonb = df_shap_nonb.abs().mean(axis=0)
	av_shap_beno.name = av_shap_nonb.name = "mean_abs_shap"
	# estimate the mean of the mean absolute SHAP values
	y_beno = map_df[[vtype, "gene"]].merge(av_shap_beno, left_on=vtype, right_index=True, how="right")
	y_av_beno = y_beno.groupby("gene").describe()
	y_nonb = map_df[[vtype, "gene"]].merge(av_shap_nonb, left_on=vtype, right_index=True, how="right")
	y_av_nonb = y_nonb.groupby("gene").describe()
	if "intergenic" in y_av_nonb.index:
		y_av_nonb.drop(index="intergenic", inplace=True)
	# subset the genes represented in the SHAP datasets that have DRF values
	x_beno = drf_beno_genes.dropna()
	x_nonb = drf_nonb_genes.dropna()
	# estimate the pearson correlation
	r_beno, p_beno = pearsonr(
		x=x_beno, y=y_av_beno.loc[x_beno.index, ("mean_abs_shap", "mean")],
		alternative="two-sided")
	r_nonb, p_nonb = pearsonr(
		x=x_nonb, y=y_av_nonb.loc[x_nonb.index, ("mean_abs_shap", "mean")],
		alternative="two-sided")
	# density plot
	fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6, 3), sharey=True)
	sns.kdeplot(x=x_beno, y=y_av_beno.loc[x_beno.index, ("mean_abs_shap", "mean")],
		levels=50, fill=True, cmap="viridis", legend=True, ax=ax[0])
	sns.kdeplot(x=x_nonb, y=y_av_nonb.loc[x_nonb.index, ("mean_abs_shap", "mean")],
		levels=50, fill=True, cmap="viridis", legend=True, ax=ax[1])
	ax[0].set_ylabel(f"mean (across {ylab}) of mean absolute SHAP", fontsize=8)
	ax[0].set_title(f"SGD benomyl benchmark genes\nr={r_beno:.2f}, p={p_beno:.2e}, n={len(x_beno)}", fontsize=8)
	ax[1].set_title(f"Non-benchmark genes\nr={r_nonb:.2f}, p={p_nonb:.2e}, n={len(x_nonb)}", fontsize=8)
	for i in [0, 1]:
		ax[i].set_xlabel(drf_type, fontsize=8)
		ax[i].tick_params(axis="both", labelsize=7)
	plt.tight_layout()
	plt.savefig(outf)
	plt.close("all")
	gc.collect()
	gc.collect()
	return [len(x_beno), r_beno, p_beno, len(x_nonb), r_nonb, p_nonb]


# Debugging args
# df_shap_beno = pav_shap_beno.copy(deep=True)
# df_shap_nonb = cnv_shap_nonb.copy(deep=True)
# drf_beno_av = drf_beno_summary[("Benomyl_abs", "mean")].copy(deep=True)
# map_df = map_orfs.copy(deep=True)
# vtype="orf"
# drf_type="Absolute DRF"
# outf="Scripts/Data_Vis/_shap_vs_abs_drf_sgd_benomyl_pav.png"

out_snp_shap = shap_vs_drf(snp_shap_beno, snp_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_snps, vtype="snp", ylab="SNPs", drf_type="Mean absolute DRF", outf="Scripts/Data_Vis/_shap_vs_abs_drf_sgd_benomyl_snp.pdf")
out_pav_shap = shap_vs_drf(pav_shap_beno, pav_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_orfs, vtype="orf", ylab="PAVs", drf_type="Mean absolute DRF", outf="Scripts/Data_Vis/_shap_vs_abs_drf_sgd_benomyl_pav.pdf")
out_cnv_shap = shap_vs_drf(cnv_shap_beno, cnv_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_orfs, vtype="orf", ylab="CNVs", drf_type="Mean absolute DRF", outf="Scripts/Data_Vis/_shap_vs_abs_drf_sgd_benomyl_cnv.pdf")
out_snp_shap2 = shap_vs_drf(snp_shap_beno, snp_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_snps, vtype="snp", ylab="SNPs", drf_type="Mean DRF", outf="Scripts/Data_Vis/_shap_vs_drf_sgd_benomyl_snp.pdf")
out_pav_shap2 = shap_vs_drf(pav_shap_beno, pav_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_orfs, vtype="orf", ylab="PAVs", drf_type="Mean DRF", outf="Scripts/Data_Vis/_shap_vs_drf_sgd_benomyl_pav.pdf")
out_cnv_shap2 = shap_vs_drf(cnv_shap_beno, cnv_shap_nonb, drf_beno_summary[("Benomyl_abs", "mean")],
			map_orfs, vtype="orf", ylab="CNVs", drf_type="Mean DRF", outf="Scripts/Data_Vis/_shap_vs_drf_sgd_benomyl_cnv.pdf")

out = pd.DataFrame({("SNP", "mean absolute DRF", "mean of mean absolute SHAP"): out_snp_shap,
				    ("PAV", "mean absolute DRF", "mean of mean absolute SHAP"): out_pav_shap,
					("CNV", "mean absolute DRF", "mean of mean absolute SHAP"): out_cnv_shap,
				    ("SNP", "mean DRF", "mean of mean absolute SHAP"): out_snp_shap2,
				    ("PAV", "mean DRF", "mean of mean absolute SHAP"): out_pav_shap2,
				    ("CNV", "mean DRF", "mean of mean absolute SHAP"): out_cnv_shap2}).T
out.columns = ["num SGD benomyl benchmark genes", "r_benchmark", "pval_benchmark",
			   "num non-benchmark genes", "r_non_benchmark", "pval_non_benchmark"]
out.index.names = ["variant type", "DRF value", "SHAP value"]
out.to_csv("Scripts/Data_Vis/_shap_vs_drf_sgd_benomyl_pearson.csv")				  


# ------------------------------------------------------------------------------
# 3. Analysis of SHAP interactions from the Figure 4A models
# ------------------------------------------------------------------------------
#--------------------------------------#
# a. Integrated feature tables derived
#    from the Fig. 4A feature sets
#--------------------------------------#
# Figure 4A feature lists
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
snp = pd.read_csv(f"{d}/Features_bench_plus_important_non_bench_YPDBENOMYL500_snp.txt", header=None)
pav = pd.read_csv(f"{d}/Features_important_non_bench_plus_bench_genes_YPDBENOMYL500_pav.txt", header=None)
cnv = pd.read_csv(f"{d}/Features_important_non_bench_plus_bench_genes_YPDBENOMYL500_cnv.txt", header=None)

# Load the genotype matrices to make feature tables for integrated models
snp_geno = pd.read_csv("Scripts/Data_Vis/S2_File.csv")
pav_geno = pd.read_excel("Scripts/Data_Vis/S5_File.xlsx", engine="openpyxl")
cnv_geno = pd.read_excel("Scripts/Data_Vis/S6_File.xlsx", engine="openpyxl")
snp_geno.set_index("ID", inplace=True)
pav_geno.set_index("ID", inplace=True)
cnv_geno.set_index("ID", inplace=True)
pav_geno = pav_geno.loc[snp_geno.index, :] # ensure isolates are in the same order
cnv_geno = cnv_geno.loc[snp_geno.index, :]

# Create the integrated feature tables
snp_pav_df = pd.concat([snp_geno[snp[0].tolist()], pav_geno[pav[0].tolist()]], axis=1)
snp_cnv_df = pd.concat([snp_geno[snp[0].tolist()], cnv_geno[cnv[0].tolist()]], axis=1)
pav_cnv_df = pd.concat([pav_geno[pav[0].tolist()].add_suffix(";PAV", axis=1), 
						cnv_geno[cnv[0].tolist()].add_suffix(";CNV", axis=1)], axis=1)
snp_pav_cnv_df = pd.concat([snp_geno[snp[0].tolist()],
							pav_geno[pav[0].tolist()].add_suffix(";PAV", axis=1), 
							cnv_geno[cnv[0].tolist()].add_suffix(";CNV", axis=1)], axis=1)
d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_fig4a_rf"
snp_pav_df.to_csv(f"{d}/integrated_snp_pav_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")
snp_cnv_df.to_csv(f"{d}/integrated_snp_cnv_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")
pav_cnv_df.to_csv(f"{d}/integrated_pav_cnv_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")
snp_pav_cnv_df.to_csv(f"{d}/integrated_snp_pav_cnv_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")

# Save the single feature type tables in the same directory too
snp_geno[snp[0].tolist()].to_csv(f"{d}/snp_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")
pav_geno[pav[0].tolist()].to_csv(f"{d}/pav_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")
cnv_geno[cnv[0].tolist()].to_csv(f"{d}/cnv_important_non_bench_plus_bench_genes_YPDBENOMYL500.txt")

#--------------------------------------#
# b. SHAP interaction values were estimated
# c. Correlate them to Costanzo GI epsilon
#--------------------------------------#
# read SNP and ORF gene maps
map_snps = pd.read_excel("Scripts/Data_Vis/S8_File.xlsx", engine="openpyxl")
map_orfs = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl")
map_snps_dict = map_snps[["snp", "gene"]].set_index("snp").to_dict(orient="dict")
map_orfs_dict = map_orfs[["orf", "gene"]].set_index("orf").to_dict(orient="dict")

# read in SHAP interaction values
path = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_fig4a_rf"

def load_shap_itx(shap_dir, map_orfs_dict=map_orfs_dict, map_snps_dict=map_snps_dict):
	if shap_dir in ["snp", "pav", "cnv"]:
		file = f"{path}/{shap_dir}/total_shap_interaction_scores_{shap_dir}_fig4a_benomyl_500ugml_interaction_summed.txt"
	else:
		file = f"{path}/{shap_dir}/total_shap_interaction_scores_integrated_{shap_dir}_fig4a_benomyl_500ugml_interaction_summed.txt"
	shap_int = pd.read_csv(file, sep="\t")
	# add column to indicate which model the scores came from
	shap_int["Model"] = shap_dir.replace("_", " + ").upper()
	# add two columns that correspond to the feature data types
	if shap_dir in ["snp", "pav", "cnv"]:
		shap_int["Feature1_Data"] = shap_int["Model"]
		shap_int["Feature2_Data"] = shap_int["Model"]
	elif shap_dir == "snp_pav":
		shap_int["Feature1_Data"] = shap_int.apply(
			lambda x: "SNP" if "chromosome" in x["Feature1"] else "PAV", axis=1)
		shap_int["Feature2_Data"] = shap_int.apply(
			lambda x: "SNP" if "chromosome" in x["Feature2"] else "PAV", axis=1)
	elif shap_dir == "snp_cnv":
		shap_int["Feature1_Data"] = shap_int.apply(
			lambda x: "SNP" if "chromosome" in x["Feature1"] else "CNV", axis=1)
		shap_int["Feature2_Data"] = shap_int.apply(
			lambda x: "SNP" if "chromosome" in x["Feature2"] else "CNV", axis=1)
	else:
		shap_int["Feature1_Data"] = shap_int.apply(
			lambda x: x["Feature1"].split(";")[-1] if ";" in x["Feature1"] else "SNP", axis=1)
		shap_int["Feature2_Data"] = shap_int.apply(
			lambda x: x["Feature2"].split(";")[-1] if ";" in x["Feature2"] else "SNP", axis=1)
	# Clean the PAV/CNV feature IDs
	def fix_orf_names(orf_pd_series, axis=0):
		'''Prepare ORF identifiers for mapping to genes:'''
		orf_pd_series = orf_pd_series.apply(lambda x: re.sub("^X", "", x))
		orf_pd_series = orf_pd_series.apply(lambda x: re.sub(";PAV$", "", x))
		orf_pd_series = orf_pd_series.apply(lambda x: re.sub(";CNV$", "", x))
		orf_pd_series = orf_pd_series.apply(lambda x: re.sub("\.", "-", x))
		return orf_pd_series
	shap_int.Feature1 = fix_orf_names(shap_int.Feature1)
	shap_int.Feature2 = fix_orf_names(shap_int.Feature2)
	# Map features to genes
	def feature2gene(feature):
		if feature.endswith("-"): # for 45-PUT_IRON_TRASP-R-14849-
			feature = feature[:-1]
		try:
			return map_orfs_dict["gene"][feature]
		except KeyError:
			try:
				return map_snps_dict["gene"][feature]
			except KeyError:
				return feature
	shap_int["Gene1"] = shap_int.Feature1.apply(lambda x: feature2gene(x))
	shap_int["Gene2"] = shap_int.Feature2.apply(lambda x: feature2gene(x))
	return shap_int

all_shap_itx = []
for shap_dir in ["snp", "pav", "cnv", "snp_pav", "snp_cnv", "pav_cnv", "snp_pav_cnv"]:
	shap_itx = load_shap_itx(shap_dir)
	print(shap_itx)
	all_shap_itx.append(shap_itx)
	gc.collect()

ts17_df = pd.concat(all_shap_itx, axis=0)
ts17_df.shape # (1633971, 633)
del all_shap_itx

# Save Table S17
out_path = "Scripts/Data_Vis/Section_6/shap_interaction/Table_S17_benomyl_fig4a_RF_models_SHAP_interactions.txt"
n_chunks = 16
chunk_size = int(np.ceil(len(ts17_df) / n_chunks))
for i, start in enumerate(range(0, len(ts17_df), chunk_size), start=1):
	chunk = ts17_df.iloc[start:start + chunk_size]
	chunk.to_csv(out_path, sep="\t", index=False, header=(i == 1), mode="w" if i == 1 else "a")
	print(f"Finished writing chunk {i}/{n_chunks}")

# Determine the overlap between SHAP interactions and Costanzo et al., 2021-----
ts17_df = dt.fread(out_path).to_pandas()
costanzo = pd.read_excel(
	"Data/Costanzo_2021/2021_Costanzo_Data File S3_Raw interaction dataset.xlsx",
	engine="openpyxl", sheet_name="Genome-scale_Benomyl")
costanzo.shape # (100237, 25)
costanzo.head()
costanzo.isna().sum() # no missing values in condition, reference condition columns

# extract gene systematic identifiers from query_orf and array_orf columns
costanzo["Gene1"] = costanzo["query_orf"].str.split("_").str[0]
costanzo["Gene2"] = costanzo["array_orf"].str.split("_").str[0]

# sort gene pairs alphabetically to facilitate comparison between datasets
def sort_genes(row, gene1_col, gene2_col):
	gene1 = row[gene1_col]
	gene2 = row[gene2_col]
	sorted_genes = sorted([gene1, gene2])
	return pd.Series(sorted_genes)

costanzo[["Gene1_sorted", "Gene2_sorted"]] = costanzo.apply(
	lambda row: sort_genes(row, "Gene1", "Gene2"), axis=1)
ts17_df[["Gene1_sorted", "Gene2_sorted"]] = ts17_df.apply(
	lambda row: sort_genes(row, "Gene1", "Gene2"), axis=1)

# number of unique Costanzo et al. 2021 GIs
beno_30_ntwk = costanzo[["Gene1_sorted", "Gene2_sorted", "rep1_condition_epsilon",
	"rep2_condition_epsilon", "mean_condition_epsilon", "sd_condition_epsilon", "condition_p_value"]]

def get_unique_gp(df, Gene1="Gene1", Gene2="Gene2"):
    '''
    Get the unique gene pairs from the dataframe
    '''
    df_gp = df.apply(lambda x: set(
        [x[Gene1], x[Gene2]]), axis=1).values  # gene pairs
    df_gp = {frozenset(sorted(set))
             for set in df_gp}  # get unique interactions
    return df_gp

beno_30_ntwk_gp = get_unique_gp(beno_30_ntwk, "Gene1_sorted", "Gene2_sorted")
len(beno_30_ntwk_gp) # 92906

beno_30_ntwk_significant = beno_30_ntwk[beno_30_ntwk["condition_p_value"] < 0.05]
beno_30_ntwk_sig_gp = get_unique_gp(beno_30_ntwk_significant, "Gene1_sorted", "Gene2_sorted")
len(beno_30_ntwk_sig_gp) # 11979
	
# overlap between Costanzo et al. 2021 and SHAP interaction data
costanzo_dict = set(costanzo[["Gene1_sorted", "Gene2_sorted"]].apply(tuple, axis=1))
len(costanzo_dict) # 92906

shap_dict = set(ts17_df[["Gene1_sorted", "Gene2_sorted"]].apply(tuple, axis=1))
len(shap_dict) # 1212842

overlap = costanzo_dict.intersection(shap_dict)
len(overlap) # 3794 overlapping unique gene pairs

# how many of the gene pairs have significant p-values
significant = costanzo.loc[(costanzo["Gene1_sorted"].isin(ts17_df["Gene1_sorted"])) & (
	costanzo["Gene2_sorted"].isin(ts17_df["Gene2_sorted"])) & (
	costanzo["condition_p_value"] < 0.05),:]
significant.shape # (1267, 29)

significant[["Gene1_sorted", "Gene2_sorted"]].drop_duplicates(keep="first").shape # (1179, 2), has 88 duplicates.
significant.loc[significant.duplicated(subset=["Gene1_sorted", "Gene2_sorted"],\
									   keep=False)].sort_values(by=["Gene1_sorted", "Gene2_sorted"]) # 165 rows

# how many of the significant benomyl interactions overlap with the individual model SHAP interaction sets?
models = ["SNP", "PAV", "CNV", "SNP + PAV", "SNP + CNV", "PAV + CNV", "SNP + PAV + CNV"]
for model in models:
	model_gp = get_unique_gp(ts17_df[ts17_df.Model==model])
	print(model, len(model_gp.intersection(beno_30_ntwk_sig_gp)))

# SNP 59
# PAV 0
# CNV 2
# SNP + PAV 65
# SNP + CNV 79
# PAV + CNV 4
# SNP + PAV + CNV 400 gene pairs

# determine the range of genetic interaction scores (epsilon) from Costanzo et al. 2021
# for the overlapping gene pairs
overlap_df = costanzo[
	costanzo[["Gene1_sorted", "Gene2_sorted"]].apply(tuple, axis=1).isin(overlap)]
overlap_df.shape # (4196, 29)

# calculate summary statistics on the Costanzo et al. 2021 columns
overlap_df_summary = overlap_df.groupby(["Gene1_sorted", "Gene2_sorted"]).describe()
overlap_df_summary = overlap_df_summary.stack(level=0, future_stack=True)
overlap_df_summary.index.names = ["Gene1_sorted", "Gene2_sorted", "Column"]
overlap_df_summary.to_csv(
	"Scripts/Data_Vis/_costanzo_fig4a_shap_interaction_overlap_summary.csv")

# Estimate the correlation between SHAP interaction scores and epsilon values---
import statsmodels.api as sm
from glob import glob
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform

# load the genetic distance matrices
dist_snp = pd.read_excel("Scripts/Data_Vis/S13_File.xlsx", engine="openpyxl")
dist_pav = pd.read_excel("Scripts/Data_Vis/S14_File.xlsx", engine="openpyxl")
dist_snp.set_index("ID", inplace=True)
dist_pav.set_index(dist_pav.columns[0], inplace=True)

# load the SNP & PAV matrices with S288C encodings
snp = dt.fread("Scripts/Data_Vis/S12_File.csv").to_pandas()
snp.set_index("ID", inplace=True)
# fix the data types
cols_tofix = snp.select_dtypes("object").columns.tolist()
from setuptools._distutils.util import strtobool
for col in cols_tofix:
	snp[col] = snp[col].apply(lambda x: strtobool(x) if isinstance(x, str) else x)
assert snp.select_dtypes("object").empty
#snp.to_csv("Scripts/Data_Vis/S12_File_fixed_dtypes.csv")
pav = pd.read_csv("Data/Peter_2018/ORFs_pres_abs_with_S288C.csv", index_col=0)
pav.index = pav.index.str.replace("^X", "", regex=True).str.replace("\.", "-", regex=True)
pav = pav.T
pav = pav.loc[snp.index, :] # ensure isolates are in the same order
pav = pav.astype(int)

# subset the benomyl genetic interactions from Costanzo et al. 2021
ben_gi = costanzo[["Gene1_sorted", "Gene2_sorted", "mean_condition_epsilon", "condition_p_value"]].dropna()
ben_gi["Gene_Pair"] = ben_gi[["Gene1_sorted", "Gene2_sorted"]].apply(tuple, axis=1)
# gi_df = ben_gi.copy(deep=True)

# functions to calculate the correlation between local SHAP interaction values 
# and epsilon values of the overlapping benomyl gene pairs
def fix_sample_names(rho, dist_snp, dist_pav):
	"""Sample names do not match between the shap interaction matrices and
	the genetic distance matrices for samples that begin with SACE_"""
	new_idx = []
	counter = 0
	for idx in rho.index:
		if idx in dist_snp["S288C"].index:
			new_idx.append(idx)
		elif idx in ["Interaction", "mean_condition_epsilon", "condition_p_value"]:
			new_idx.append(idx)
		else:
			if f"SACE_{idx}" in dist_snp["S288C"].index:
				new_idx.append(f"SACE_{idx}")
				counter += 1
	tmp = pd.merge(rho, dist_snp["S288C"], left_index=True, right_index=True, how="left")
	tmp = pd.merge(tmp, dist_pav["S288C"], left_index=True, right_index=True, how="left")
	try:
		assert counter == tmp.isna().sum()["S288C_x"] - 1 # Interaction
		assert counter == tmp.isna().sum()["S288C_y"] - 1
	except AssertionError:
		assert counter == tmp.isna().sum()["S288C_x"] - 3 # Interaction, mean_condition_epsilon, and condition_p_value
		assert counter == tmp.isna().sum()["S288C_y"] - 3
	return new_idx


def add_genetic_distances(rho, name, dist_snp, dist_pav):
	"""Add columns of SNP-based and PAV-based genetic distances between the
	isolates and S288C (Euclidean distances). These distances were calculated
	using the 118k SNPs and 7.7k PAVs"""
	rho.name = name
	rho.index = fix_sample_names(rho, dist_snp, dist_pav)
	rho = pd.merge(rho, dist_snp["S288C"], left_index=True, right_index=True, how="left")
	rho = pd.merge(rho, dist_pav["S288C"], left_index=True, right_index=True, how="left")
	rho.rename(columns={"S288C_x": "All_SNPs_dist", "S288C_y": "All_PAVs_dist"}, inplace=True)
	return rho


def fit_line(data, x, y):
	"""Fit a line to the data and return the slope, intercept, and R² value"""
	data_nona = data.dropna(subset=[x, y])
	X = sm.add_constant(data_nona[x])
	model = sm.OLS(data_nona[y], X).fit()
	m = model.params[1]
	b = model.params[0]
	adj_r2 = model.rsquared_adj
	# return f"slope: {m:.4f}, intercept: {b:.4f}, adj R²: {adj_r2:.4f}"
	return {x: {"y": y, "m": m, "b": b, "adj_r2": adj_r2}}


def calc_local_shap_itx_corr(ts17_df, shap_dir, wdir, gi_df, dist_snp, dist_pav, snp, pav, save_name=None):
	# Subset the relevant SHAP interaction matrix for the given model type
	shap_itx = ts17_df.loc[ts17_df["Model"] == shap_dir.replace("_", " + ").upper(),:].copy(deep=True)
	shap_itx.index = pd.MultiIndex.from_frame(shap_itx.iloc[:,:2])
	shap_itx.drop(columns=shap_itx.columns[:2], inplace=True)
	
	# 1. Determine overlap with significant Costanzo et al. 2021 benomyl genetic interactions
	shap_itx["Gene_Pair"] = shap_itx[["Gene1_sorted", "Gene2_sorted"]].apply(tuple, axis=1)
	overlap = gi_df.loc[(gi_df["Gene_Pair"].isin(shap_itx["Gene_Pair"])) & (
		gi_df["condition_p_value"] < 0.05),:]
	shap_itx_overlap = shap_itx.loc[shap_itx["Gene_Pair"].isin(overlap["Gene_Pair"]),:]
	if len(shap_itx_overlap) == 0:
		return None
	
	# 2. Calculate correlations between epsilon and local SHAP interaction values
	# check for duplicate gene pairs before setting Gene_Pair as index
	if overlap.Gene_Pair.nunique() == len(overlap):
		overlap.set_index("Gene_Pair", inplace=True)
	else:
		print("Duplicate gene pairs in overlap, cannot set index to Gene_Pair. Calculating mean.")
		overlap = overlap[["mean_condition_epsilon", "condition_p_value", "Gene_Pair"]].groupby("Gene_Pair").mean()
	if shap_itx_overlap.Gene_Pair.nunique() == len(shap_itx_overlap):
		shap_itx_overlap.reset_index(inplace=True)
		shap_itx_overlap.set_index("Gene_Pair", inplace=True)
		# ensure indices are the same order for both dataframes
		shap_itx_overlap = shap_itx_overlap.loc[overlap.index,:]
		# calculate spearman and pearson correlations
		#for top_itx in []: # does correlation increase when only the most important overlapping gene pairs are considered?
		rho = shap_itx_overlap.select_dtypes("number").corrwith(overlap["mean_condition_epsilon"], method="spearman", axis=0)
		pcc = shap_itx_overlap.select_dtypes("number").corrwith(overlap["mean_condition_epsilon"], method="pearson", axis=0)
	else:
		print("Duplicate gene pairs in shap_itx_overlap, cannot set index to Gene_Pair.")
		shap_itx_overlap = pd.merge(shap_itx_overlap.reset_index(), overlap,
			left_on="Gene_Pair", right_on="Gene_Pair", how="left")
		# calculate spearman and pearson correlations
		rho = shap_itx_overlap.select_dtypes("number").corrwith(shap_itx_overlap["mean_condition_epsilon"], method="spearman", axis=0)
		pcc = shap_itx_overlap.select_dtypes("number").corrwith(shap_itx_overlap["mean_condition_epsilon"], method="pearson", axis=0)
	
	# Load the feature importance values from the corresponding RF model
	if shap_dir in ["snp", "pav", "cnv"]:
		imp = pd.read_csv(f"{wdir}/{shap_dir}/SHAP_values_sorted_average_{shap_dir}_fig4a_benomyl_500ugml_training.txt", sep="\t", header=0, index_col=0)
	elif shap_dir in ["snp_pav", "snp_cnv", "pav_cnv", "snp_pav_cnv"]:
		imp = pd.read_csv(f"{wdir}/{shap_dir}/SHAP_values_sorted_average_integrated_{shap_dir}_fig4a_benomyl_500ugml_training.txt", sep="\t", header=0, index_col=0)
	imp = imp.loc[imp["0"] > 0, :].sort_values("0", ascending=False)
	
	# 3. Estimate genetic distances based on top feature subsets
	for n_imp in [700, 600, 500, 400, 300, 200, 100, 50, 10]:
		# Add columns of genetic distances between the isolates and S288C based on the SNPs and PAVs
		if n_imp == 700:
			rho = add_genetic_distances(rho, "rho", dist_snp, dist_pav)
			pcc = add_genetic_distances(pcc, "r", dist_snp, dist_pav)
		if n_imp > len(imp):
			continue
		if shap_dir == "snp":
			subset = snp[imp.iloc[:n_imp].index.tolist()]
		elif shap_dir in ["pav", "cnv"]:
			pav_cols = imp.iloc[:n_imp].index.tolist()
			pav_cols = [re.sub("^X", "", col) for col in pav_cols]
			pav_cols = [re.sub("\.", "-", col) for col in pav_cols]
			subset = pav[pav_cols]
		elif shap_dir in ["snp_pav", "snp_cnv", "snp_pav_cnv"]:
			snp_cols = [col for col in imp.iloc[:n_imp].index if "chromosome" in col]
			pav_cols = [col for col in imp.iloc[:n_imp].index if ";" in col]
			snp_scaled = StandardScaler().fit_transform(snp[snp_cols])
			subset = snp_scaled.copy()
			if len(pav_cols) > 0:
				pav_cols = set([re.sub(";CNV", "", re.sub(";PAV", "", re.sub("\.", "-", re.sub("^X", "", c)))) for c in pav_cols])
				pav_scaled = StandardScaler().fit_transform(pav[list(pav_cols)])
				subset = np.hstack([subset, pav_scaled])
			subset = pd.DataFrame(subset, index=snp.index)
		elif shap_dir in ["pav_cnv"]:
			cols = [re.sub(";CNV", "", re.sub(";PAV", "", re.sub("\.", "-", re.sub("^X", "", c)))) for c in imp.iloc[:n_imp].index]
			subset = pav[cols]
		# calculate Euclidean distance
		dist = pd.DataFrame(squareform(pdist(subset, metric="euclidean")), index=subset.index, columns=subset.index)
		# insert column into rho
		rho.insert(3, f"{shap_dir}_top{n_imp}_based_dist", rho.apply(lambda row: dist.loc[row.name, "S288C"] if row.name in dist.index else np.nan, axis=1))
	
	# 4. Save results and plot histogram of correlations
	assert sum(rho.index==pcc.index) == len(rho)
	rho.insert(1, "r", pcc["r"])
	rho.to_csv(save_name.replace(".pdf", "_with_genet_dist.csv"))
	
	# calculate regression lines
	results = {}
	for col in rho.columns[2:]:
		out = fit_line(rho, col, "r")
		results.update(out)
		del out
	pd.DataFrame.from_dict(results, orient="index").to_csv(save_name.replace(".pdf", "_regression_results_r.csv"))
	
	fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
	try:
		ax[0].hist(rho.drop(["Interaction", "mean_condition_epsilon", "condition_p_value"])["rho"], bins=50)
		ax[1].hist(pcc.drop(["Interaction", "mean_condition_epsilon", "condition_p_value"])["r"], bins=50)
	except KeyError:
		ax[0].hist(rho.drop("Interaction")["rho"], bins=50)
		ax[1].hist(pcc.drop("Interaction")["r"], bins=50)
	ax[0].set_xlabel("Spearman's rho", fontsize=8)
	ax[0].set_ylabel("Number of isolates", fontsize=8)
	ax[0].tick_params(axis="both", labelsize=7)
	ax[1].set_xlabel("Pearson's r", fontsize=8)
	ax[1].set_ylabel("Number of isolates", fontsize=8)
	ax[1].tick_params(axis="both", labelsize=7)
	plt.suptitle(f"Correlation btwn local SHAP interaction values ({len(shap_itx_overlap)} feature pairs) &\n Costanzo benomyl GI scores ({len(overlap)} gene pairs); ({shap_dir} model)", fontsize=8)
	plt.tight_layout()
	plt.savefig(save_name)
	plt.close("all")
	
	del shap_itx, shap_itx_overlap, overlap, rho, pcc, results, imp
	gc.collect()


wdir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/benomyl_shap_int_fig4a_rf"
shap_dirs = ["snp", "pav", "cnv", "snp_pav", "snp_cnv", "pav_cnv", "snp_pav_cnv"]
for shap_dir in shap_dirs:
	save_name = f"Scripts/Data_Vis/_fig4a_local_shap_interaction_corrs_benomyl_costanzo_mean_{shap_dir}.pdf"
		# before, I calculated the median for the gene pairs in the 'overlap' dataframe
		# the file using median GI scores is called _fig4a_local_shap_interaction_corrs_benomyl_costanzo_{shap_dir}.pdf
	calc_local_shap_itx_corr(ts17_df, shap_dir, wdir, ben_gi, dist_snp, dist_pav, snp, pav, save_name)

# Plot regression results
# files = glob("_local_*_regression_results_rho.csv")
files = glob("Scripts/Data_Vis/_fig4a_local_*mean*_regression_results_r.csv")
reg_results = pd.concat([pd.read_csv(file).assign(
	Model=file.split("costanzo_")[1].split("_reg")[0]) for file in files], axis=0)
reg_results["Genetic_distance_type"] = reg_results["Unnamed: 0"].str.replace("^([snpavc_])+", "", regex=True)
reg_results = reg_results.pivot(index="Model", columns="Genetic_distance_type", values="adj_r2")
reg_results = reg_results[["All_PAVs_dist", "All_SNPs_dist"] + [f"top{n}_based_dist" for n in [700, 600, 500, 400, 300, 200, 100, 50, 10]]]

sns.heatmap(reg_results.T, annot=True, cmap="RdBu_r", vmin=-0.002, vmax=0.32)
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/_fig4a_local_shap_interaction_corrs_benomyl_costanzo_mean_regression_results_r_heatmap.pdf")
# the file using the median GI scores is called _fig4a_local_shap_interaction_corrs_benomyl_costanzo_regression_results_r_heatmap.pdf
plt.close("all")

# What are the mean and standard deviations of the rho and r values between local SHAP interaction values and epsilon for the different models?
for shap_dir in shap_dirs:
	print(shap_dir)
	res = pd.read_csv(f"Scripts/Data_Vis/_fig4a_local_shap_interaction_corrs_benomyl_costanzo_{shap_dir}_with_genet_dist.csv")
	res[["rho", "r"]].describe().to_dict()

""" using the mean GI scores:
snp
{'rho': {'count': 626.0, 'mean': -0.007388108891349521, 'std': 0.11271769896564873, 'min': -0.2165984804208065, '25%': -0.0820572764465225, '50%': -0.03702513150204555, '75%': 0.0469316189362945, 'max': 0.2967855055523086}, 'r': {'count': 626.0, 'mean': -0.009216861906812656, 'std': 0.09184294908237517, 'min': -0.2745529677918091, '25%': -0.07355768764257002, '50%': -0.04116199718894885, '75%': 0.0363301753926538, 'max': 0.2366517446475052}}
pav
{'rho': {'count': 0.0, 'mean': nan, 'std': nan, 'min': nan, '25%': nan, '50%': nan, '75%': nan, 'max': nan}, 'r': {'count': 0.0, 'mean': nan, 'std': nan, 'min': nan, '25%': nan, '50%': nan, '75%': nan, 'max': nan}}
cnv
{'rho': {'count': 626.0, 'mean': -0.645367412140575, 'std': 0.7644831618762685, 'min': -1.0, '25%': -1.0, '50%': -1.0, '75%': -1.0, 'max': 1.0}, 'r': {'count': 626.0, 'mean': -0.645367412140575, 'std': 0.7644831618762685, 'min': -1.0, '25%': -1.0, '50%': -1.0, '75%': -1.0, 'max': 1.0}}
snp_pav
{'rho': {'count': 565.0, 'mean': -0.18687758941957217, 'std': 0.26713899470069, 'min': -0.5603715170278638, '25%': -0.4241486068111454, '50%': -0.26109391124871, '75%': 0.0257997936016511, 'max': 0.5067079463364293}, 'r': {'count': 565.0, 'mean': -0.11910896562992336, 'std': 0.15200974793378572, 'min': -0.4544154827158169, '25%': -0.2030012379367969, '50%': -0.1636579572107573, '75%': -0.0090888684111931, 'max': 0.3392008841734866}}
snp_cnv
{'rho': {'count': 560.0, 'mean': 0.08419450728959732, 'std': 0.17548962257603865, 'min': -0.6300751879699248, '25%': 0.02819548872180445, '50%': 0.1353383458646616, '75%': 0.18796992481203, 'max': 0.4285714285714286}, 'r': {'count': 560.0, 'mean': 0.20003809623672542, 'std': 0.24031922564499977, 'min': -0.5801079838550148, '25%': 0.0436397911619838, '50%': 0.2993383803352372, '75%': 0.3644688663818251, 'max': 0.548882570376176}}
pav_cnv
{'rho': {'count': 626.0, 'mean': -0.4025559105431309, 'std': 0.031897599737434515, 'min': -0.7999999999999999, '25%': -0.3999999999999999, '50%': -0.3999999999999999, '75%': -0.3999999999999999, 'max': -0.3999999999999999}, 'r': {'count': 626.0, 'mean': -0.8472535150512799, 'std': 0.02400944360496991, 'min': -0.961109807320122, '25%': -0.8471530457693914, '50%': -0.8471530457693914, '75%': -0.841475959994029, 'max': -0.7096202498870641}}
snp_pav_cnv
{'rho': {'count': 433.0, 'mean': 0.03627645217952414, 'std': 0.06362519229798458, 'min': -0.0785730291597561, '25%': 0.0018515285579862, '50%': 0.0363936474767917, '75%': 0.0659662794170369, 'max': 0.9999999999999998}, 'r': {'count': 433.0, 'mean': 0.07077656353431802, 'std': 0.06100434425044568, 'min': -0.0764195885886929, '25%': 0.0481234733153076, '50%': 0.074815111232697, '75%': 0.0938742011839898, 'max': 1.0}}
"""

""" using the median GI scores:
snp
{'rho': {'count': 626.0, 'mean': -0.007388108891349521, 'std': 0.11271769896564873, 'min': -0.2165984804208065, '25%': -0.0820572764465225, '50%': -0.03702513150204555, '75%': 0.0469316189362945, 'max': 0.2967855055523086},
 'r': {'count': 626.0, 'mean': -0.009216861906812656, 'std': 0.09184294908237517, 'min': -0.2745529677918091, '25%': -0.07355768764257002, '50%': -0.04116199718894885, '75%': 0.0363301753926538, 'max': 0.2366517446475052}}
pav
{'rho': {'count': 0.0, 'mean': nan, 'std': nan, 'min': nan, '25%': nan, '50%': nan, '75%': nan, 'max': nan},
 'r': {'count': 0.0, 'mean': nan, 'std': nan, 'min': nan, '25%': nan, '50%': nan, '75%': nan, 'max': nan}}
cnv
{'rho': {'count': 626.0, 'mean': -0.645367412140575, 'std': 0.7644831618762685, 'min': -1.0, '25%': -1.0, '50%': -1.0, '75%': -1.0, 'max': 1.0},
 'r': {'count': 626.0, 'mean': -0.645367412140575, 'std': 0.7644831618762685, 'min': -1.0, '25%': -1.0, '50%': -1.0, '75%': -1.0, 'max': 1.0}}
snp_pav
{'rho': {'count': 565.0, 'mean': -0.18687758941957217, 'std': 0.26713899470069, 'min': -0.5603715170278638, '25%': -0.4241486068111454, '50%': -0.26109391124871, '75%': 0.0257997936016511, 'max': 0.5067079463364293},
 'r': {'count': 565.0, 'mean': -0.11910896562992336, 'std': 0.15200974793378572, 'min': -0.4544154827158169, '25%': -0.2030012379367969, '50%': -0.1636579572107573, '75%': -0.0090888684111931, 'max': 0.3392008841734866}}
snp_cnv
{'rho': {'count': 560.0, 'mean': 0.08419450728959732, 'std': 0.17548962257603865, 'min': -0.6300751879699248, '25%': 0.02819548872180445, '50%': 0.1353383458646616, '75%': 0.18796992481203, 'max': 0.4285714285714286},
 'r': {'count': 560.0, 'mean': 0.20003809623672542, 'std': 0.24031922564499977, 'min': -0.5801079838550148, '25%': 0.0436397911619838, '50%': 0.2993383803352372, '75%': 0.3644688663818251, 'max': 0.548882570376176}}
pav_cnv
{'rho': {'count': 626.0, 'mean': -0.4025559105431309, 'std': 0.031897599737434515, 'min': -0.7999999999999999, '25%': -0.3999999999999999, '50%': -0.3999999999999999, '75%': -0.3999999999999999, 'max': -0.3999999999999999},
 'r': {'count': 626.0, 'mean': -0.8472535150512799, 'std': 0.02400944360496991, 'min': -0.961109807320122, '25%': -0.8471530457693914, '50%': -0.8471530457693914, '75%': -0.841475959994029, 'max': -0.7096202498870641}}
snp_pav_cnv
{'rho': {'count': 433.0, 'mean': 0.03627645217952414, 'std': 0.06362519229798458, 'min': -0.0785730291597561, '25%': 0.0018515285579862, '50%': 0.0363936474767917, '75%': 0.0659662794170369, 'max': 0.9999999999999998},
 'r': {'count': 433.0, 'mean': 0.07077656353431802, 'std': 0.06100434425044568, 'min': -0.0764195885886929, '25%': 0.0481234733153076, '50%': 0.074815111232697, '75%': 0.0938742011839898, 'max': 1.0}}
"""

# ------------------------------------------------------------------------------
"""4. Examining CNVs (S12 Table)
- How do the effects of CNV features (SHAP values) differ across environments?
	- A feature's CNVs across isolates--> compare to (pearson)--> local SHAP vals
		- calculate for all 625 training isolates combined
		- calculate for each SHAP cluster of isolates
	- Compare differences in pearson correlations across environments
	- Correlate CNV values with fitness values across isolates (pearson)
		- How do the correlations change with median absolute feature importance?
"""
# ------------------------------------------------------------------------------

import re, gc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statsmodels.api as sm
from glob import glob
from scipy.stats import pearsonr
from matplotlib import ticker

# Subset training CNV data
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
cnv = pd.read_excel("Scripts/Data_Vis/S6_File.xlsx", engine="openpyxl")
cnv.set_index("ID", inplace=True)
cnv = cnv[~cnv.index.isin(test[0])]
cnv_var = cnv.var()

# Local SHAP values from the optimized RF models
shap_dir = "/mnt/ufs18/rs-049/glbrc_group/shiulab/kenia/yeast_project/SHAP/CNV/fs/"
shap_files = sorted(glob(f"{shap_dir}/SHAP_values_sorted_YP*"))

# Fitness data (for plotting cnv values, local shap values, and fitness)
pheno = pd.read_excel("Scripts/Data_Vis/S1_File.xlsx", engine="openpyxl")
pheno.set_index("ID", inplace=True)
pheno = pheno[~pheno.index.isin(test[0])]
pheno = pheno.loc[cnv.index,:]

# Calculate the pearson correlation between each CNV feature's values & SHAP values
def shap_subplots(top20, shap_df, save_name):
	fig, ax = plt.subplots(nrows=5, ncols=4, figsize=(12, 11), sharex=False, sharey=True)
	for i, feat_tuple in enumerate(top20.iterrows()):
		row_idx = i // 4
		col_idx = i % 4
		# create colorbar gradient
		g = ax[row_idx, col_idx].scatter(
			x=cnv[feat_tuple[0]], y=pheno[env], c=shap_df[feat_tuple[0]],
			cmap="RdBu_r", s=40, alpha=0.6, marker=".", linewidths=.8,
			edgecolors="black")
		cbar = fig.colorbar(g, ax=ax[row_idx, col_idx], shrink=.6)
		cbar.set_label(label="SHAP", y=1.15, rotation=0, labelpad=-22)
		cbar.ax.tick_params(labelsize=6)
		cbar.ax.minorticks_on()
		ax[row_idx, col_idx].set_title(feat_tuple[0], fontsize=7, fontweight="bold", fontstyle="italic")
		# ax[row_idx, col_idx].set_xscale("log")
		# ax[row_idx, col_idx].set_xticks(sorted(cnv[feat_tuple[0]].unique()))
		ax[row_idx, col_idx].get_xaxis().set_major_formatter(ticker.ScalarFormatter())
		ax[row_idx, col_idx].get_xaxis().set_minor_formatter(ticker.ScalarFormatter())
		ax[row_idx, col_idx].tick_params(axis="both", labelsize=7)
		# ax[row_idx, col_idx].tick_params(axis="y", labelsize=7)
		# ax[row_idx, col_idx].tick_params(axis="x", which="minor", labelsize=6, rotation=45)
		ax[row_idx, col_idx].set_xlabel("CNV value", fontsize=7, fontweight="bold", fontstyle="italic")
		ax[row_idx, col_idx].set_ylabel("Fitness", fontsize=7, fontweight="bold", fontstyle="italic")
		ax[row_idx, col_idx].set_box_aspect(1)
	
	plt.tight_layout()
	plt.savefig(save_name)
	plt.close('all')
	gc.collect()
	gc.collect()


def calc_shap_cnv_corr(shap_file, cnv=cnv, pheno=pheno, plot=True):
	# read in local shap values
	shap_df = pd.read_csv(shap_file, sep="\t", index_col=0)
	# ensure isolates are in the same order as in the CNV matrix
	shap_df = shap_df.loc[cnv.index, :]
	# for each feature in shap_df, calculate pearson r with CNV values
	r_res = shap_df.apply(lambda f: pearsonr(f, cnv[f.name]), axis=0).T
	r_res.columns = ["r_shap_vs_cnv", "pvalue_shap_vs_cnv"]
	# insert summary statistics of the local shap values for each feature
	r_res = r_res.merge(shap_df.describe().T.add_suffix("_shap"),
						left_index=True, right_index=True)
	r_res.insert(10, "median_abs_shap", shap_df.abs().median())
	r_res.insert(11, "mean_abs_shap", shap_df.abs().mean())
	# r_res.insert(12, "variance_cnv", r_res.apply(
	# 	lambda f: cnv[f.name].var(), axis=1))
	r_res.insert(12, "variance_cnv", [cnv[col].var() for col in shap_df.columns])
	r_res = r_res.merge(cnv[shap_df.columns].describe().T.add_suffix("_cnv"),
						left_index=True, right_index=True)
	# sort features by median absolute shap values
	r_res.sort_values(by="median_abs_shap", ascending=False)
	# get environment to subset the fitness values for plotting
	env = re.search("YP[A-Z0-9]+(?!cnv)", shap_file).group()
	# add the correlation between CNV values and fitness into r_res
	r_res2 = shap_df.apply(lambda f: pearsonr(cnv[f.name], pheno[env]), axis=0).T
	r_res2.columns = ["r_cnv_vs_fitness", "pvalue_cnv_vs_fitness"]
	r_res = pd.concat([r_res, r_res2], ignore_index=False, axis=1)
	# plot the cnv values vs local SHAP values scatterplots of the top 20 features
	if plot:
		top20 = r_res.iloc[:20,:]
		shap_subplots(top20, shap_df, f"Scripts/Data_Vis/_cnv_{env}_shap_scatter.pdf")
	# clear memory
	del r_res2, shap_df, env
	gc.collect()
	return r_res


target_envs = ["YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"]
results = []
for shap_file in shap_files:
	env = re.search("YP[A-Z0-9]+(?!cnv)", shap_file).group()
	if env not in target_envs:
		env_r_res = calc_shap_cnv_corr(shap_file, plot=False)
	if env in target_envs:
		# calculate the correlations between cnv values, local SHAP values, and fitness
		env_r_res = calc_shap_cnv_corr(shap_file, plot=False)
	# How do the correlations change with median absolute feature importance?
	env_r_res["median_abs_shap_bins"] = pd.qcut(env_r_res.median_abs_shap, q=20, retbins=True)[0]
	out = env_r_res.groupby("median_abs_shap_bins", observed=False)[["variance_cnv",
		"r_shap_vs_cnv", "pvalue_shap_vs_cnv", "r_cnv_vs_fitness", "pvalue_cnv_vs_fitness"]].describe()
	env_r_res["environment"] = env
	results.append(env_r_res)
	# save env_r_res and out into an excel file
	# with pd.ExcelWriter(f"Scripts/Data_Vis/_cnv_{env}_shap_corrwith_cnv_and_fitness.xlsx", engine="openpyxl") as writer:
	# 	env_r_res.to_excel(writer, sheet_name="shap_corr_summary", index=True)
	# 	out.to_excel(writer, sheet_name="shap_corr_summary_binned", index=True)
	# 	env_r_res[['r_shap_vs_cnv', 'median_abs_shap', 'mean_abs_shap', 'variance_cnv',
	# 			   'mean_cnv', '50%_cnv', 'max_cnv', 'r_cnv_vs_fitness']].corr().to_excel(
	# 		writer, sheet_name="shap_corr_summary_comparisons", index=True)
	del env_r_res, out


results_df = pd.concat(results, axis=0)

with pd.ExcelWriter(f"Scripts/Data_Vis/_cnv_vs_fitness_pearson.xlsx", engine="openpyxl") as writer:
	results_df.to_excel(writer, sheet_name="cnv_vs_fitness_all_envs", index=True)
	results_df.groupby("environment")[["r_cnv_vs_fitness", "pvalue_cnv_vs_fitness"]].describe().to_excel(
		writer, sheet_name="summary statistics", index=True) # S14 Table

results_df.groupby("environment")[["r_cnv_vs_fitness", "pvalue_cnv_vs_fitness"]].aggregate(["mean", "std", "min", "max"])
#                       r_cnv_vs_fitness                               pvalue_cnv_vs_fitness                                  
#                                   mean       std       min       max                  mean       std           min       max
# environment                                                                                                                 
# YPACETATE                     0.001162  0.097582 -0.093477  0.223448              0.323529  0.381611  1.636838e-08  0.977106
# YPD14                        -0.075967  0.133716 -0.310318  0.219432              0.105984  0.244577  2.033989e-15  0.990387
# YPD40                        -0.033876  0.216792 -0.260619  0.252865              0.018705  0.074607  3.662590e-11  0.298480
# YPD42                         0.049867  0.161522 -0.248786  0.298752              0.092003  0.199861  2.370238e-14  0.972746
# YPD6AU                       -0.015802  0.153078 -0.304373  0.335742              0.159851  0.275572  6.218272e-18  0.998713
# YPDANISO10                   -0.025098  0.167891 -0.304130  0.266187              0.100365  0.227654  7.671055e-15  0.976281
# YPDANISO20                    0.008931  0.217690 -0.391218  0.417515              0.095627  0.211557  9.284337e-28  0.912037
# YPDANISO50                    0.064372  0.167427 -0.313668  0.383518              0.135914  0.249781  2.476241e-23  0.988522
# YPDBENOMYL200                 0.036781  0.126635 -0.191752  0.319408              0.184791  0.274052  2.732648e-16  0.941938
# YPDBENOMYL500                 0.065957  0.225607 -0.301105  0.565374              0.085786  0.212258  4.438741e-54  0.999928
# YPDCAFEIN40                   0.104765  0.297378 -0.472642  0.518136              0.028439  0.127050  3.143833e-44  0.843116
# YPDCAFEIN50                   0.119301  0.297999 -0.478640  0.537109              0.024987  0.122008  5.345420e-48  0.755850
# YPDCHX05                     -0.033692  0.158884 -0.265770  0.320876              0.057377  0.159833  1.963417e-16  0.887620
# YPDCHX1                       0.003136  0.114423 -0.225092  0.308095              0.180794  0.277976  3.289135e-15  0.998023
# YPDCUSO410MM                  0.067088  0.273934 -0.541250  0.619061              0.124354  0.244784  2.156298e-67  0.949785
# YPDDMSO                      -0.028509  0.120207 -0.218070  0.166365              0.060244  0.165792  3.645070e-08  0.879718
# YPDETOH                      -0.033577  0.165127 -0.231873  0.281375              0.039597  0.134923  7.746282e-13  0.704323
# YPDFLUCONAZOLE               -0.008478  0.118415 -0.237368  0.221839              0.106687  0.191410  1.873747e-09  0.794451
# YPDFORMAMIDE4                -0.046035  0.110707 -0.210931  0.193937              0.095540  0.206858  1.022298e-07  0.895878
# YPDFORMAMIDE5                -0.027498  0.104491 -0.189569  0.207270              0.149474  0.248459  1.711043e-07  0.954230
# YPDHU                        -0.094950  0.129923 -0.215446  0.230194              0.017809  0.055416  5.825682e-09  0.300872
# YPDKCL2M                      0.005950  0.170329 -0.265307  0.283545              0.060723  0.158758  5.078104e-13  0.856463
# YPDLICL250MM                  0.040645  0.177282 -0.251764  0.409786              0.067116  0.141353  1.045325e-26  0.552439
# YPDMV                        -0.037208  0.173752 -0.270379  0.287200              0.071658  0.189436  2.472407e-13  0.777709
# YPDNACL15M                    0.043694  0.224504 -0.308195  0.256621              0.099784  0.272667  3.219247e-15  0.803231
# YPDNACL1M                    -0.003930  0.171785 -0.276321  0.270098              0.160356  0.299737  2.042466e-12  0.979839
# YPDNYSTATIN                   0.150796  0.186313 -0.286132  0.261442              0.004101  0.011595  3.054108e-13  0.032796
# YPDSDS                        0.054607  0.199274 -0.212683  0.279352              0.000137  0.000201  1.144726e-12  0.000595
# YPDSODIUMMETAARSENITE         0.131393  0.200554 -0.125306  0.493666              0.062027  0.141939  1.049859e-39  0.547814
# YPETHANOL                     0.032142  0.167553 -0.341599  0.371457              0.212951  0.286544  6.982756e-22  0.943015
# YPGALACTOSE                   0.016493  0.158776 -0.289950  0.296440              0.135088  0.258856  3.822645e-14  0.981352
# YPGLYCEROL                   -0.009019  0.159300 -0.238632  0.268076              0.068320  0.163824  9.518919e-12  0.645228
# YPRIBOSE                      0.065532  0.258356 -0.173300  0.305547              0.000096  0.000184  5.676133e-15  0.000372
# YPSORBITOL                    0.054920  0.130318 -0.159921  0.279801              0.116672  0.234143  1.049938e-12  0.965073
# YPXYLOSE                      0.101082  0.205689 -0.101295  0.304233              0.058603  0.109813  7.506014e-15  0.223129


# Aggregate the results from "shap_corr_summary_comparisons" for all env
excel_files = glob("Scripts/Data_Vis/_cnv_*_shap_corrwith_cnv_and_fitness.xlsx")
df_list = []
for file in excel_files:
	env = re.search("YP[A-Z0-9]+(?!cnv)", file).group()
	df = pd.read_excel(file, sheet_name="shap_corr_summary_comparisons", index_col=0)
	# pivot into long format
	df_long = df.stack().reset_index()
	df_long.columns = ["var1", "var2", "correlation"]
	df_long["env"] = env
	df_list.append(df_long)


# combine all environments
all_env_corr = pd.concat(df_list, ignore_index=True)
# pivot to compare envs side-by-side
all_env_corr = all_env_corr.pivot_table(
	index=["var1", "var2"], columns="env", values="correlation")
all_env_corr.to_csv("Scripts/Data_Vis/_cnv_all_envs_shap_corr_summary_comparisons.csv", index=True)


# Aggregate the r_x_vs_y results from "shap_corr_summary" for all env (S14 Table)
# Also fit a trendline to see how the r_cnv_vs_fitness association changes with feature importance
df_list = []
line_fits = []
for file in excel_files:
	env = re.search("YP[A-Z0-9]+(?!cnv)", file).group()
	df = pd.read_excel(file, sheet_name="shap_corr_summary", index_col=0)
	# Fit a line to compare r_cnv_vs_fitness as median_abs_shap_bins increases
	medians = df.groupby("median_abs_shap_bins")["r_cnv_vs_fitness"].median() # boxplot medians
	x = np.arange(1, len(medians) + 1)  # boxplot positions are 1-based
	X = sm.add_constant(x)
	model = sm.OLS(medians.values, X).fit()
	b = model.params[0]
	m = model.params[1]
	r2_adj = model.rsquared_adj
	pval = model.pvalues[1] # slope pvalue
	y_pred = model.predict(X)
	line_fits.append([env, m, b, r2_adj, pval])
	if env in target_envs:
		# plot r_cnv_vs_fitness correlations based on bins of median abs shap
		ax = df.boxplot(column="r_cnv_vs_fitness", by="median_abs_shap_bins", fontsize=7)
		ax.set_ylabel("pearson r (CNV value vs fitness)")
		ax.set_xlabel("median absolute SHAP value bin")
		# plot the linear trendline
		ax.plot(x, y_pred, color="red", linestyle="--", linewidth=1) # plot
		eq_text = f"y = {m:.3f}x + {b:.3f}"
		ax.text(0.05, 0.95, eq_text, transform=ax.transAxes, fontsize=8, verticalalignment='top')
		plt.xticks(rotation=45, ha="right")
		plt.title("")
		plt.suptitle("")
		plt.tight_layout()
		plt.savefig(f"Scripts/Data_Vis/_cnv_{env}_r_cnv_vs_fitness_by_shap_bin.pdf")
		plt.close("all")
		gc.collect()
	df_sub = df[["r_shap_vs_cnv", "r_cnv_vs_fitness"]].copy()
	# preserve feature identity (index → column)
	df_sub = df_sub.reset_index().rename(columns={"index": "feature"})
	# pivot into long format
	df_long = df_sub.melt(
		id_vars="feature", value_vars=["r_shap_vs_cnv", "r_cnv_vs_fitness"], var_name="metric", value_name="correlation")
	df_long["env"] = env
	df_list.append(df_long)

# combine all envs
line_fits = pd.DataFrame(line_fits, columns=["environment", "coef", "intercept", "adjusted_r2", "pvalue"])
all_env_corr = pd.concat(df_list, ignore_index=True)
wide_df = all_env_corr.pivot_table(
	index=["feature", "metric"], columns="env", values="correlation")

# save
with pd.ExcelWriter("Scripts/Data_Vis/_cnv_all_envs_shap_corr_summary.xlsx", engine="openpyxl") as writer:
	all_env_corr.to_excel(writer, sheet_name="long_format", index=False)
	wide_df.to_excel(writer, sheet_name="wide_format")
	wide_df.groupby("metric").describe().stack().T.to_excel(writer, sheet_name="mean_r_across_features")
	line_fits.to_excel(writer, sheet_name="linear_fit_r_cnv_vs_fitness")
