"""
1. Estimate the amount of genetic variation
	- bi-allelic SNP allele and genotype frequencies
	- CNV frequency per gene
	- PAV frequency per gene
	- Number of variants per gene
	
2. What clades do the diploid isolates belong to?
	- What clades are represented among the training and test sets?

3. How much genetic variation do the benchmark genes have?
	- Compare the amount of genetic variation to that of the important non-benchmark genes
	- Lack of genetic variation can come in two forms: Are the benchmark genes
	present in most or less isolates than the important non-benchmark genes?
"""

import os, gc
import datatable as dt
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# 1. Estimate the amount of genetic variation ----------------------------------
# calculate the genotype & allele frequencies from the SNPs -------
snp = dt.fread("Scripts/Data_Vis/S12_File_fixed_dtypes.csv").to_pandas() # S12 File (0,1,2 encoding)
snp.set_index("ID", inplace=True)
snp.drop(index="S288C", inplace=True)
snp = snp.astype("int")
snp2 = dt.fread("Scripts/Data_Vis/S2_File.csv").to_pandas()
snp2.set_index("ID", inplace=True)
snp2 = snp2.astype("int")
assert sum(snp.columns==snp2.columns) == 118382 # pass!

# read in test data (to separate calculations for training & test sets)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
snp_test = snp.loc[snp.index.isin(test[0]), :]
snp2_test = snp2.loc[snp2.index.isin(test[0]), :]
snp_train = snp.loc[~snp.index.isin(test[0]), :]
snp2_train = snp2.loc[~snp2.index.isin(test[0]), :]

def calc_geno_and_allele_freq(snp, snp2):
	""" Calculate genotype and allele frequencies
	Args:
		snp (pd.DataFrame): SNPs as columns, isolates as rows; 0, 1, 2 encoding
		snp2 (pd.DataFrame): SNPs as columns, isolates as rows; -1, 0, 1, encoding
	"""
	genotype_freq = {}
	# calculate genotype frequencies
	for i, row in tqdm(snp.T.iterrows(), total=snp.shape[1]):
		freqs = row.value_counts(normalize=True)
		if len(freqs) < 3:
			# add the frequencies of the missing genotypes
			for geno in [0, 1, 2]:
				if geno not in freqs.index:
					freqs[geno] = 0.0
		freqs.rename(index={0: "homozygous_ref", 1: "heterozygous_ref_alt", 2:"homozygous_alt"}, inplace=True)
		genotype_freq[i] = freqs.to_dict()
	
	for i, row in tqdm(snp2.T.iterrows(), total=snp2.shape[1]):
		freqs = row.value_counts(normalize=True)
		# add the frequencies of the missing genotypes
		for geno in [-1, 0, 1]:
			if geno not in freqs.index:
				freqs[geno] = 0.0
		freqs.rename(index={-1: "homozygous_maj", 0:"heterozygous_maj_min", 1:"homozygous_min"}, inplace=True)
		genotype_freq[i].update(freqs.to_dict())
	
	genotype_freq = pd.DataFrame(genotype_freq).T
	
	# calculate allele frequencies
	genotype_freq["ref_allele_freq"] = genotype_freq.apply(lambda x: x["homozygous_ref"] + (x["heterozygous_ref_alt"]/2), axis=1)
	genotype_freq["alt_allele_freq"] = genotype_freq.apply(lambda x: x["homozygous_alt"] + (x["heterozygous_ref_alt"]/2), axis=1)
	genotype_freq["maj_allele_freq"] = genotype_freq.apply(lambda x: x["homozygous_maj"] + (x["heterozygous_maj_min"]/2), axis=1)
	genotype_freq["min_allele_freq"] = genotype_freq.apply(lambda x: x["homozygous_min"] + (x["heterozygous_maj_min"]/2), axis=1)
	return genotype_freq


geno_freq_all = calc_geno_and_allele_freq(snp, snp2)
geno_freq_train = calc_geno_and_allele_freq(snp_train, snp2_train)
geno_freq_test = calc_geno_and_allele_freq(snp_test, snp2_test)

geno_freq_all = geno_freq_all.add_suffix("_all")
geno_freq_train = geno_freq_train.add_suffix("_train")
geno_freq_test = geno_freq_test.add_suffix("_test")

genotype_freq = pd.concat([geno_freq_all, geno_freq_train, geno_freq_test], ignore_index=False, axis=1)
genotype_freq.to_csv("Scripts/Data_Vis/reviewer_2_analysis/_genotype_allele_freq_snp_data.csv")

'''genotype_freq descriptive statistics before I converted the snp and snp2 to "int" and put in missing genotypes for each SNP:
genotype_freq.describe().T
#                                count      mean       std       min       25%       50%       75%       max
# homozygous_ref_all          118382.0  0.678533  0.266926  0.010667  0.494667  0.774667  0.910667  0.989333
# homozygous_alt_all          118330.0  0.279291  0.254715  0.001333  0.066667  0.173333  0.426667  0.984000
# heterozygous_ref_alt_all    115876.0  0.043213  0.038918  0.001333  0.017333  0.038667  0.060000  0.850667
# homozygous_maj_all          118382.0  0.778238  0.145080  0.113333  0.657333  0.814667  0.914667  0.989333
# homozygous_min_all          118330.0  0.179543  0.124739  0.001333  0.064000  0.141333  0.284000  0.788000
# heterozygous_maj_min_all    115876.0  0.043213  0.038918  0.001333  0.017333  0.038667  0.060000  0.850667
# ref_allele_freq_all         115876.0  0.694581  0.260356  0.013333  0.522667  0.792000  0.918667  0.988667
# alt_allele_freq_all         115824.0  0.305530  0.260361  0.011333  0.081333  0.208000  0.477333  0.986667
# maj_allele_freq_all         115876.0  0.796292  0.133633  0.198000  0.684667  0.830667  0.922667  0.988667
# min_allele_freq_all         115824.0  0.203774  0.133627  0.011333  0.078000  0.169333  0.316000  0.802000
# homozygous_ref_train        118382.0  0.677887  0.265974  0.009600  0.494400  0.774400  0.910400  0.990400
# heterozygous_ref_alt_train  115286.0  0.043925  0.039137  0.001600  0.018000  0.040000  0.060800  0.854400
# homozygous_alt_train        118328.0  0.279464  0.253606  0.001600  0.068800  0.174400  0.427200  0.985600
# homozygous_maj_train        118382.0  0.776550  0.145770  0.118400  0.654400  0.812800  0.913600  0.990400
# heterozygous_maj_min_train  115286.0  0.043925  0.039137  0.001600  0.018000  0.040000  0.060800  0.854400
# homozygous_min_train        118328.0  0.180756  0.125278  0.001600  0.065600  0.142400  0.284800  0.793600
# ref_allele_freq_train       115286.0  0.692998  0.259321  0.012000  0.520000  0.788800  0.916800  0.991200
# alt_allele_freq_train       115232.0  0.307119  0.259325  0.008800  0.083200  0.211200  0.480000  0.988000
# maj_allele_freq_train       115286.0  0.794094  0.134177  0.192800  0.681600  0.828000  0.920800  0.991200
# min_allele_freq_train       115232.0  0.205975  0.134169  0.008800  0.079200  0.172800  0.318400  0.807200
# homozygous_ref_test         118371.0  0.681826  0.272895  0.008000  0.496000  0.784000  0.912000  1.000000
# homozygous_alt_test         118120.0  0.278946  0.261427  0.008000  0.064000  0.168000  0.432000  0.984000
# heterozygous_ref_alt_test   105862.0  0.044628  0.040686  0.008000  0.016000  0.040000  0.056000  0.848000
# homozygous_maj_test         118382.0  0.786678  0.144026  0.088000  0.672000  0.824000  0.920000  1.000000
# homozygous_min_test         118109.0  0.173815  0.124206  0.008000  0.064000  0.136000  0.272000  0.800000
# heterozygous_maj_min_test   105862.0  0.044628  0.040686  0.008000  0.016000  0.040000  0.056000  0.848000
# ref_allele_freq_test        105851.0  0.677849  0.267133  0.012000  0.468000  0.764000  0.908000  0.996000
# alt_allele_freq_test        105611.0  0.322862  0.267179  0.012000  0.092000  0.240000  0.532000  0.988000
# maj_allele_freq_test        105862.0  0.792357  0.132398  0.200000  0.688000  0.816000  0.912000  0.996000
# min_allele_freq_test        105600.0  0.208033  0.132307  0.012000  0.088000  0.184000  0.316000  0.800000

genotype_freq.isna().sum() # there are NAs, hence why the SNP counts above are not all 118,382
genotype_freq[genotype_freq['homozygous_alt_all'].isna()] # look at some examples to figure out why there is an NA
snp['chromosome2_7725'].value_counts(normalize=True)
# chromosome2_7725
# False    0.925333 # some SNPs don't have all 3 genotypes represented, datatable converted these to True/False, convert datatypes to 'int'
# True     0.074667
# Name: proportion, dtype: float64
genotype_freq[genotype_freq['homozygous_min_all'].isna()]
snp2['chromosome2_7725'].value_counts(normalize=True)
# chromosome2_7725
# -1    0.925333 # I need to add a for loop to add the missing genotypes
#  0    0.074667
# Name: proportion, dtype: float64
'''

'''genotype_freq descriptive statistics after solving the issues described in the docstring above:
genotype_freq.describe().T
#                                count      mean       std       min       25%       50%       75%       max
# homozygous_ref_all          118382.0  0.678533  0.266926  0.010667  0.494667  0.774667  0.910667  0.989333
# homozygous_alt_all          118382.0  0.279169  0.254727  0.000000  0.066667  0.173333  0.426667  0.984000
# heterozygous_ref_alt_all    118382.0  0.042299  0.039003  0.000000  0.016000  0.037333  0.058667  0.850667
# homozygous_maj_all          118382.0  0.778238  0.145080  0.113333  0.657333  0.814667  0.914667  0.989333
# homozygous_min_all          118382.0  0.179464  0.124768  0.000000  0.064000  0.141333  0.282667  0.788000
# heterozygous_maj_min_all    118382.0  0.042299  0.039003  0.000000  0.016000  0.037333  0.058667  0.850667
# ref_allele_freq_all         118382.0  0.699682  0.260168  0.013333  0.534667  0.800667  0.922000  0.989333
# alt_allele_freq_all         118382.0  0.300318  0.260168  0.010667  0.078000  0.199333  0.465333  0.986667
# maj_allele_freq_all         118382.0  0.799387  0.133893  0.198000  0.686667  0.836000  0.925333  0.989333
# min_allele_freq_all         118382.0  0.200613  0.133893  0.010667  0.074667  0.164000  0.313333  0.802000
# homozygous_ref_train        118382.0  0.677887  0.265974  0.009600  0.494400  0.774400  0.910400  0.990400
# heterozygous_ref_alt_train  118382.0  0.042777  0.039253  0.000000  0.016000  0.038400  0.060800  0.854400
# homozygous_alt_train        118382.0  0.279337  0.253619  0.000000  0.067200  0.174400  0.425600  0.985600
# homozygous_maj_train        118382.0  0.776550  0.145770  0.118400  0.654400  0.812800  0.913600  0.990400
# heterozygous_maj_min_train  118382.0  0.042777  0.039253  0.000000  0.016000  0.038400  0.060800  0.854400
# homozygous_min_train        118382.0  0.180674  0.125309  0.000000  0.065600  0.142400  0.284800  0.793600
# ref_allele_freq_train       118382.0  0.699275  0.259128  0.012000  0.534400  0.799200  0.921600  0.991200
# alt_allele_freq_train       118382.0  0.300725  0.259128  0.008800  0.078400  0.200800  0.465600  0.988000
# maj_allele_freq_train       118382.0  0.797938  0.134501  0.192800  0.684800  0.834400  0.924800  0.991200
# min_allele_freq_train       118382.0  0.202062  0.134501  0.008800  0.075200  0.165600  0.315200  0.807200
# homozygous_ref_test         118382.0  0.681763  0.272961  0.000000  0.496000  0.784000  0.912000  1.000000
# homozygous_alt_test         118382.0  0.278329  0.261467  0.000000  0.064000  0.168000  0.424000  0.984000
# heterozygous_ref_alt_test   118382.0  0.039908  0.040849  0.000000  0.016000  0.032000  0.056000  0.848000
# homozygous_maj_test         118382.0  0.786678  0.144026  0.088000  0.672000  0.824000  0.920000  1.000000
# homozygous_min_test         118382.0  0.173414  0.124342  0.000000  0.064000  0.136000  0.272000  0.800000
# heterozygous_maj_min_test   118382.0  0.039908  0.040849  0.000000  0.016000  0.032000  0.056000  0.848000
# ref_allele_freq_test        118382.0  0.701717  0.266494  0.012000  0.536000  0.808000  0.924000  1.000000
# alt_allele_freq_test        118382.0  0.298283  0.266494  0.000000  0.076000  0.192000  0.464000  0.988000
# maj_allele_freq_test        118382.0  0.806632  0.132986  0.200000  0.696000  0.844000  0.928000  1.000000
# min_allele_freq_test        118382.0  0.193368  0.132986  0.000000  0.072000  0.156000  0.304000  0.800000

genotype_freq.var()
# homozygous_ref_all            0.071250
# homozygous_alt_all            0.064886
# heterozygous_ref_alt_all      0.001521
# homozygous_maj_all            0.021048
# homozygous_min_all            0.015567
# heterozygous_maj_min_all      0.001521
# ref_allele_freq_all           0.067687
# alt_allele_freq_all           0.067687
# maj_allele_freq_all           0.017927
# min_allele_freq_all           0.017927
# homozygous_ref_train          0.070742
# heterozygous_ref_alt_train    0.001541
# homozygous_alt_train          0.064322
# homozygous_maj_train          0.021249
# heterozygous_maj_min_train    0.001541
# homozygous_min_train          0.015702
# ref_allele_freq_train         0.067147
# alt_allele_freq_train         0.067147
# maj_allele_freq_train         0.018090
# min_allele_freq_train         0.018090
# homozygous_ref_test           0.074508
# homozygous_alt_test           0.068365
# heterozygous_ref_alt_test     0.001669
# homozygous_maj_test           0.020744
# homozygous_min_test           0.015461
# heterozygous_maj_min_test     0.001669
# ref_allele_freq_test          0.071019
# alt_allele_freq_test          0.071019
# maj_allele_freq_test          0.017685
# min_allele_freq_test          0.017685
'''


# plot genotype and allele frequencies
def plot_geno_and_allele_freq(genotype_freq, suffix, save_name):
	'''
	Args:
		genotype_freq (pd.DataFrame): object returned by calc_geno_and_allele_freq()
		suffix (str): suffix added to genotype_freq column names
		save_name (str): path and file name to save plots as
	'''
	fig, ax = plt.subplots(2, 2, figsize=(7.5, 7.5), sharey=True)
	ax[0, 0].hist(genotype_freq[[f"homozygous_ref_{suffix}", f"heterozygous_ref_alt_{suffix}", f"homozygous_alt_{suffix}"]],
				bins=100, alpha=0.5, rwidth=1, label=["Homozygous reference", "Heterozygous", "Homozygous alternative"])
	ax[0, 1].hist(genotype_freq[[f"ref_allele_freq_{suffix}", f"alt_allele_freq_{suffix}"]],
			   bins=100, alpha=0.5, rwidth=1, label=["Reference", "Alternative"])
	ax[1, 0].hist(genotype_freq[[f"homozygous_maj_{suffix}", f"heterozygous_maj_min_{suffix}", f"homozygous_min_{suffix}"]],
				bins=100, alpha=0.5, rwidth=1, label=["Homozygous major", "Heterozygous", "Homozygous minor"])
	ax[1, 1].hist(genotype_freq[[f"maj_allele_freq_{suffix}", f"min_allele_freq_{suffix}"]],
			   bins=100, alpha=0.5, rwidth=1, label=["Major", "Minor"])
	for i in range(4):
		r_idx = i // 2
		c_idx = i % 2
		ax[r_idx, c_idx].set_ylabel("Frequency (number of SNPs)")
		ax[r_idx, c_idx].minorticks_on()
		ax[r_idx, c_idx].legend(prop={"size": 7}, loc="upper center")
		if i < 2:
			ax[i, 0].set_title("Distribution of genotype frequencies")
			ax[i, 1].set_title("Distribution of allele frequencies")
			ax[i, 0].set_xlabel("Genotype frequency (number of isolates with the genotype / total number of isolates)")
			ax[i, 1].set_xlabel("Allele frequency")
	fig.tight_layout()
	plt.savefig(save_name)
	plt.close("all")
	gc.collect()


plot_geno_and_allele_freq(geno_freq_all, "all", "Scripts/Data_Vis/reviewer_2_analysis/_genotype_allele_freq_snp_distributions_all.pdf")
plot_geno_and_allele_freq(geno_freq_train, "train", "Scripts/Data_Vis/reviewer_2_analysis/_genotype_allele_freq_snp_distributions_train.pdf")
plot_geno_and_allele_freq(geno_freq_test, "test", "Scripts/Data_Vis/reviewer_2_analysis/_genotype_allele_freq_snp_distributions_test.pdf")


# count the number of SNPs per gene -------------------------------
# Note: use "gene" region, which bundles up gene's transcripts & regulatory elements)
snp_map = pd.read_excel("Scripts/Data_Vis/S8_File.xlsx", engine="openpyxl") # SNP to gene mapping
snp_map["gene"] = snp_map["gene"].str.split(",") # split up the genes into a list
snp_map = snp_map.explode(column="gene", ignore_index=True) # each gene gets its own row
snp_counts = snp_map.groupby("gene")["snp"].nunique().sort_values(ascending=False)
snp_counts.to_csv("Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_snps.csv")
snp_counts.drop("intergenic").hist(bins=16, figsize=(5, 4))
plt.xlim(1, 170)
plt.xlabel("Number of SNPs per gene")
plt.ylabel("Frequency (number of genes)")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_snps_distribution.pdf")
plt.close("all")

snp_counts = pd.read_csv("Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_snps.csv", index_col=0)
snp_counts.drop("intergenic").describe()
#                snp
# count  6285.000000
# mean     11.973906
# std      11.571628
# min       1.000000
# 25%       4.000000
# 50%       9.000000
# 75%      16.000000
# max     168.000000

snp_counts.drop("intergenic").var() # 133.902565

# calculate CNV frequency per gene --------------------------------
# Number of ORFs per gene
orf_map = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl") # ORF to gene mapping
orf_map.groupby("gene")["orf"].nunique().sort_values(ascending=False).to_csv("Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_orfs.csv")
# note: 29 genes have 2 ORFs.

# read in the CNV data
cnv = pd.read_excel("Scripts/Data_Vis/S6_File.xlsx", engine="openpyxl") # cnv data
cnv.set_index("ID", inplace=True)

def fix_orf_ids(orf_df):
	"""Fix the ORF feature names in the SHAP files to match the actual ORF
	identifiers"""
	orf_df.columns = orf_df.columns.str.replace("^X", "", regex=True).str.replace("\.", "-", regex=True)
	return orf_df

cnv = fix_orf_ids(cnv)
cnv = cnv.T
cnv["gene"] = cnv.index
cnv["gene"] = cnv["gene"].map(lambda x: orf_map.set_index("orf")["gene"].get(x, x))

# how to determine CNV variation for genes with multiple ORFs?
cnv.loc[cnv.gene == "YIL068C", :].apply(lambda x: x.value_counts().to_dict(), axis=1)
# note: one orf has more CNV variation than the other

# summary statistics per ORF:
cnv_summary = cnv.iloc[:,:-1].T.describe().T
cnv_summary["gene"] = cnv_summary.index.map(lambda x: orf_map.set_index("orf")["gene"].get(x, x))
cnv_summary.insert(8, "variance", cnv.iloc[:,:-1].T.var())
cnv_summary.index.name = "orf"
cnv_summary["mean"].describe() # mean( the mean CNV value in an ORF ) : mean of means across ORFs
# count    7708.000000
# mean        0.882930
# std         2.332259
# min         0.000000
# 25%         0.836000
# 50%         1.002000
# 75%         1.007333
# max        99.172667

# get the number of isolates with each CNV value for each ORF
cnv_counts = cnv.iloc[:,:-1].apply(lambda x: x.value_counts().to_dict(), axis=1)
cnv_counts = pd.DataFrame(cnv_counts.tolist(), index=cnv.index)
cnv_counts = cnv_counts[sorted(cnv_counts.columns.to_list())]
cnv_counts.fillna(0.0, inplace=True)
cnv_counts.index.name = "orf"
cnv_counts_summary = cnv_counts.describe().T
cnv_counts_summary.index.name = "cnv_value"

with pd.ExcelWriter("Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_cnv_isolates.xlsx") as writer:
	cnv_summary.to_excel(writer, sheet_name="summary_stats_per_orf")
	cnv_counts.to_excel(writer, sheet_name="num_isolates_per_cnv_value")
	cnv_counts_summary.to_excel(writer, sheet_name="summary_stats_cnv_counts")

# plot CNV summary stats
cnv_summary = pd.read_excel(
	"Scripts/Data_Vis/reviewer_2_analysis/_count_per_gene_cnv_isolates.xlsx",
	sheet_name="summary_stats_per_orf", engine="openpyxl")

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
ax[0].hist(cnv_summary["variance"], bins=50, color="lightblue", label="CNV variance across isolates")
ax[1].hist(cnv_summary["mean"], bins=50, color="green", label="mean CNV value across isolates")
for i in [0, 1]:
	ax[i].set_yscale("log")
	ax[i].set_ylabel("Number of ORFs")
	ax[i].set_box_aspect(1)
	ax[i].minorticks_on()

ax[0].legend(loc="upper left")
ax[1].legend(loc="upper right")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/reviewer_2_analysis/_cnv_isolates_var_mean.pdf")

plt.figure(figsize=(8, 3))
plt.errorbar(x=cnv_counts_summary.index, y=cnv_counts_summary["mean"],
	yerr=cnv_counts_summary["std"], fmt="o", color="blue", ecolor="lightblue",
	elinewidth=0.6, capsize=1, ms=1) # cnv_value = 0 & 0.5 have the most number of orfs
plt.grid(True, linestyle="--", alpha=0.5)
plt.minorticks_on()
plt.yscale("log")
plt.ylabel("Number of ORFs")
plt.xlabel("ORF copy number")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/reviewer_2_analysis/_cnv_values_mean_orf_num.pdf")
plt.close("all")
gc.collect()

# calculate PAV percent presence ----------------------------------
orf_map = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl") # ORF to gene mapping
pav = pd.read_excel("Scripts/Data_Vis/S5_File.xlsx", engine="openpyxl") # pav data
pav.set_index("ID", inplace=True)
pav = fix_orf_ids(pav)
pav_presence = pd.DataFrame(pav.sum(axis=0)/750)
pav_presence.columns = ["percent_presence"]
pav_presence["gene"] = pav_presence.index.map(lambda x: orf_map.set_index("orf")["gene"].get(x, x))
pav_presence.to_csv("Scripts/Data_Vis/reviewer_2_analysis/_pav_percent_presence.csv")

fig, ax = plt.subplots(figsize=(3, 3))
ax.hist(pav_presence["percent_presence"], bins=100, color="lightblue", label="ORF % presence")
ax.set_yscale("log")
ax.minorticks_on()
ax.legend()
ax.set_box_aspect(1)
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/reviewer_2_analysis/_pav_percent_presence.pdf")
plt.close("all")

pav_presence = pd.read_csv("Scripts/Data_Vis/reviewer_2_analysis/_pav_percent_presence.csv")
pav_presence["percent_presence"].describe()
# count    7708.000000
# mean        0.780006
# std         0.401554
# min         0.000000
# 25%         0.966667
# 50%         1.000000
# 75%         1.000000
# max         1.000000
pav_presence["percent_presence"].var() # 0.16

# How does the percent presence of benchmark genes compare to that of important non-benchmark genes?
target_envs = ["YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", "YPDCUSO410MM", "YPDSODIUMMETAARSENITE"]

# read in the feature tables from Fig4A to get the benchmark and important non-benchmark feature sets
from glob import glob
data_dir = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP/"
snp_shap_files = glob(os.path.join(data_dir, "SNP/fs/SHAP_values_sorted_Y*.txt"))
pav_shap_files = glob(os.path.join(data_dir, "PAV/fs/SHAP_values_sorted_Y*.txt"))
cnv_shap_files = glob(os.path.join(data_dir, "CNV/fs/SHAP_values_sorted_Y*.txt"))
snp_shap_file = [f for f in snp_shap_files if any(env in f for env in target_envs)]
pav_shap_file = [f for f in pav_shap_files if any(env in f for env in target_envs)]
cnv_shap_file = [f for f in cnv_shap_files if any(env in f for env in target_envs)]

# subset the PAV matrix to only include the training isolates
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)
pav_train = pav[~pav.index.isin(test[0])]

# subset
results = []
for env, bench_env in [("YPDBENOMYL500", "Benomyl"), ("YPDCAFEIN40", "Caffeine"),
					   ("YPDCAFEIN50", "Caffeine"), ("YPDCUSO410MM", "CuSO4"),
					   ("YPDSODIUMMETAARSENITE", "Sodium_meta-arsenite")]:
	for vtype in ["pav", "cnv"]:
		# Load Figure 4A feature lists
		d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
		orf = pd.read_csv(f"{d}/Features_important_non_bench_plus_bench_genes_{env}_{vtype}.txt", header=None)
		orf["clean_id"] = fix_orf_ids(orf[0])
		
		# subset benchmark gene features and non-benchmark features
		feat1_bench = orf_map[(orf_map.orf.isin(orf["clean_id"])) & (orf_map[bench_env] == 1)]
		allfeat_bench = orf_map.loc[orf_map.gene.isin(feat1_bench.gene.unique()),:]
		feat1_non_bench = orf_map[(orf_map.orf.isin(orf["clean_id"])) & (orf_map[bench_env] == 0)]
		
		# subset the genotype matrices
		feat1_bench_df = pav_train.loc[:, fix_orf_ids(feat1_bench.orf, clean=False)]
		allfeat_bench_df = pav_train.loc[:, fix_orf_ids(allfeat_bench.orf, clean=False)]
		feat1_non_bench_df = pav_train.loc[:, fix_orf_ids(feat1_non_bench.orf, clean=False)]
		
		# calculate percent presence for each feature set
		feat1_bench_percent_presence = feat1_bench_df.sum(axis=0)/625
		allfeat_bench_percent_presence = allfeat_bench_df.sum(axis=0)/625
		feat1_non_bench_percent_presence = feat1_non_bench_df.sum(axis=0)/625
		
		# statistically estimate how different the mean percent presence values are between the feature sets using a Welch's t-test
		from scipy.stats import ttest_ind
		t_stat_g_feat1, p_val_g_feat1 = ttest_ind(
			feat1_bench_percent_presence, feat1_non_bench_percent_presence,
			equal_var=False, alternative="greater")
		t_stat_l_feat1, p_val_l_feat1 = ttest_ind(
			feat1_bench_percent_presence, feat1_non_bench_percent_presence,
			equal_var=False, alternative="less")
		t_stat_g_all, p_val_g_all = ttest_ind(
			allfeat_bench_percent_presence, feat1_non_bench_percent_presence,
			equal_var=False, alternative="greater")
		t_stat_l_all, p_val_l_all = ttest_ind(
			allfeat_bench_percent_presence, feat1_non_bench_percent_presence,
			equal_var=False, alternative="less")
		
		results.append([env, vtype, feat1_bench_percent_presence.mean(),
					feat1_bench_percent_presence.std(),
					feat1_bench_percent_presence.var(),
					feat1_non_bench_percent_presence.mean(),
					feat1_non_bench_percent_presence.std(),
					feat1_non_bench_percent_presence.var(),
					allfeat_bench_percent_presence.mean(),
					allfeat_bench_percent_presence.std(),
					allfeat_bench_percent_presence.var(),
					t_stat_g_feat1, p_val_g_feat1, t_stat_l_feat1, p_val_l_feat1,
					t_stat_g_all, p_val_g_all, t_stat_l_all, p_val_l_all])


results_df = pd.DataFrame(results, columns=[
	"environment", "variant_type",
	"feat1_bench_mean_percent_presence", "feat1_bench_std_percent_presence",
	"feat1_bench_var_percent_presence", "feat1_non_bench_mean_percent_presence",
	"feat1_non_bench_std_percent_presence", "feat1_non_bench_var_percent_presence",
	"allfeat_bench_mean_percent_presence", "allfeat_bench_std_percent_presence",
	"allfeat_bench_var_percent_presence", "t_stat_g_feat1", "p_val_g_feat1",
	"t_stat_l_feat1", "p_val_l_feat1", "t_stat_g_all", "p_val_g_all", "t_stat_l_all",
	"p_val_l_all"])
results_df.to_csv("Scripts/Data_Vis/reviewer_2_analysis/_benchmark_vs_non_bench_percent_presence_stats.csv", index=False)


# 2. Diploid isolate clade representation --------------------------------------
pheno = pd.read_excel("S1_File.xlsx", engine="openpyxl")
meta = pd.read_excel("../../Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls",
	sheet_name="Table S1", skiprows=3)
test = pd.read_csv("../../Data/Peter_2018/Test.txt", header=None)

meta_diploids = meta.loc[meta["Standardized name"].isin(pheno["ID"])]
meta_diploids.shape[0] == 750 # True
meta_diploids["train_test_assignment"] = (meta_diploids["Standardized name"].isin(test[0]).map({True: "Test", False: "Train"}))
meta_diploids["train_test_assignment"].value_counts() # looks good

# What clades are represented among the training and test sets? 26 in Test, 36 in Train
pd.set_option("display.max_rows", 62)
meta_diploids.groupby("train_test_assignment")["Clades"].value_counts()
# train_test_assignment  Clades                       
# Test                   1. Wine/European                  43
#                        M3. Mosaic region 3               10
#                        8. Mixed origin                    8
#                        25. Sake                           7
#                        5. French dairy                    5
#                        1. Wine/European (subclade 3)      4
#                        1. Wine/European (subclade 4)      4
#                        M2. Mosaic region 2                4
#                        10. French Guiana human            3
#                        23. North American oak             3
#                        26. Asian fermentation             3
#                        3. Brazilian bioethanol            3
#                        9. Mexican agave                   3
#                        1. Wine/European (subclade 1)      2
#                        13. African palm wine              2
#                        2. Alpechin                        2
#                        22. Far East Russian               2
#                        M3. Mosaic region 3                2
#                        1. Wine/European (subclade 2)      1
#                        11. Ale beer                       1
#                        12. West African cocoa             1
#                        18. Far East Asia                  1
#                        21. Ecuadorean                     1
#                        24. Asian islands                  1
#                        4. Mediterranean oak               1
#                        M1. Mosaic region 1                1
# Train                  1. Wine/European                 184
#                        M3. Mosaic region 3               50
#                        25. Sake                          34
#                        8. Mixed origin                   29
#                        10. French Guiana human           27
#                        26. Asian fermentation            24
#                        3. Brazilian bioethanol           24
#                        1. Wine/European (subclade 4)     23
#                        5. French dairy                   23
#                        1. Wine/European (subclade 3)     20
#                        1. Wine/European (subclade 1)     14
#                        12. West African cocoa            12
#                        13. African palm wine             12
#                        2. Alpechin                       10
#                        7. Mosaic beer                    10
#                        M2. Mosaic region 2               10
#                        M1. Mosaic region 1                9
#                        1. Wine/European (subclade 2)      8
#                        21. Ecuadorean                     8
#                        23. North American oak             8
#                        M3. Mosaic region 3                8
#                        18. Far East Asia                  7
#                        4. Mediterranean oak               7
#                        24. Asian islands                  6
#                        19. Malaysian                      5
#                        6. African beer                    5
#                        9. Mexican agave                   4
#                        17. Taiwanese                      3
#                        14. CHNIII                         2
#                        15. CHNII                          2
#                        20. CHN V                          2
#                        22. Far East Russian               2
#                        11. Ale beer                       1
#                        16. CHNI                           1
#                        M1. Mosaic region 1                1
#                        M2. Mosaic region 2                1
# Name: count, dtype: int64

meta_diploids.groupby("train_test_assignment")["Ecological origins"].value_counts()
# train_test_assignment  Ecological origins
# Test                   Wine                   37
#                        Tree                   15
#                        Human, clinical        10
#                        Nature                  7
#                        Sake                    7
#                        Fruit                   6
#                        Soil                    6
#                        Bakery                  4
#                        Dairy                   4
#                        Bioethanol              3
#                        Distillery              3
#                        Human                   3
#                        Insect                  3
#                        Palm wine               3
#                        Water                   3
#                        Cider                   2
#                        Fermentation            2
#                        Industrial              2
#                        Probiotic               2
#                        Beer                    1
#                        Flower                  1
#                        Unknown                 1
# Train                  Wine                  167
#                        Human, clinical        66
#                        Sake                   33
#                        Tree                   33
#                        Nature                 31
#                        Fermentation           27
#                        Fruit                  27
#                        Human                  25
#                        Soil                   25
#                        Bioethanol             22
#                        Bakery                 21
#                        Dairy                  20
#                        Distillery             18
#                        Unknown                17
#                        Beer                   15
#                        Insect                 15
#                        Cider                  13
#                        Flower                 13
#                        Water                  13
#                        Industrial             11
#                        Palm wine              11
#                        Lab strain              2
# Name: count, dtype: int64


# 3. Amount of genetic variation of benchmark genes ----------------------------
# Genetic variant matrices
# snp_geno = dt.fread("Scripts/Data_Vis/S12_File_fixed_dtypes.csv").to_pandas() # S12 File (0, 1, 2 encoding)
snp_geno = dt.fread("Scripts/Data_Vis/S2_File.csv").to_pandas()
pav_geno = pd.read_excel("Scripts/Data_Vis/S5_File.xlsx", engine="openpyxl")
cnv_geno = pd.read_excel("Scripts/Data_Vis/S6_File.xlsx", engine="openpyxl")
snp_geno.set_index("ID", inplace=True)
# snp_geno.drop(index="S288C", inplace=True)
pav_geno.set_index("ID", inplace=True)
cnv_geno.set_index("ID", inplace=True)
pav_geno = pav_geno.loc[snp_geno.index, :] # ensure isolates are in the same order
cnv_geno = cnv_geno.loc[snp_geno.index, :]

# SNP and ORF gene mappings
map_snps = pd.read_excel("Scripts/Data_Vis/S8_File.xlsx", engine="openpyxl")
map_orfs = pd.read_excel("Scripts/Data_Vis/S9_File.xlsx", engine="openpyxl")
map_snps.set_index("snp", inplace=True)
map_orfs.set_index("orf", inplace=True)

# Test set (to subset training data only)
test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)

# function to fix ORF IDs
def fix_orf_ids(pd_series, clean=True):
	if clean:
		return pd_series.str.replace("^X", "", regex=True).str.replace("\.", "-", regex=True)
	else:
		return "X" + pd_series.str.replace("-", ".", regex=True)


# function to estimate genotype and allele frequencies
def calc_geno_freqs(df, minor=True):
	genotype_freq = {}
	for i, row in df.T.iterrows():
		genotype_freq[i] = row.value_counts(normalize=True).to_dict()
	genotype_freq = pd.DataFrame(genotype_freq).T
	# calculate allele frequencies
	if minor: # -1, 0, 1 encoding
		genotype_freq["minor_allele_freq"] = genotype_freq.apply(lambda x: x[1] + x[0]/2, axis=1)
		genotype_freq["major_allele_freq"] = genotype_freq.apply(lambda x: x[-1] + x[0]/2, axis=1)
	else:
		genotype_freq["ref_allele_freq"] = genotype_freq.apply(lambda x: x[0] + (x[1]/2), axis=1)
		genotype_freq["alt_allele_freq"] = genotype_freq.apply(lambda x: x[2] + (x[1]/2), axis=1)
	return genotype_freq


# Estimate the amount of genetic variation among the benchmark gene features
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon, mannwhitneyu
results = []
for env, bench_env in [("YPDBENOMYL500", "Benomyl"), ("YPDCAFEIN40", "Caffeine"),
					   ("YPDCAFEIN50", "Caffeine"), ("YPDCUSO410MM", "CuSO4"),
					   ("YPDSODIUMMETAARSENITE", "Sodium_meta-arsenite")]:
	for vtype in ["snp", "pav", "cnv"]:
		# Load Figure 4A feature lists
		d = "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SHAP_Interaction/RF"
		if vtype == "snp":
			snp = pd.read_csv(f"{d}/Features_bench_plus_important_non_bench_{env}_snp.txt", header=None)
		else:
			orf = pd.read_csv(f"{d}/Features_important_non_bench_plus_bench_genes_{env}_{vtype}.txt", header=None)
			orf["clean_id"] = fix_orf_ids(orf[0])
		
		# subset benchmark gene features and non-benchmark features
		if vtype == "snp":
			feat1_bench = map_snps[(map_snps.index.isin(snp[0])) & (map_snps[bench_env] == 1)] # 1 feature per gene; used to train Fig 4A models
			allfeat_bench = map_snps.loc[map_snps.gene.isin(feat1_bench.gene.unique()),:] # all features per gene
			feat1_non_bench = map_snps[(map_snps.index.isin(snp[0])) & (map_snps[bench_env] == 0)] # important non-benchmark features
			# subset the genotype matrices
			feat1_bench_df = snp_geno.loc[~snp_geno.index.isin(test[0]), feat1_bench.index]
			allfeat_bench_df = snp_geno.loc[~snp_geno.index.isin(test[0]), allfeat_bench.index]
			feat1_non_bench_df = snp_geno.loc[~snp_geno.index.isin(test[0]), feat1_non_bench.index]
		else:
			feat1_bench = map_orfs[(map_orfs.index.isin(orf["clean_id"])) & (map_orfs[bench_env] == 1)]
			allfeat_bench = map_orfs.loc[map_orfs.gene.isin(feat1_bench.gene.unique()),:]
			feat1_non_bench = map_orfs[(map_orfs.index.isin(orf["clean_id"])) & (map_orfs[bench_env] == 0)]
			# subset the genotype matrices
			if vtype == "pav":
				feat1_bench_df = pav_geno.loc[~pav_geno.index.isin(test[0]), fix_orf_ids(feat1_bench.index, clean=False)]
				allfeat_bench_df = pav_geno.loc[~pav_geno.index.isin(test[0]), fix_orf_ids(allfeat_bench.index, clean=False)]
				feat1_non_bench_df = pav_geno.loc[~pav_geno.index.isin(test[0]), fix_orf_ids(feat1_non_bench.index, clean=False)]
			else:
				feat1_bench_df = cnv_geno.loc[~cnv_geno.index.isin(test[0]), fix_orf_ids(feat1_bench.index, clean=False)]
				allfeat_bench_df = cnv_geno.loc[~cnv_geno.index.isin(test[0]), fix_orf_ids(allfeat_bench.index, clean=False)]
				feat1_non_bench_df = cnv_geno.loc[~cnv_geno.index.isin(test[0]), fix_orf_ids(feat1_non_bench.index, clean=False)]
		
		# Calculate the Pearson correlation and Euclidean distance for isolate pairs
		r_feat1_bench = feat1_bench_df.T.corr(method="pearson")
		e_feat1_bench = pd.DataFrame(squareform(pdist(feat1_bench_df, metric="euclidean")),
			index=feat1_bench_df.index, columns=feat1_bench_df.index)
		r_allfeat_bench = allfeat_bench_df.T.corr(method="pearson")
		e_allfeat_bench = pd.DataFrame(squareform(pdist(allfeat_bench_df, metric="euclidean")),
			index=allfeat_bench_df.index, columns=allfeat_bench_df.index)
		r_feat1_non_bench = feat1_non_bench_df.T.corr(method="pearson")
		e_feat1_non_bench = pd.DataFrame(squareform(pdist(feat1_non_bench_df, metric="euclidean")),
			index=feat1_non_bench_df.index, columns=feat1_non_bench_df.index)
		
		# Compare the Euclidean distance distributions with Wilcoxon signed-rank test
		def melt_df(df):
			# Reshape the dataframe to long format for comparing distributions
			out = pd.DataFrame(np.triu(df, k=1), index=df.index, columns=df.columns)
			out_melt = out.reset_index().melt(id_vars="ID")
			out_melt.columns = ["Isolate1", "Isolate2", "value"]
			return out_melt.set_index(["Isolate1", "Isolate2"])
		
		r_feat1_bench_melt = melt_df(r_feat1_bench)
		r_allfeat_bench_melt = melt_df(r_allfeat_bench)
		r_feat1_non_bench_melt = melt_df(r_feat1_non_bench)
		e_feat1_bench_melt = melt_df(e_feat1_bench)
		e_allfeat_bench_melt = melt_df(e_allfeat_bench)
		e_feat1_non_bench_melt = melt_df(e_feat1_non_bench)
		
		w_r_f1_bench, w_p_r_f1_bench = wilcoxon(x=r_feat1_bench_melt - r_feat1_non_bench_melt,
			zero_method="wilcox", alternative="less", nan_policy="omit")
		w_r_af_bench, w_p_r_af_bench = wilcoxon(x=r_allfeat_bench_melt - r_feat1_non_bench_melt,
			zero_method="wilcox", alternative="less", nan_policy="omit")
		w_e_f1_bench, w_p_e_f1_bench = wilcoxon(x=e_feat1_bench_melt - e_feat1_non_bench_melt,
			zero_method="wilcox", alternative="less", nan_policy="omit")
		w_e_af_bench, w_p_e_af_bench = wilcoxon(x=e_allfeat_bench_melt - e_feat1_non_bench_melt,
			zero_method="wilcox", alternative="less", nan_policy="omit")
		
		if vtype == "snp":
			# Calculate the genotype and minor/major allele frequencies
			freq_feat1_bench = calc_geno_freqs(feat1_bench_df, minor=True)
			freq_allfeat_bench = calc_geno_freqs(allfeat_bench_df, minor=True)
			freq_feat1_non_bench = calc_geno_freqs(feat1_non_bench_df, minor=True)
			
			# Compare the allele frequency distributions with Wilcoxon signed-rank test
			u_minor_1feat, u_p_minor_1feat = mannwhitneyu(
				x=freq_feat1_bench["minor_allele_freq"],
				y=freq_feat1_non_bench["minor_allele_freq"],
				alternative="less", nan_policy="omit")
			u_minor_allfeat, u_p_minor_allfeat = mannwhitneyu(
				x=freq_allfeat_bench["minor_allele_freq"],
				y=freq_feat1_non_bench["minor_allele_freq"],
				alternative="less", nan_policy="omit")
			u_major_1feat, u_p_major_1feat = mannwhitneyu(
				x=freq_feat1_bench["major_allele_freq"],
				y=freq_feat1_non_bench["major_allele_freq"],
				alternative="less", nan_policy="omit")
			u_major_allfeat, u_p_major_allfeat = mannwhitneyu(
				x=freq_allfeat_bench["major_allele_freq"],
				y=freq_feat1_non_bench["major_allele_freq"],
				alternative="less", nan_policy="omit")
			results.append([env, bench_env, vtype, w_r_f1_bench, w_p_r_f1_bench,
				w_r_af_bench, w_p_r_af_bench, w_e_f1_bench, w_p_e_f1_bench,
				w_e_af_bench, w_p_e_af_bench, u_minor_1feat, u_p_minor_1feat,
				u_minor_allfeat, u_p_minor_allfeat, u_major_1feat, u_p_major_1feat,
				u_major_allfeat, u_p_major_allfeat, len(freq_feat1_bench),
				len(freq_feat1_non_bench), len(freq_allfeat_bench)])
			del freq_feat1_bench, freq_allfeat_bench, freq_feat1_non_bench
			del u_minor_1feat, u_p_minor_1feat, u_minor_allfeat, u_p_minor_allfeat
			del u_major_1feat, u_p_major_1feat, u_major_allfeat, u_p_major_allfeat
		else:
			results.append([env, bench_env, vtype, w_r_f1_bench, w_p_r_f1_bench,
				w_r_af_bench, w_p_r_af_bench, w_e_f1_bench, w_p_e_f1_bench,
				w_e_af_bench, w_p_e_af_bench] + [None]*11)
		
		# clear memory
		del feat1_bench, allfeat_bench, feat1_non_bench
		del feat1_bench_df, allfeat_bench_df, feat1_non_bench_df
		del r_feat1_bench, r_feat1_bench_melt, r_allfeat_bench, r_allfeat_bench_melt, r_feat1_non_bench, r_feat1_non_bench_melt
		del e_feat1_bench, e_feat1_bench_melt, e_allfeat_bench, e_allfeat_bench_melt, e_feat1_non_bench, e_feat1_non_bench_melt
		del w_r_f1_bench, w_p_r_f1_bench, w_r_af_bench, w_p_r_af_bench, w_e_f1_bench, w_p_e_f1_bench, w_e_af_bench, w_p_e_af_bench


results_df = pd.DataFrame(results, columns = ["Environment", "Benchmark_gene_set",
	"Variant_type", "w_r_1feat_bench", "w_pval_r_1feat_bench", "w_r_allfeat_bench",
	"w_pval_r_allfeat_bench", "w_e_1feat_bench", "w_pval_e_1feat_bench",
	"w_e_allfeat_bench", "w_pval_e_allfeat_bench", "u_minor_1feat",
	"u_p_minor_1feat", "u_minor_allfeat", "u_p_minor_allfeat", "u_major_1feat",
	"u_p_major_1feat", "u_major_allfeat", "u_p_major_allfeat", "Num_1feat_bench",
	"Num_1feat_non_bench", "Num_allfeat_bench"])
for col in results_df.columns[3:]:
	results_df[col] = results_df.apply(lambda r: r[col][0], axis=1)

results_df.to_csv("Scripts/Data_Vis/reviewer_2_analysis/_genetic_variance_bench_vs_non_bench.csv", index=False)







