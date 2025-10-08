#!/usr/bin/env python3
################################################################################
# Figure 1B-E
################################################################################
import os
import math
import json
import pandas as pd
import datatable as dt
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib.patches import Rectangle
from matplotlib import cm
from matplotlib.colors import Normalize, ListedColormap, rgb2hex


def plot_colortable(colors, *, ncols=4, sort_colors=True):
    """
    Plot a table of colors with their names
    source: https://matplotlib.org/stable/gallery/color/named_colors.html 
    """
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    n = len(colors)
    nrows = math.ceil(n / ncols)
    width = cell_width * ncols + 2 * margin
    height = cell_height * nrows + 2 * margin
    dpi = 72
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-margin)/height)
    ax.set_xlim(0, cell_width * ncols)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    for i, name in enumerate(colors):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
        swatch_start_x = cell_width * col
        text_pos_x = cell_width * col + swatch_width + 7
        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')
        ax.add_patch(
            Rectangle(xy=(swatch_start_x, y-9), width=swatch_width,
                      height=18, facecolor=colors[i], edgecolor='0.7')
        )
    return fig


os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Environment truncated and full names mapping
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

pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
pheno.rename(columns=mapping, inplace=True)  # full environment names
# fitness correlations between environments
pCorEnvs = pheno.corr(method="pearson")

# S1b-e. Heatmaps of isolate pair data correlations
kinship = pd.read_csv("Data/Peter_2018/kinship.csv", index_col=0)
snps = dt.fread("Data/Peter_2018/geno.csv").to_pandas()
snps.set_index("ID", inplace=True)
pavs = pd.read_csv("Data/Peter_2018/ORFs_pres_abs.csv", index_col=0)
pavs = pavs.loc[kinship.index, :]  # sort ORFs by kinship isolate order
cnvs = pd.read_csv("Data/Peter_2018/ORFs_no_NA.csv", index_col=0)
cnvs = cnvs.loc[kinship.index, :]  # sort CNVs by kinship isolate order
clades = pd.read_excel("Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls",
                       sheet_name="Table S1", skiprows=3, nrows=1011)  # isolate clades

# summary statistics for each isolate and each clade
# fitness, snp heterozygosity, pav stats, cnv stats, number of isolates
clades = clades[["Standardized name", "Clades"]]  # subset relevant columns
pheno.merge(clades, left_index=True, right_on="Standardized name").\
    groupby("Clades").describe().to_csv(
        "Scripts/Data_Vis/Section_1/pheno_clades_stats.csv")

# cluster data
kin_dendrogram = dendrogram(linkage(kinship, method="average"))
kin_row_order = kin_dendrogram["leaves"]
pheno_corr = pheno.T.corr(method="pearson")
pavs_corr = pavs.T.corr(method="pearson")
cnvs_corr = cnvs.T.corr(method="pearson")
pheno_corr.to_csv("Data/Peter_2018/fitness_correlations_isolate_pair.csv")
pavs_corr.to_csv("Data/Peter_2018/ORFs_pres_abs_correlations_isolate_pair.csv")
cnvs_corr.to_csv("Data/Peter_2018/ORFs_copy_num_correlations_isolate_pair.csv")
sum(kinship.index == pheno_corr.index)
sum(kinship.index == pavs_corr.index)
sum(kinship.index == cnvs_corr.index)
sum(kinship.columns == pheno_corr.columns)
sum(kinship.columns == pavs_corr.columns)
sum(kinship.columns == cnvs_corr.columns)

# Map clades to isolates and create a colormap
clades = clades[["Standardized name", "Clades"]]  # subset relevant columns
clades = clades.loc[clades["Standardized name"].isin(
    kinship.index)]  # diploid isolates
clades.set_index("Standardized name", inplace=True)
clades = clades.loc[kinship.index, :]  # reorder by original order
clades = clades.iloc[kin_row_order, :]  # reorder by dendrogram order
# replace NaN with Unknown
clades.loc[clades.Clades.isna(), "Clades"] = "Unknown"

cmap = cm.get_cmap('hsv', clades.Clades.nunique())  # sample hsv color palette
color_list = {clades.Clades.unique()[i]: rgb2hex(
    cmap(i)[:3]) for i in range(cmap.N)}  # assign colors to clades
row_colors = clades['Clades'].map(color_list)  # assign colors to isolates

plot_colortable(row_colors.unique())  # plot color table
plt.savefig("Scripts/Data_Vis/Section_1/clade_colors.png")
plt.close()
with open("Scripts/Data_Vis/Section_1/clade_colors.json", "w") as f:
    json.dump(color_list, f, indent=4)  # save color dictionary

# plot with clades (Figure 1C)
fig, ax = plt.subplots(2, 2)
kin_ordered = kinship.iloc[kin_row_order, kin_row_order]
a = sns.heatmap(kin_ordered, cmap="RdBu_r", center=0, square=True, ax=ax[0][0],
                xticklabels=False, yticklabels=False, vmin=kin_ordered.min().min(),
                vmax=kin_ordered.max().max(), cbar_kws={'orientation': 'horizontal'})
ax[0][0].tick_params(axis="y", which="major", pad=20,
                     length=0)  # extra padding for row colors
for i, color in enumerate(row_colors):
    ax[0][0].add_patch(plt.Rectangle(xy=(-0.05, i), width=0.05, height=1, color=color, lw=0,
                                     transform=ax[0][0].get_yaxis_transform(), clip_on=False))

# Figure 1B,D,E
pheno_corr_ordered = pheno_corr.iloc[kin_row_order, kin_row_order]
b = sns.heatmap(pheno_corr_ordered, cmap="RdBu_r", center=pheno_corr_ordered.mean().mean(),
                square=True, ax=ax[0][1], xticklabels=False, yticklabels=False,
                vmin=pheno_corr_ordered.min().min(), vmax=pheno_corr_ordered.max().max(),
                cbar_kws={'orientation': 'horizontal'})
pavs_corr[pavs_corr < 0.8] = 0.8  # set minimum value to 0.8
pav_corr_ordered = pavs_corr.iloc[kin_row_order, kin_row_order]
c = sns.heatmap(pav_corr_ordered, cmap="RdBu_r", center=0.9, square=True,
                ax=ax[1][0], xticklabels=False, yticklabels=False,
                vmin=pav_corr_ordered.min().min(), vmax=pav_corr_ordered.max().max(),
                cbar_kws={'orientation': 'horizontal'})
cnvs_corr[cnvs_corr < 0.4] = 0.4
cnv_corr_ordered = cnvs_corr.iloc[kin_row_order, kin_row_order]
d = sns.heatmap(cnv_corr_ordered, cmap="RdBu_r", center=0.8, square=True,
                ax=ax[1][1], xticklabels=False, yticklabels=False,
                vmin=cnv_corr_ordered.min().min(), vmax=cnv_corr_ordered.max().max(),
                cbar_kws={'orientation': 'horizontal'})
fig.tight_layout()
fig.savefig(
    "Scripts/Data_Vis/Section_1/Figure_1b-e_data_correlations5.pdf", dpi=300)
plt.close()
