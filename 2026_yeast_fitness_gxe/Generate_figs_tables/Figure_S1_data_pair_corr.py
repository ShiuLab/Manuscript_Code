#!/usr/bin/env python3
################################################################################
# Figure S1
################################################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from scipy.stats import spearmanr, ttest_ind

# Read in isolate pair correlations
clades = pd.read_excel("Data/Peter_2018/0_raw_data/Peter_2018_Supplementary_Tables.xls",
                       sheet_name="Table S1", skiprows=3, nrows=1011)
kinship = pd.read_csv("Data/Peter_2018/kinship.csv", index_col=0)
pheno_corr = pd.read_csv(
    "Data/Peter_2018/fitness_correlations_isolate_pair.csv", index_col=0)
pavs_corr = pd.read_csv(
    "Data/Peter_2018/ORFs_pres_abs_correlations_isolate_pair.csv", index_col=0)
cnvs_corr = pd.read_csv(
    "Data/Peter_2018/ORFs_copy_num_correlations_isolate_pair.csv", index_col=0)

# Remove duplicate pairs ex. AAC AAC; AVD APE and APE AVD
kin_melt = kinship.where(np.triu(kinship, k=1).astype(bool)).stack()
pheno_melt = pheno_corr.where(np.triu(pheno_corr, k=1).astype(bool)).stack()
pavs_melt = pavs_corr.where(np.triu(pavs_corr, k=1).astype(bool)).stack()
cnvs_melt = cnvs_corr.where(np.triu(cnvs_corr, k=1).astype(bool)).stack()

# Sanity check (order matters for spearman's rank correlation)
sum(kin_melt.index == pheno_melt.index)
sum(kin_melt.index == pavs_melt.index)
sum(kin_melt.index == cnvs_melt.index)
kin_melt = kin_melt.reset_index()
pheno_melt = pheno_melt.reset_index()
pavs_melt = pavs_melt.reset_index()
cnvs_melt = cnvs_melt.reset_index()
kin_melt["rank"] = kin_melt.rank(numeric_only=True, ascending=False)
pheno_melt["rank"] = pheno_melt.rank(numeric_only=True, ascending=False)
pavs_melt["rank"] = pavs_melt.rank(numeric_only=True, ascending=False)
cnvs_melt["rank"] = cnvs_melt.rank(numeric_only=True, ascending=False)
kin_melt.to_csv("Scripts/Data_Vis/kin_melt.csv")
pheno_melt.to_csv("Scripts/Data_Vis/pheno_melt.csv")
pavs_melt.to_csv("Scripts/Data_Vis/pavs_melt.csv")
cnvs_melt.to_csv("Scripts/Data_Vis/cnvs_melt.csv")


def normalized_cbar(z):
    """ Normalize colorbar to density values """
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    return norm


def density_scatter(x, y, ax=None, fig=None, sort=True, bins=20):
    """
    Calculate density between x and y vectors and plot as a scatter plot
    source: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density/53865762#53865762
    """
    if ax is None:
        fig, ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])),
                data, np.vstack([x, y]).T, method="splinef2d", bounds_error=False)
    z[np.where(np.isnan(z))] = 0.0  # to be sure to plot all data
    if sort:  # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    # return z
    if ax is not None:
        ax.scatter(x, y, c=z, s=0.3)
        fig.colorbar(cm.ScalarMappable(
            norm=normalized_cbar(z), cmap="viridis"), ax=ax)
        return ax
    else:
        plt.scatter(x, y, c=z, s=0.3)
        fig.colorbar(cm.ScalarMappable(
            norm=normalized_cbar(z), cmap="viridis"))
        return fig


data_pair_corr = pd.DataFrame(columns=["comparison", "rho", "p-value", "note"])
# Plot
fig, ax = plt.subplots(2, 3, figsize=(12, 6))
ax[0][0] = density_scatter(kin_melt[0], pheno_melt[0], ax[0][0], fig, bins=50)
# calculate spearman's rho
rho, p = spearmanr(kin_melt['rank'], pheno_melt['rank'])
ttest_ind(kin_melt[0], pheno_melt[0], equal_var=False,
          alternative='two-sided')  # p-values are all 0
# rho = .269, p-val = 0.000E+00 ; in R, cor.test p-values are all < 2.2e-16
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))

ax[0][1] = density_scatter(pavs_melt[0], pheno_melt[0], ax[0][1], fig, bins=50)
rho, p = spearmanr(pavs_melt['rank'], pheno_melt['rank'])
ttest_ind(pavs_melt[0], pheno_melt[0],
          equal_var=False, alternative='two-sided')
# rho = .299, p-val = 0.000E+00
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))

ax[0][2] = density_scatter(cnvs_melt[0], pheno_melt[0], ax[0][2], fig, bins=50)
rho, p = spearmanr(cnvs_melt['rank'], pheno_melt['rank'])
ttest_ind(cnvs_melt[0], pheno_melt[0],
          equal_var=False, alternative='two-sided')
# rho = 0.138, p-val = 0.000E+00
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))

ax[1][0] = density_scatter(pavs_melt[0], kin_melt[0], ax[1][0], fig, bins=50)
rho, p = spearmanr(pavs_melt['rank'], kin_melt['rank'])
ttest_ind(pavs_melt[0], kin_melt[0],
          equal_var=False, alternative='two-sided')
# rho = 0.475, p-val = 0.000E+00
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))

ax[1][1] = density_scatter(cnvs_melt[0], kin_melt[0], ax[1][1], fig, bins=50)
rho, p = spearmanr(cnvs_melt['rank'], kin_melt['rank'])
ttest_ind(cnvs_melt[0], kin_melt[0],
          equal_var=False, alternative='two-sided')
# rho = 0.096, p-val = 0.000E+00
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))

ax[1][2] = density_scatter(cnvs_melt[0], pavs_melt[0], ax[1][2], fig, bins=50)
rho, p = spearmanr(cnvs_melt['rank'], pavs_melt['rank'])
ttest_ind(cnvs_melt[0], pavs_melt[0],
          equal_var=False, alternative='two-sided')
# rho = 0.172, p-val = 0.000E+00
print('rho = {:.3f}\np-value = {:.3E}'.format(rho, p))
fig.tight_layout()
fig.savefig(
    "Scripts/Data_Vis/Section_1/Figure_S1_data_pair_correlations.png", dpi=300)
plt.close()
