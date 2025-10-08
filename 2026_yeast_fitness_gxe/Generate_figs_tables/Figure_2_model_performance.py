#!/usr/bin/env python3
############################################################################
# Figure 2
############################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Single environment RF (after feature selection) model performances
snp = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_SNPs_FS.txt", sep="\t")
pav = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_PAVs_FS.txt", sep="\t")
cnv = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_CNVs_FS.txt", sep="\t")
pc = pd.read_csv(
    "Scripts/Data_Vis/Section_2/RESULTS_RF_PCs_sorted.txt", sep="\t")

snp = snp.set_index('new_cond').loc[pc['new_cond']]
pav = pav.set_index('new_cond').loc[pc['new_cond']]
cnv = cnv.set_index('new_cond').loc[pc['new_cond']]

pc.r2_test.mean()  # 0.19543497114629035
pc.r2_test.median()  # 0.1463634586478904
cnv.r2_test.median()  # 0.2064282441291797
snp.r2_test.median()  # 0.2007607748189971
pav.r2_test.median()  # 0.184326612703762
cnv.r2_test.max()  # 0.7021635825306303

# How many envs were SNPs, PAVs, and CNVs better than PCs?
snp_diff = snp.r2_test - pc.set_index("new_cond").r2_test
pav_diff = pav.r2_test - pc.set_index("new_cond").r2_test
cnv_diff = cnv.r2_test - pc.set_index("new_cond").r2_test
sum(snp_diff > 0)  # 25
sum(pav_diff > 0)  # 22
sum(cnv_diff > 0)  # 22

# How many envs did CNVs perform better than SNPs, PAVs, and PCs?
sum((cnv.r2_test - snp.r2_test) > 0)  # 18
sum((cnv.r2_test - pav.r2_test) > 0)  # 18
sum((cnv.r2_test - pc.set_index('new_cond').r2_test) > 0)  # 22
# 1.561012 for sodium meta-arsenite and 1.661651 for cuso4
cnv.r2_test / pc.set_index('new_cond').r2_test
cnv.loc[(cnv.r2_test > pav.r2_test) & (cnv.r2_test > snp.r2_test) &
        (cnv.r2_test > pc.set_index('new_cond').r2_test), :].shape[0]  # 13 out of 35

# How many envs did SNPs perform better than PAVs, CNVs, and PCs?
snp.loc[(snp.r2_test > pav.r2_test) & (snp.r2_test > cnv.r2_test) &
        (snp.r2_test > pc.set_index('new_cond').r2_test), :].shape[0]  # 12

# How many envs did PAVs perform better than SNPs, CNVs, and PCs?
pav.loc[(pav.r2_test > snp.r2_test) & (pav.r2_test > cnv.r2_test) &
        (pav.r2_test > pc.set_index('new_cond').r2_test), :].shape[0]  # 8

# How many envs did PCs perform better than SNPs, PAVs, and CNVs?
pc.set_index('new_cond').loc[(pc.set_index('new_cond').r2_test > snp.r2_test) &
                             (pc.set_index('new_cond').r2_test > pav.r2_test) &
                             (pc.set_index('new_cond').r2_test > cnv.r2_test), :].shape[0]  # 2

# How many envs did SNPs and PAVs outperform CNVs?
snp.loc[(snp.r2_test > cnv.r2_test)].shape[0]  # 17
pav.loc[(pav.r2_test > cnv.r2_test)].shape[0]  # 17
cnv.loc[(cnv.r2_test > snp.r2_test)].shape[0]  # 18
cnv.loc[(cnv.r2_test > pav.r2_test)].shape[0]  # 18

# Do environments that clustered together by fitness values have similar
# model performance?
pc.set_index("new_cond", inplace=True)
r2_clust = hierarchy.linkage(pdist(pd.concat(
    [pc.r2_test, snp.r2_test, pav.r2_test, cnv.r2_test],
    ignore_index=False, axis=1), metric='euclidean'), method='complete')

ordered_linkage = hierarchy.optimal_leaf_ordering(r2_clust, pdist(pd.concat(
    [pc.r2_test, snp.r2_test, pav.r2_test, cnv.r2_test],
    ignore_index=False, axis=1), metric='euclidean'))

hierarchy.dendrogram(ordered_linkage, labels=pc.index.tolist(),
                     leaf_font_size=8, color_threshold=0.5, orientation='right')
plt.title("Hierarchical Clustering of Model Performance")
plt.xlabel("Euclidean Distance")
plt.tight_layout()
plt.savefig(
    "Scripts/Data_Vis/Section_2/model_performance_dendrogram.pdf")
plt.close()

############################################################################
# Figure 2
############################################################################
h2 = pd.read_csv("Data/Peter_2018/Heritability_h2_H2_sommer.csv")
h2 = pd.concat([h2.set_index("Conditions"), pc.reset_index().set_index("Y")["new_cond"]],
               axis=1, ignore_index=False).reset_index()

h2.set_index("new_cond", inplace=True)

# sort values by the max r2_test of at least one variant type
order = pd.concat([pc.r2_test, snp.r2_test, pav.r2_test, cnv.r2_test], axis=1,
                  ignore_index=False).max(axis=1).sort_values(ascending=False)
order.to_csv("Scripts/Data_Vis/Section_2/Figure_2_r2_test_v5_env_order.txt",
             sep="\t", header=False)

h2 = h2.loc[order.index]
pc = pc.loc[h2.index]
snp = snp.loc[h2.index]
pav = pav.loc[h2.index]
cnv = cnv.loc[h2.index]

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5.5, 7.5))
y_pos = np.arange(35)
bar_height = 0.35
ax.barh(y=y_pos - bar_height/2,
        width=h2["h2"], height=bar_height, color="white", edgecolor="black", label="h²")
ax.barh(y=y_pos + bar_height/2, width=pc.r2_test, height=bar_height,
        color="lightgray", edgecolor="black", label="Test R² (PCs)")
ax.scatter(x=snp.r2_test, y=pc.index, color="tomato", label="SNP")
ax.scatter(x=pav.r2_test, y=pc.index, color="deepskyblue", label="PAV")
ax.scatter(x=cnv.r2_test, y=pc.index, color="mediumorchid", label="CNV")
ax.set_xlabel("Value")
ax.legend(loc="upper right")
plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/Figure_2_h2_explained_v5.pdf")
plt.close()
