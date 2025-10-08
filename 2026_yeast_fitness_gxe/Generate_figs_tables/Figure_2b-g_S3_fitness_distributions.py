#!/bin/usr/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

pheno = pd.read_csv("Data/Peter_2018/pheno.csv", index_col=0)  # fitness data
sorted_cols = sorted(pheno.columns.tolist())
pheno = pheno[sorted_cols]

test = pd.read_csv("Data/Peter_2018/Test.txt", header=None)  # test set data
pheno_train = pheno.loc[~pheno.index.isin(test[0].tolist())]
pheno_test = pheno.loc[pheno.index.isin(test[0].tolist())]

# Figures 2B-G and S3. Training and test set distributions on one plot
fig, ax = plt.subplots(7, 5, figsize=(10.5, 13), sharex=True)
ax = ax.flatten()
for i, col in enumerate(pheno.columns):
    pheno_train[col].plot.hist(
        bins=30, alpha=0.5, ax=ax[i], label="Training set", color="#38cbff")
    pheno_test[col].plot.hist(bins=30, alpha=0.5, ax=ax[i], label="Test set",
                              color="#0053a1")
    ax[i].set_xlabel("Fitness")
    ax[i].set_ylabel("Frequency")
    ax[i].set_title(col)


plt.tight_layout()
plt.savefig("Scripts/Data_Vis/Section_2/Figure_S3_train_test_pheno.pdf",
            bbox_inches="tight")
plt.close()
