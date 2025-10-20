#!/usr/bin/env python3

import os
import datatable as dt
import pandas as pd
from sklearn.decomposition import PCA # sklearn v1.2.2

# Load in snp data
os.chdir('/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project')
geno = dt.fread('Data/Peter_2018/geno.csv').to_pandas() # S2 File SNP data
test = pd.read_csv('Data/Peter_2018/Test.txt', header=None)

# Split train and test data; fit PCA on train data
train = geno[~geno['ID'].isin(test[0])].set_index('ID')
test = geno[geno['ID'].isin(test[0])].set_index('ID')

pca = PCA(n_components=624)
pca.fit(train)
pc_train = pca.transform(train)
pc_test = pca.transform(test)

# Save PCs to file
out = pd.concat([
    pd.DataFrame(pc_train, index=train.index, columns=[f'PC{i+1}' for i in range(pc_train.shape[1])]),
    pd.DataFrame(pc_test, index=test.index, columns=[f'PC{i+1}' for i in range(pc_test.shape[1])])],
    axis=0)
out.to_csv('Data/Peter_2018/PCs.csv')
out.iloc[:,:5].to_csv('Data/Peter_2018/PCs_first5.csv')
sum(pca.explained_variance_ratio_[:5]) # 0.5911324069476012

pd.DataFrame(pca.explained_variance_ratio_,
    index=[f'PC{i+1}' for i in range(pca.explained_variance_ratio_.shape[0])],
    columns=['Explained Variance Ratio']).to_csv('Data/Peter_2018/PCs_variance.csv')
