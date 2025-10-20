# !/usr/bin/env python3
# Filter phenotype data to only include diploid isolates

import pandas as pd

# Read in the raw fitness data
pheno = pd.read_csv('phenoMatrix_35ConditionsNormalizedByYPD.tab', sep='\t')
pheno = pheno.set_index(pheno.iloc[:, 0])  # set isolate IDs as index

# Read in the ploidy information
labels = pd.read_excel(
    'Peter_2018_Supplementary_Tables.xls', sheet_name='Table S1', skiprows=3)

# select on diploid isolates
diploid = labels.loc[labels['Ploidy'] == 2]
print("diploid", diploid.shape)  # 787 isolates

# filter diploids with phenotype data
pheno_diploid = pheno[pheno.index.isin(diploid['Standardized name'])]  # 750
pheno_diploid.index.name = 'ID'
pd.Series(pheno_diploid.index).to_csv(
    'diploid_750_IDs.txt', index=False, header=False)
pheno_diploid.iloc[:, 1:].to_csv('pheno.csv')
