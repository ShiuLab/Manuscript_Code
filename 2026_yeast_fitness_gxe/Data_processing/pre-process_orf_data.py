#!/usr/bin/env python3

'''Pre-processing of raw open reading frame data from Peter et al., 2018'''

import os
import pandas as pd
import matplotlib.pyplot as plt

# Keep only isolates with phenotype information
cnv = pd.read_csv('Data/Peter_2018/0_raw_data/genesMatrix_CopyNumber.tab_diploid.txt', sep='\t')
pheno = pd.read_csv('Data/Peter_2018/pheno.csv') # S1 File fitness data
out = cnv[cnv.STD_name.isin(pheno['ID'])]
out.rename(columns={'STD_name': 'ID'}, inplace=True)

# Remove CNVs with missing values (out of 7,796 ORFs, kept 7,708)
cnvs_to_drop = out.columns[out.isna().sum() == 750]

# S6 File CNV data
out.dropna(axis=1).to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv")

# Remove PAVs with missing values
pav = pd.read_csv('/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/0_raw_data/genesMatrix_PresenceAbsence.tab', sep='\t')
pav = pav[pav.STD_name.isin(pheno['ID'])] # diploid isolates only
pav.rename(columns={'STD_name': 'ID'}, inplace=True)
pav.isna().sum().apply(lambda x: x/750).sum() # 88 PAVs missing in all 750 isolates
pav.columns[pav.isna().sum() == 750].isin(cnvs_to_drop).sum() # the same CNVs are missing in PAVs

# S5 File PAV data
pav.dropna(axis=1).to_csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv", index=False)
