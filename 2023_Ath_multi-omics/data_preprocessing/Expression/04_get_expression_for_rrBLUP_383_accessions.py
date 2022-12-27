import h5py, numpy
import pandas as pd
import sys,os
import math
with open('../Models_for_Grimm_pheno/Common_accessions_between_SNP_Methy_Exp_Grimm_pheno.txt') as f:
	accession = f.read().splitlines()

df = pd.read_csv('GSE80744_TPM.txt',header=0,index_col=0,sep='\t')
df = df.T
df.index = pd.Series(df.index).str.replace("X","")
#df.to_csv('GSE80744_TPM_727_accessions.csv',index=True, header=True,sep=",")
# df_res = df +1
# df_res = df_res.apply(numpy.log)
# df_res.to_csv('GSE80744_TPM_727_accessions_loge_plus1.csv',index=True, header=True,sep=",")
# df_res = df_res.loc[accession,:]
# df_res.to_csv('GSE80744_TPM_383_accessions_loge_plus1.csv',index=True, header=True,sep=",")

df_res = df.loc[accession,:]
df_res.to_csv('GSE80744_TPM_383_accessions.csv',index=True, header=True,sep=",")


