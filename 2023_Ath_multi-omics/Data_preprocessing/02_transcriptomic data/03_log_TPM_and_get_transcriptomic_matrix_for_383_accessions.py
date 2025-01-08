import h5py, numpy
import pandas as pd
import sys,os
import math


with open('Common_accessions_between_SNP_Methy_Exp_Grimm_pheno.txt') as f:
	accession = f.read().splitlines()

df = pd.read_csv('GSE80744_TPM.txt',header=0,index_col=0,sep='\t')
df = df.T
df.index = pd.Series(df.index).str.replace("X","")
df_res = df.loc[accession,:]
df_res.to_csv('GSE80744_TPM_%s_accessions.csv'%len(accession),index=True, header=True,sep=",")

df_res = df +1
df_res = df_res.apply(numpy.log)
df_res.to_csv('GSE80744_TPM_%s_accessions_loge_plus1.csv'%len(accession),index=True, header=True,sep=",")



