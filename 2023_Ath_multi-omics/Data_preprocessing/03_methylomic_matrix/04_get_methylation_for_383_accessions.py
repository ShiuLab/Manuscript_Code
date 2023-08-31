import h5py, numpy
import pandas as pd
import sys,os

# the 383 accession ids can be found in the Datasets folder 
with open('Common_accessions_between_SNP_Methy_Exp_Grimm_pheno.txt') as f:
	accession = f.read().splitlines()

# the 1107 accession ids for methylomic data can be found in the Datasets folder 
with open('Accession_ID_for_methylation_1107.txt') as f:
	Methy_accession = f.read().splitlines()

df = pd.read_csv('Araport11_GB_mCG_strict.tsv',header=None,index_col=0,sep='\t')
df.columns = Methy_accession
rep = df.columns.value_counts()
df_new = pd.DataFrame(index=df.index, columns=['0'])
for col in rep.index:
	subdf = df.loc[:,col]
	if rep.loc[col] > 1:
		rowmedian = pd.DataFrame(subdf.median(axis=1,skipna=True))
		rowmedian.columns = [col]
	else:
		rowmedian = subdf
	df_new = pd.concat([df_new,rowmedian],axis=1)
	
df_new = df_new.drop(labels='0',axis=1)

df_new.to_csv('Araport11_GB_mCG_strict_1028_accessions.csv',index=True, header=True,sep=",")

df_res = df_new.loc[:,accession]
df_res = df_res.T
df_res = df_res.add_prefix('GB_mCG_')
df_res.to_csv('Araport11_GB_mCG_strict_383_accessions.csv',index=True, header=True,sep=",")

df = pd.read_csv('Araport11_GB_mCHG_strict.tsv',header=None,index_col=0,sep='\t')
df.columns = Methy_accession
rep = df.columns.value_counts()
df_new = pd.DataFrame(index=df.index, columns=['0'])
for col in rep.index:
	subdf = df.loc[:,col]
	if rep.loc[col] > 1:
		rowmedian = pd.DataFrame(subdf.median(axis=1,skipna=True))
		rowmedian.columns = [col]
	else:
		rowmedian = subdf
	df_new = pd.concat([df_new,rowmedian],axis=1)
	
df_new = df_new.drop(labels='0',axis=1)

df_new.to_csv('Araport11_GB_mCHG_strict_1028_accessions.csv',index=True, header=True,sep=",")

df_res = df_new.loc[:,accession]
df_res = df_res.T
df_res = df_res.add_prefix('GB_mCHG_')
df_res.to_csv('Araport11_GB_mCHG_strict_383_accessions.csv',index=True, header=True,sep=",")

df = pd.read_csv('Araport11_GB_mCHH_strict.tsv',header=None,index_col=0,sep='\t')
df.columns = Methy_accession
rep = df.columns.value_counts()
df_new = pd.DataFrame(index=df.index, columns=['0'])
for col in rep.index:
	subdf = df.loc[:,col]
	if rep.loc[col] > 1:
		rowmedian = pd.DataFrame(subdf.median(axis=1,skipna=True))
		rowmedian.columns = [col]
	else:
		rowmedian = subdf
	df_new = pd.concat([df_new,rowmedian],axis=1)
	
df_new = df_new.drop(labels='0',axis=1)

df_new.to_csv('Araport11_GB_mCHH_strict_1028_accessions.csv',index=True, header=True,sep=",")

df_res = df_new.loc[:,accession]
df_res = df_res.T
df_res = df_res.add_prefix('GB_mCHH_')
df_res.to_csv('Araport11_GB_mCHH_strict_383_accessions.csv',index=True, header=True,sep=",")


mCG = pd.read_csv('Araport11_GB_mCG_strict_383_accessions.csv',header=0,index_col=0,sep=',')
mCHG = pd.read_csv('Araport11_GB_mCHG_strict_383_accessions.csv',header=0,index_col=0,sep=',')
mCHH = pd.read_csv('Araport11_GB_mCHH_strict_383_accessions.csv',header=0,index_col=0,sep=',')
res = pd.concat([mCG,mCHG],axis=1)
res = pd.concat([res,mCHH],axis=1)
res.to_csv('Araport11_GB_methylation_383_accessions.csv',index=True, header=True,sep=",")














