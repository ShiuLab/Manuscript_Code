import h5py, numpy
import pandas as pd
import sys,os

# the SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz
f = h5py.File('1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')

# the 383 accession ids can be found in the Datasets folder 
with open('Common_accessions_between_SNP_Methy_Exp_Grimm_pheno.txt') as F:
	accession = F.read().splitlines()
	
accs = f['accessions'][:]
accs_new = []
for a in accs:
	accs_new.append(a.decode('UTF-8'))
	
accs_new = numpy.array(accs_new)
ix_of_acc = []
for a in accession:
	ix_of_acc.append(numpy.where(accs_new==a)[0][0])
	
# Get all SNP positions for all chromosomes (len=10709949)
positions = f['positions'][:]
# Array of tupels with start/stop indices for each chromosome
chr_regions = f['positions'].attrs['chr_regions']
snp = f['snps'][:]
res = snp[:,ix_of_acc]
df = pd.DataFrame(res)
df.columns = accession
pos = pd.Series([])
for i in range(0,5):
	chr = chr_regions[i]
	position_chr = positions[chr[0]:chr[1]]
	loc = pd.Series(position_chr).apply(lambda x: 'Chr%s_%s'%(i+1,x))
	pos = pos.append(loc)
	
df.index = pos
df.to_csv('SNP_binary_matrix_%s_accessions.csv'%len(accession),index=True, header=True,sep=",")
#drop rows with all zero
df = df.loc[~(df==0).all(axis=1)]
df.to_csv('SNP_binary_matrix_%s_accessions_drop_all_zero.csv'%len(accession),index=True, header=True,sep=",")
#drop rows with MAF < 0.05
summ = df.sum(axis=1)
remain = summ[summ >= len(accession)*0.05]
remain = remain[remain <= len(accession)*0.95]
df_MAF = df.loc[remain.index,:]
df_MAF.to_csv('SNP_binary_matrix_%s_accessions_drop_all_zero_MAF_larger_than_0.05.csv'%len(accession),index=True, header=True,sep=",")
# replace more common allele as 1, less common allele as -1
index = remain.index.tolist()
df_convert = df_MAF.copy(deep=True)
for i in range(0,remain.shape[0]):
for i in range(0,remain.shape[0]):
	if remain.iloc[i] <= len(accession)*0.5: # more common is 0, less common is 1
		df_convert.loc[index[i],df_convert.loc[index[i],:]==1] = -1
		df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = 1
	else: # more common is 1, less common is 0
		df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = -1
df_convert = df_convert.T
df_convert.to_csv('SNP_binary_matrix_%s_accessions_drop_all_zero_MAF_larger_than_0.05_converted.csv'%len(accession),index=True, header=True,sep=",")
