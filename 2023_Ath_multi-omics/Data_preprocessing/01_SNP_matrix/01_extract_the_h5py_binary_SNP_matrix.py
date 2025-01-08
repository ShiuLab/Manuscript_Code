import h5py, numpy
import pandas as pd
# the SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz
f = h5py.File('1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')
df = pd.DataFrame(f['snps'][:])
df.columns = f['accessions']
pos = pd.Series([])
for i in range(0,5):
	chr = chr_regions[i]
	position_chr = positions[chr[0]:chr[1]]
	loc = pd.Series(position_chr).apply(lambda x: 'Chr%s_%s'%(i+1,x))
	pos = pos.append(loc)

df.index = pos
df.to_csv('SNP_binary_matrix_all.csv',index=True, header=True,sep=",")
df.iloc[0,:].to_csv('All_SNPs_list.txt',index=False, header=False,sep=",")


