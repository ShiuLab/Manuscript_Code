import h5py, numpy
import pandas as pd
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
df.to_csv('SNP_binary_matrix.csv',index=True, header=True,sep=",")


