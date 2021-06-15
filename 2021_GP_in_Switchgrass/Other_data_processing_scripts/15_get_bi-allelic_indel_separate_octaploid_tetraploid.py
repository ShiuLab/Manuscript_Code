import sys,os
import pandas as pd
'''
input1: classification of SNPs
input2: file about ploidy of individuals
input3: SNP matrix
'''
type = open(sys.argv[1],'r').readlines()
ploidy = open(sys.argv[2],'r').readlines()
inp = open(sys.argv[3],'r').readlines()
P = {}
for inl in ploidy:
	tem = inl.split('\t')
	ind = tem[0]
	p = int(tem[2])
	P[ind] = p

D = {}	
for inl in type[1:]:
	tem = inl.split('\t')
	chr = tem[0]
	pos = tem[1]
	type = tem[4]
	allelic = int(tem[5])
	if type == 'indel' and allelic == 2:
		D[chr + '-' + pos] = 1

out = open(sys.argv[3].split('.txt')[0] + '_biallelic_indel.txt','w')
out.write(inp[0])
for inl in inp[1:]:
	tem = inl.split('\t')
	chr = tem[0]
	pos = tem[1]
	if chr + '-' + pos in D:
		out.write(inl)

out.close()
		
df = pd.read_csv(sys.argv[3].split('.txt')[0] + '_biallelic_indel.txt', sep='\t', index_col = 0, header = 0)
colname = df.columns.tolist()
df_4 = df.copy()
df_8 = df.copy()
for col in colname[3:]:
	if P[col] == 8:
		del df_4[col]
	if P[col] == 4:
		del df_8[col]
		
df_4.to_csv(sys.argv[3].split('.txt')[0] + '_biallelic_indel_tetraploid.txt', index=True, header=True,sep="\t")
df_8.to_csv(sys.argv[3].split('.txt')[0] + '_biallelic_indel_octaploid.txt', index=True, header=True,sep="\t")