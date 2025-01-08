import pandas as pd
import sys,os
import numpy as np
sys.stdout.flush()
import pickle
file = sys.argv[1]
df = pd.read_csv(file,index_col=0,header=0,sep=',')
df = df.T
df.to_csv('%s_transposed'%file,index=True, header=True,sep=",")

with open('/mnt/home/peipeiw/Documents/Ath_GS/Models_for_Grimm_pheno/Methylation_data_all_sites/Methylation_genome_wide_383_accessions_overlapping_with_SNP_targeted.pkl', 'rb') as f1:
	D = pickle.load(f1)


inp =  open('%s_transposed'%file,'r')
out = open('%s_transposed_filled_with_NaN.csv'%file,'w')
inl = inp.readline()
out.write(inl)
out.flush()
L = {} # methylation sites
tem = inl.strip().split(',')
for i in range(1,len(tem)):
	L[i] = tem[i]

inl = inp.readline()
while inl:
	tem = inl.strip().split(',')
	accession = tem[0]
	with open('%s_met.txt'%accession) as f:
		r = f.read().splitlines()
	R = dict.fromkeys(r,'1')
	print(len(R), flush=True)
	for i in range(1,len(tem)):
		if tem[i] == '0':
			if L[i] not in R: # not in the original methylation files
				if L[i] in D:
					if tem[0] not in D[L[i]]:
						tem[i] = ''
				else:
					tem[i] = ''
		if (i%1000 == 0):
			print(i,flush=True)
	out.write(','.join(tem) + '\n')
	out.flush()
	inl = inp.readline()

out.close()


