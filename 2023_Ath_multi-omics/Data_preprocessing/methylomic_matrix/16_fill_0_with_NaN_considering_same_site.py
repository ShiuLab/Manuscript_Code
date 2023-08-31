import pandas as pd
import sys,os
import numpy as np
sys.stdout.flush()
import pickle
file = sys.argv[1]
# df = pd.read_csv(file,index_col=0,header=0,sep=',')
# df = df.T
# df.to_csv('%s_transposed'%file,index=True, header=True,sep=",")

# with open('/mnt/home/peipeiw/Documents/Ath_GS/Models_for_Grimm_pheno/Methylation_data_all_sites/Methylation_genome_wide_383_accessions_overlapping_with_SNP_targeted.pkl', 'rb') as f1:
	# D = pickle.load(f1)
'''
rep = open('/mnt/home/peipeiw/Documents/Ath_GS/Models_for_Grimm_pheno/Methylation_data_all_sites/Methylation_sites_listgenome_wide_383_accessions_ordered.txt','r').readlines()
Rep = {}
for rep_l in rep:
	tem = rep_l.strip().split('\t')
	if '%s_%s'%(tem[1],tem[2]) not in Rep:
		Rep['%s_%s'%(tem[1],tem[2])] = []
	Rep['%s_%s'%(tem[1],tem[2])].append(tem[0])
	
import pickle
f = open("Methylation_sites_listgenome_wide_383_accessions_ordered.pkl","wb")
pickle.dump(Rep,f)
f.close()
'''

import pickle
with open('Methylation_sites_listgenome_wide_383_accessions_ordered.pkl', 'rb') as f1:
	Rep = pickle.load(f1)


inp =  open('%s_transposed_filled_with_NaN.csv'%file,'r')
out = open('%s_transposed_filled_with_NaN_considering_same_site.csv'%file,'w')
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
	#os.system('awk "NR==FNR { lines[$0]=1; next } $0 in lines" %s_selected_columns %s_transposed > %_%_met.txt'%(accession,file,accession,file.split('_')[-1]))
	with open('%s_met.txt'%accession) as f:
		r = f.read().splitlines()
	R = dict.fromkeys(r,'1')
	print(len(R), flush=True)
	for i in range(1,len(tem)):
		if tem[i] == '':
			loc = '%s_%s'%(L[i].split('_')[0],L[i].split('_')[1])
			if len(Rep[loc]) > 1:
				for site in Rep[loc]:
					if site in R and site != L[i]:
						tem[i] = '0'
						print([accession,L[i]],flush=True)
		if (i%1000 == 0):
			print(i,flush=True)
	out.write(','.join(tem) + '\n')
	out.flush()
	inl = inp.readline()
	print(accession,flush=True)

out.close()


