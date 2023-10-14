import pandas as pd
import sys,os
import numpy as np
sys.stdout.flush()

n = 0
for files in os.listdir('./'):
	if files.endswith('_methylated.txt'):
		df = pd.read_csv(files,index_col=0,header=0,sep='\t')
		if n == 0:
			res = df
			n = n + 1
		else:
			res = pd.concat([res,df], axis=1, join='outer')
			n = n + 1
		print(files, flush=True)
		print(res.shape,flush=True)

res = res.fillna(0)
res = res.astype(int)
res.to_csv('Methylation_genome_wide_383_accessions.csv',index=True, header=True,sep=",")
