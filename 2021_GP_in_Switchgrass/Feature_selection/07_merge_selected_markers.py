import sys,os
import pandas as pd
n = 0
for files in os.listdir('./'):
	if (files.startswith('GBS') or files.startswith('Exome')) and files.endswith('csv'):
		df = pd.read_csv(files, sep=',', index_col = 0, header = 0)  
		df = df.add_prefix('_'.join(files.split('_')[0:3]) + '_')
		if n == 0:
			res = df
		else:
			res = res.merge(df, left_index=True, right_index=True)
		n += 1

res.to_csv('geno.csv', index=True, header=True,sep=",")