import sys,os
import pandas as pd
import numpy as np
import pickle
files = sys.argv[1]
df = pd.read_csv(files,index_col=0,header=0,sep=',')
for accession in df.index.tolist():
	with open('%s_met.txt.pkl'%accession, 'rb') as f1:
		D = pickle.load(f1)
	for marker in df.columns.tolist():
		# load dictionary
		if marker in D:
			df.loc[accession,marker] = D[marker]

df = df.round(3)
df.to_csv(files + '_proportion',index=True, header=True,sep=",")

