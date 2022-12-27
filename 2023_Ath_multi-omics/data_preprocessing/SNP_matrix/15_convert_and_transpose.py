import pandas as pd
import sys,os

df_MAF = pd.read_csv(sys.argv[1],header=0,index_col=0,sep=',')
remain = df_MAF.sum(axis=1)
index = remain.index.tolist()
df_convert = df_MAF.copy(deep=True)
for i in range(0,remain.shape[0]):
	if remain.iloc[i] <= 383*0.5: # more common is 0, less common is 1
		df_convert.loc[index[i],df_convert.loc[index[i],:]==1] = -1
		df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = 1
	else: # more common is 1, less common is 0
		df_convert.loc[index[i],df_convert.loc[index[i],:]==0] = -1
df_convert = df_convert.T
df_convert.to_csv(sys.argv[1].split('.csv')[0] + '_converted.csv',index=True, header=True,sep=",")


