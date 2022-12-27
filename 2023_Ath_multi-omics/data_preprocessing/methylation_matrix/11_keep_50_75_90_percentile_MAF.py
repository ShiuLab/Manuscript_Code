import pandas as pd
import sys,os
import numpy as np
sys.stdout.flush()
file = sys.argv[1]
df =  pd.read_csv(file,index_col=0,header=0,sep=',')
MAF = []
for i in range(0,df.shape[1]):
	table = df.iloc[:,i].value_counts()
	if max(table) <= 383*0.95:
		MAF.append(df.columns[i])

df2 = df.loc[:,MAF]

num = df2.isna().sum()
num50 = num[num<=38]
num75 = num[num<=8]
num90 = num[num<=2]
df50 = df2.loc[:,num50.index]
df75 = df2.loc[:,num75.index]
df90 = df2.loc[:,num90.index]
df50.to_csv(file + '_50per',index=True, header=True,sep=",")
df75.to_csv(file + '_75per',index=True, header=True,sep=",")
df90.to_csv(file + '_90per',index=True, header=True,sep=",")

