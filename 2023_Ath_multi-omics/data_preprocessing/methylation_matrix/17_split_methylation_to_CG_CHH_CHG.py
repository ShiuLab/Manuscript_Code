import sys,os
import pandas as pd
sys.stdout.flush()
file = sys.argv[1]
df =  pd.read_csv(file,index_col=0,header=0,sep=',')
col = df.columns.tolist()
CG = []
CHH = []
CHG = []
for c in col:
	met = c.split('_')[2]
	if len(met) == 2 and met == 'CG':
		CG.append(c)
	if len(met) == 3 and met[2] == 'G' and met[1] != 'G':
		CHG.append(c)
	if len(met) == 3 and met[2] != 'G' and met[1] != 'G':
		CHH.append(c)

if len(CG) > 0:
	df_CG = df.loc[:,CG]
	df_CG = df_CG.round(3)
	df_CG.to_csv(file + '_CG',index=True, header=True,sep=",")
if len(CHH) > 0:
	df_CHH = df.loc[:,CHH]
	df_CHH = df_CHH.round(3)
	df_CHH.to_csv(file + '_CHH',index=True, header=True,sep=",")
if len(CHG) > 0:
	df_CHG = df.loc[:,CHG]
	df_CHG = df_CHG.round(3)
	df_CHG.to_csv(file + '_CHG',index=True, header=True,sep=",")
