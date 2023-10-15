import pandas as pd
import sys,os
import numpy as np
sys.stdout.flush()
file = sys.argv[1]
df =  pd.read_csv(file,index_col=0,header=0,sep=',')
num = df.isna().sum()
num = pd.DataFrame(num)
num.to_csv(file + '_count',index=True, header=None,sep="\t")

