import pandas as pd
import numpy as np
geno_file = sys.argv[1]
df = pd.read_csv(geno_file,header=0,index_col=0,sep=',')
df = df.loc[:, (df != df.iloc[0]).any()] 
df.to_csv(geno_file.split('.csv')[0] + '_non_constant.csv',header=True, index=True, sep=',')
