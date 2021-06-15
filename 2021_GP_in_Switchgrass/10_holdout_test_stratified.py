import sys,os
import numpy as np
import random
import pandas as pd
import imp
from sklearn.model_selection import StratifiedKFold

file = sys.argv[1] # your pheno matrix
trait = sys.argv[2] # the target trait
number = int(sys.argv[3]) # the fold you want to hold out as test set
Y = pd.read_csv(file,header=0,index_col=0,sep=',')
Y_tem = Y.copy(deep=True)
Y_tem.index = range(1,(Y.shape[0]+1))
Y_tem['Class'] = 0
Y['Class'] = 1
df = pd.concat([Y,Y_tem],axis=0)
Tr_te = StratifiedKFold(n_splits=number, random_state=42,shuffle=True)
nn = 0
for train_index, test_index in Tr_te.split(df.loc[:,trait],df.Class):
	if nn == 0:
		train_ID = df.iloc[train_index]
		test_ID = df.iloc[test_index]
		nn = 1

test = test_ID[test_ID.Class==1].index.to_list()
out = open('Test.txt','w')
for t in test:
	out.write('%s\n'%t)

out.close()
