import sys,os
import pandas as pd
import datatable as dt
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.impute import KNNImputer
class KNNImputer_Ks(BaseEstimator, TransformerMixin):
	def __init__(self, *Ks):
		self.Ks = Ks
	def fit(self, X,Ks):
		D_imputer = {}		
		for k in [3,4,5,6,7]:
			imputer = KNNImputer(n_neighbors=k)
			D_imputer[k] = imputer.fit(X)			  
		return D_imputer
	def transform(self, X):
		Impute_train = {}
		for k in [3,4,5,6,7]:
			Impute_train[k] = pd.DataFrame(D_imputer[k].transform(X))
			Impute_train[k].index = X.index
			Impute_train[k].columns = X.columns 
			if k == 3:
				Imputed = Impute_train[k].copy(deep=True)
				Imputed.loc[:,:] = 0
			Imputed = Imputed.add(Impute_train[k],fill_value=0)
		return round(Imputed/5,0)

imputer_knn = KNNImputer_Ks()

files = sys.argv[1]
test_file = sys.argv[2]

df = dt.fread(files,header=True,sep=',')
df = df.to_pandas()
df = df.set_index(df.columns[0], drop=True)
df = df.dropna(axis=1,how='all')

df_all = df.copy()
print('Removing test instances to apply model on later...')
with open(test_file) as test:
	test_instances = test.read().splitlines()
	
try:
	test_df = df.loc[test_instances, :]
	df_training = df.drop(test_instances)
except:
	test_instances = [int(x) for x in test_instances]
	test_df = df.loc[test_instances, :]
	df_training = df.drop(test_instances)

df_training = df_training.apply(pd.to_numeric)
df_training = df_training.dropna(axis=1,how='all')
test_df = test_df.apply(pd.to_numeric)
test_df = test_df.loc[:,df.columns]

imputer_knn = KNNImputer_Ks()
D_imputer = imputer_knn.fit(df_training, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df_training)
df_knn = df_knn.round(4)
test_df_knn = imputer_knn.transform(test_df)
DF = pd.concat([df_knn,test_df_knn],axis=0)
DF.to_csv(files + '_imputed',index=True, header=True,sep=",")

MAF = []
for i in range(0,DF.shape[1]):
	table = DF.iloc[:,i].value_counts()
	if max(table) <= 383*0.95:
		MAF.append(DF.columns[i])

if len(MAF) > 0:
	df2 = DF.loc[:,MAF]
	df2.to_csv(sys.argv[1] + '_imputed_MAF',index=True, header=True,sep=",")

num = df_all.isna().sum()
num75 = num[num<=14]
num90 = num[num<=3]
if len(num75) > 0:
	df75 = df2.loc[:,set(num75.index) & set(MAF)]
	df75.to_csv(sys.argv[1].split('_50per')[0] + '_75per_imputed_MAF',index=True, header=True,sep=",")
if len(num90) > 0:
	df90 = df2.loc[:,set(num90.index) & set(MAF)]
	df90.to_csv(sys.argv[1].split('_50per')[0] + '_90per_imputed_MAF',index=True, header=True,sep=",")
