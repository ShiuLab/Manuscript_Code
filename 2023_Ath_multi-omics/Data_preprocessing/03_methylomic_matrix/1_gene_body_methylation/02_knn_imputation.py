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
		return Imputed/5

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
	df = df.drop(test_instances)
except:
	test_instances = [int(x) for x in test_instances]
	test_df = df.loc[test_instances, :]
	df = df.drop(test_instances)

df = df.apply(pd.to_numeric)
df = df.dropna(axis=1,how='all')
test_df = test_df.apply(pd.to_numeric)
test_df = test_df.loc[:,df.columns]

imputer_knn = KNNImputer_Ks()
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn = df_knn.round(4)
test_df_knn = imputer_knn.transform(test_df)
DF = pd.concat([df_knn,test_df_knn],axis=0)
DF.to_csv(files + '_imputed',index=True, header=True,sep=",")

