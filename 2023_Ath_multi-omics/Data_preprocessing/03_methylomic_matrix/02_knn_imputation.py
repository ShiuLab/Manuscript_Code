import pandas as pd
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

imputer_knn = KNNImputer_Ks()

df = pd.read_csv('Araport11_GB_methylation_383_accessions.csv', sep=',', index_col = 0, header = 0,low_memory=False)
df = df.dropna(axis=1,how='all')
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn.to_csv('Araport11_GB_methylation_383_accessions_knn_imputed.csv',index=True, header=True,sep=",")

df = pd.read_csv('Araport11_GB_mCG_strict_383_accessions.csv', sep=',', index_col = 0, header = 0,low_memory=False)
df = df.dropna(axis=1,how='all')
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn.to_csv('Araport11_GB_mCG_strict_383_accessions_knn_imputed.csv',index=True, header=True,sep=",")

df = pd.read_csv('Araport11_GB_mCHG_strict_383_accessions.csv', sep=',', index_col = 0, header = 0,low_memory=False)
df = df.dropna(axis=1,how='all')
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn.to_csv('Araport11_GB_mCHG_strict_383_accessions_knn_imputed.csv',index=True, header=True,sep=",")

df = pd.read_csv('Araport11_GB_mCHH_strict_383_accessions.csv', sep=',', index_col = 0, header = 0,low_memory=False)
df = df.dropna(axis=1,how='all')
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn.to_csv('Araport11_GB_mCHH_strict_383_accessions_knn_imputed.csv',index=True, header=True,sep=",")



