import sys,os
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
		return round(Imputed/5,0)

imputer_knn = KNNImputer_Ks()


df = pd.read_csv(sys.argv[1], sep=',', index_col = 0, header = 0,low_memory=False)

df = df.dropna(axis=1,how='all')
D_imputer = imputer_knn.fit(df, Ks="3,4,5,6,7")
df_knn = imputer_knn.transform(df)
df_knn.to_csv(sys.argv[1] + '_imputed',index=True, header=True,sep=",")

MAF = []
for i in range(0,df_knn.shape[1]):
	table = df_knn.iloc[:,i].value_counts()
	if max(table) <= 383*0.95:
		MAF.append(df_knn.columns[i])

if len(MAF) > 0:
	df2 = df_knn.loc[:,MAF]
	df2.to_csv(sys.argv[1] + '_imputed_MAF',index=True, header=True,sep=",")

num = df.isna().sum()
num75 = num[num<=8]
num90 = num[num<=2]
if len(num75) > 0:
	df75 = df2.loc[:,set(num75.index) & set(MAF)]
	df75.to_csv(sys.argv[1].split('_50per')[0] + '_75per_imputed_MAF',index=True, header=True,sep=",")
if len(num90) > 0:
	df90 = df2.loc[:,set(num90.index) & set(MAF)]
	df90.to_csv(sys.argv[1].split('_50per')[0] + '_90per_imputed_MAF',index=True, header=True,sep=",")


