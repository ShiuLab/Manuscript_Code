########################################################################################
# Calculate SHAP values of features to trait prediction, and feature interaction values.
# 
# Note that the feature interaction values were calculated for each instance, 
# i.e., accession in our study
# 
# Written by: Peipei Wang
# Modified by: Kenia Segura Abá
########################################################################################

import shap
import sys, os, argparse, time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datatable as dt
import joblib
import pickle
from sklearn.ensemble import RandomForestRegressor
from scipy import sparse

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(
		description='Calculate SHAP values, i.e., the contribution of each feature to the predictio of each instance; and feature interactions.')

	### Input arguments ###
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-df', help='Feature & class dataframe for ML, (example: example_binary.txt) ', required=True)
	
	# Optional
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-df2', help='Class data (if not in -df). Need to provide -y_name', default='')
	inp_group.add_argument('-sep', help='Deliminator', default='\t')
	inp_group.add_argument('-y_name', help='Name of column in df2 to predict', default='Y')
	inp_group.add_argument('-test', help='File with test intances', default='')
	inp_group.add_argument('-feat', help='File with list of features (from x) to include', default='all')
	inp_group.add_argument('-model', help='saved model', default='')
	inp_group.add_argument('-top', help='top number of features you want to show in the figures', default=20)
	inp_group.add_argument('-interaction', help='save interaction figures or not', default='n')
	inp_group.add_argument('-interaction_score', help='save interaction scores or not', default='n')

	# Output arguments
	out_group = parser.add_argument_group(title='OUTPUT OPTIONS')
	out_group.add_argument('-save', help='prefix for output files. CAUTION: will overwrite!', default='')

	# Default Hyperparameters for RF models
	params_group = parser.add_argument_group(title='DEFINE HYPERPARAMETERS')
	params_group.add_argument('-n_estimators', help='RF/GB parameter.', default=100)
	params_group.add_argument('-max_depth', help='RF/GB parameter. Grid Search [3, 5, 10]', default=5)
	params_group.add_argument('-max_features', help='RF/GB parameter. Grid Search [0.1, 0.25, 0.5, 0.75, sqrt, log2, None]', default='sqrt')
	params_group.add_argument('-n_jobs', '-p', help='Number of processors for parallel computing (max for HPCC = 14)', type=int, default=1)


	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()

	try:
		args.max_features = float(args.max_features)
	except:
		args.max_features = args.max_features

	args.n_estimators = int(args.n_estimators)
	args.max_depth = int(args.max_depth)

	####### Load Dataframe & Pre-process #######
	df = dt.fread(args.df,header=True,sep=args.sep)
	df = df.to_pandas()
	df = df.set_index(df.columns[0], drop=True)

	# If features  and class info are in separate files, merge them: 
	if args.df2 != '':
		start_dim = df.shape
		df_class = pd.read_csv(args.df2, sep=args.sep, index_col = 0)
		df = pd.concat([df_class[args.y_name], df], axis=1, join='inner')
		print('Merging the feature & class dataframes changed the dimensions from %s to %s (instance, features).' 
			% (str(start_dim), str(df.shape)))

	# Specify Y column - default = Class
	if args.y_name != 'Y':
		df = df.rename(columns = {args.y_name:'Y'})


	# Filter out features not in feat file given - default: keep all
	if args.feat != 'all':
		print('Using subset of features from: %s' % args.feat)
		with open(args.feat) as f:
			features = f.read().strip().splitlines()
			features = ['Y'] + features
		df = df.loc[:,features]

	# Check for Nas
	if df.isnull().values.any() == True:
		if args.drop_na.lower() in ['t', 'true']:
			start_dim = df.shape
			df = df.dropna(axis=0)
			print('Dropping rows with NA values changed the dimensions from %s to %s.' 
				% (str(start_dim), str(df.shape)))
		else:
			print(df.columns[df.isnull().any()].tolist())
			print('There are Na values in your dataframe.\n Impute them or add -drop_na True to remove rows with nas')
			quit()

	# Make sure Y is datatype numeric
	df['Y'] = pd.to_numeric(df['Y'], errors = 'raise')

	# Set up dataframe of test instances that the final models will be applied to
	if args.test !='':
		df_all = df.copy()
		print('Removing test instances to apply model on later...')
		with open(args.test) as test_file:
			test_instances = test_file.read().splitlines()
		try:
			test_df = df.loc[test_instances, :]
			df = df.drop(test_instances)
		except:
			test_instances = [int(x) for x in test_instances]
			test_df = df.loc[test_instances, :]
			df = df.drop(test_instances)
	else:
		test_df = 'None'
		test_instances = 'None'

	X_train = df.drop(['Y'], axis=1)
	y_train = df['Y']

	X_test = test_df.drop(['Y'], axis=1)
	y_test = test_df['Y']

	print("Model")
	if args.model != '':
			my_model = joblib.load(args.model)
	else:
		my_model = RandomForestRegressor(
			n_estimators=int(args.n_estimators),
			max_depth=args.max_depth,
			max_features=args.max_features,
			criterion='mse',
			random_state=42,
			n_jobs=args.n_jobs)
		my_model.fit(X_train, y_train)

	print("Calculating shap values")
	# calculate and save the SHAP value matrix
	explainer = shap.Explainer(my_model)
	# prediction on train
	shap_values = explainer(X_train)
	# plt.clf()
	# shap.plots.beeswarm(shap_values,max_display=21,show=False)
	# plt.savefig("SHAP_beeswarm_%s_training_top20.pdf"%args.save)
	values = pd.DataFrame(shap_values.values)
	values.index = X_train.index
	values.columns = X_train.columns
	#values.to_csv('SHAP_values_%s.txt'%args.save,sep='\t',header=True,index=True)
	
	# calculate the average absolute SHAP value of a column, and then rank the features based on this average abs SHAP value
	values_abs = abs(values)
	values_abs_mean = values_abs.mean(axis=0) # commented out by Kenia 01/30/2024
	values_abs_mean_sorted = values_abs_mean.sort_values(ascending = False) # commented out by Kenia 01/30/2024
	values_abs_mean_sorted.to_csv('SHAP_values_sorted_average_%s_training.txt'%args.save,sep='\t',header=True,index=True)

	# Added 01/30/2024: Calculate the median absolute SHAP value of a column, and then rank the features based on this median abs SHAP value
	# values_abs_med = values_abs.median(axis=0)
	# values_abs_med_sorted = values_abs_med.sort_values(ascending = False)

	# sort the column in the SHAP value matrix according to the average absolute values
	values_sorted = values[values_abs_mean_sorted.index]
	# values_sorted = values[values_abs_med_sorted.index] # Added by Kenia 01/30/2024
	values_sorted.to_csv('SHAP_values_sorted_%s_training.txt'%args.save,sep='\t',header=True,index=True)

	print("Reorder and recalculate shap")
	# reorder the X_training, and remake the shap_values. Otherwise in the following figure, the feature names would be wrong.
	X_train = X_train[values_abs_mean_sorted.index]
	# X_train = X_train[values_abs_med_sorted.index] # Added by Kenia 01/30/2024
	explainer = shap.Explainer(my_model)
	shap_values = explainer(X_train)

	# get rid of the features beyond top number of features
	values_top = values_sorted.iloc[:,0:int(args.top)]
	Data = pd.DataFrame(shap_values.data)
	Data_top = Data.iloc[:,0:int(args.top)]
	base_values_top = shap_values.base_values[0:int(args.top)]
	shap_values.values = np.array(values_top)
	shap_values.base_values = np.array(base_values_top)
	shap_values.data = np.array(Data_top)
	plt.clf()
	shap.plots.beeswarm(shap_values,max_display=int(args.top) + 1,show=False)
	plt.savefig("SHAP_beeswarm_%s_training_top%s_without_other_features.pdf"%(args.save,args.top))

	# min-max normalization of shap values
	values_normalized = values_top.copy(deep=True) # modified by Kenia 01/30/2024
	for col in values_top.columns.to_list():
		if sum(abs(values_top[col])) != 0:
			#if col in values_abs_mean_sorted.index.to_list()[0:int(args.top)]:
			max_shap = np.max(values_top[col])
			min_shap = np.min(values_top[col])
			max_abs = np.max([abs(max_shap),abs(min_shap)])
			values_normalized[col] = values_top[col] / max_abs

	# draw normalized SHAP values
	values_normalized_np = np.array(values_normalized)
	shap_values.values = values_normalized_np
	plt.clf()
	shap.plots.beeswarm(shap_values,max_display=int(args.top) + 1,show=False)
	plt.savefig("SHAP_beeswarm_%s_training_normalized_top%s.pdf"%(args.save,args.top))

	if args.interaction == 'y':
		# get the interaction between top features
		explainer_interaction = shap.TreeExplainer(my_model)
		shap_interaction = explainer_interaction.shap_values(X_train)
		for i in range(0,int(args.top)):
			for j in range(0,int(args.top)):
				if i != j:
					iloc = X_train.columns.get_loc(values_abs_mean_sorted.index.tolist()[i])
					jloc = X_train.columns.get_loc(values_abs_mean_sorted.index.tolist()[j])
					# iloc = X_train.columns.get_loc(values_abs_med_sorted.index.tolist()[i]) # Added by Kenia 01/30/2024
					# jloc = X_train.columns.get_loc(values_abs_med_sorted.index.tolist()[j]) # Added by Kenia 01/30/2024
					shap.dependence_plot(iloc, shap_interaction, X_train, display_features=X_train,interaction_index=jloc)
					plt.savefig("shap_interaction_%s_%s_%s_interaction.pdf"%(args.save,values_abs_mean_sorted.index.tolist()[i],values_abs_mean_sorted.index.tolist()[j]))
					# plt.savefig("shap_interaction_%s_%s_%s_interaction.pdf"%(args.save,values_abs_med_sorted.index.tolist()[i],values_abs_med_sorted.index.tolist()[j])) # Added by Kenia 01/30/2024
					plt.clf()

	if args.interaction_score == 'y':
		print("Calculating shap interaction scores")
		explainer_interaction = shap.TreeExplainer(my_model)
		print('**********DEBUG CORE DUMP**********')
		with open(f'shap_interaction_scores_{args.save}_explainer.pkl', 'wb') as file:
			pickle.dump(explainer_interaction, file)
		shap_interaction_values = explainer.shap_interaction_values(X_train)
		print('SHAP INTERACTION VALUES SUCCESSFULLY CALCULATED')
		print('***********************************')
		# save the ordered interaction
		dim = shap_interaction_values.shape
		shap_interaction_values_2d = np.reshape(np.ravel(shap_interaction_values), (dim[0], dim[1]*dim[2]))
		# make all pairs of features
		x2 = X_train[np.repeat(X_train.columns.tolist(), len(X_train.columns))]
		x2.columns =  [i+":"+str(j) for i in X_train.columns for j in X_train.columns]
		shap.summary_plot(shap_interaction_values_2d, x2, sort=True)
		plt.savefig("shap_interaction_scores_ordered_%s.pdf"%(args.save))
		plt.clf()
		shap.summary_plot(shap_interaction_values, X_train)
		plt.savefig("shap_interaction_scores_%s.pdf"%(args.save))
		plt.clf()
		shap_interaction_values_abs = abs(shap_interaction_values)
		shap_interaction_values_abs_mean = sum(shap_interaction_values_abs)/shap_interaction_values.shape[0]
		shap_interaction_values_abs_mean_df = pd.DataFrame(shap_interaction_values_abs_mean)
		shap_interaction_values_abs_mean_df.index = X_train.columns
		shap_interaction_values_abs_mean_df.columns = X_train.columns
		shap_interaction_values_abs_mean_df.to_csv("shap_interaction_scores_%s.txt"%(args.save),sep='\t',header=True,index=True)
		for j in range(0,shap_interaction_values.shape[0]):
			shap_interaction_values_j = pd.DataFrame(shap_interaction_values[j])
			shap_interaction_values_j.index = X_train.columns
			shap_interaction_values_j.columns = X_train.columns
			# shap_interaction_values_j.to_csv("shap_interaction_scores_%s_%s.txt"%(args.save,X_train.index.tolist()[j]),sep='\t',header=True,index=True) # Commented out by Kenia 02/02/2024
			# Added by Kenia 02/03/2024: Convert to scipy sparse matrix and save
			shap_interaction_values_j = sparse.csr_matrix(shap_interaction_values_j)
			sparse.save_npz("shap_interaction_scores_%s_%s.npz" % \
				(args.save, X_train.index.tolist()[j]), shap_interaction_values_j, compressed=True)

if __name__ == '__main__':
	main()


