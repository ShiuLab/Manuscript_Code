import sys, os, argparse, time
from irf import irf_utils
from irf.ensemble import RandomForestClassifierWithWeights
from irf.ensemble import RandomForestRegressorWithWeights
from sklearn.ensemble import RandomForestRegressor
import pandas as pd
import numpy as np
from datetime import datetime
start_total_time = time.time()
import datatable as dt
from sklearn import metrics
from scipy.stats.stats import pearsonr

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():

	########################
	### Parse Input Args ###
	########################
	parser = argparse.ArgumentParser(
		description='Machine learning regression pipeline using tools from Scikit-Learn. \
			See README.md for more information about the pipeline and preprocessing/post-analysis tools \
			available. All required packages \
			are available on MSUs HPCC: "export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH"',
		epilog='https://github.com/ShiuLab')
	
	
	### Input arguments ###
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-df', help='Feature & class dataframe for ML, (example: example_binary.txt) ', required=True)
	req_group.add_argument('-save_name', help='short name for the outputs', required=True)
	req_group.add_argument('-test', help='File with testing lines', required=True)
	# Optional
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-df2', help='Class data (if not in -df). Need to provide -y_name', default='')
	inp_group.add_argument('-sep', help='Deliminator', default=',')
	inp_group.add_argument('-y_name', help='Name of column in Y_file to predict', default='Y')
	inp_group.add_argument('-feat', help='File with list of features (from x) to include', default='all')
	inp_group.add_argument('-rf', help='classification or regression', default=RandomForestRegressorWithWeights)
	inp_group.add_argument('-rf_bootstrap', help='RandomForest model to fit to the bootstrap samples', default=None)
	inp_group.add_argument('-K', help='The number of iterations in iRF',type=int, default=7)
	inp_group.add_argument('-B', help='The number of bootstrap samples',type=int, default=10)
	inp_group.add_argument('-n_estimators', help='The number of trees in the random forest when computing weights',type=int, default=20)
	inp_group.add_argument('-propn_n_samples', help='The proportion of samples drawn for bootstrap',type=float, default=0.2)
	inp_group.add_argument('-bin_class_type', help='regression or classification, regression is 0, binary classification is 1',type=int, default=0)
	inp_group.add_argument('-max_depth', help='The built tree will never be deeper than `max_depth`',type=int, default=2)
	inp_group.add_argument('-num_splits', help='At each node, the maximum number of children to be added`',type=int, default=2)
	inp_group.add_argument('-n_estimators_bootstrap', help='The number of trees in the random forest when fitting to bootstrap samples',type=int, default=5)
	inp_group.add_argument('-save_model', help='Y,N', default='Y')
	inp_group.add_argument('-save_predictions', help='Y,N', default='Y')
	inp_group.add_argument('-save_importance', help='Y,N', default='Y')

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()

	####### Load Dataframe & Pre-process #######
	
	#df = pd.read_csv(args.df, sep=args.sep, index_col = 0)
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

	# assign function to parameter
	if args.rf == 'RandomForestClassifierWithWeights':
		args.rf = RandomForestClassifierWithWeights

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

	print("Snapshot of data being used:")
	print(df.iloc[:5, :5])

	n_features = len(list(df)) - 1

	X_train=df.drop(['Y'], axis=1)
	X_test=test_df.drop(['Y'],axis=1)
	y_train=df['Y']
	y_test=test_df['Y']
	# run the iRF to get the feature interactions
	# the output stability_score represents the proportion of times (out of B bootstrap samples) an interaction appears as an output of the RIT
	all_rf_weights, all_K_iter_rf_data, \
		all_rf_bootstrap_output, all_rit_bootstrap_output, \
		stability_score = irf_utils.run_iRF(X_train=X_train,
		X_test=X_test,
		y_train=y_train,
		y_test=y_test,
		K=int(args.K),                          # number of iteration
		rf = args.rf(n_estimators=int(args.n_estimators)),
		B=int(args.B),
		random_state_classifier=2018, # random seed
		propn_n_samples=float(args.propn_n_samples),
		bin_class_type=int(args.bin_class_type),
		M=20,
		max_depth=int(args.max_depth),
		noisy_split=False,
		num_splits=int(args.num_splits),
		n_estimators_bootstrap=int(args.n_estimators_bootstrap))
	
	# save the prediction and PCC 
	pred = {}
	for k in range(1,int(args.K)+1):
		y_pred = all_K_iter_rf_data['rf_iter%s'%k]['rf_validation_metrics']['prediction']
		y_pred = pd.DataFrame(y_pred)
		y_pred.index = y_test.index
		y_pred['True_y'] = y_test
		if k not in pred:
			pred[k] = {}
		pred[k] = y_pred
	
	out = open(args.save_name+'_performance.txt','w')
	out.write('K\tMse_loss\tPCC\n')
	for k in range(1,int(args.K)+1):
		mse_loss = metrics.mean_squared_error(pred[k]['True_y'], pred[k][0])
		PCC = pearsonr(pred[k]['True_y'], pred[k][0])[0]
		out.write('%s\t%s\t%s\n'%(k,mse_loss,PCC))
	out.close()
	# save the stability score
	colnames = df.drop(['Y'], axis=1).columns.tolist()
	stability_score_new = {}
	for keys in stability_score:
		cols = keys.split('_')
		keys_new = ''
		for col in cols:
			keys_new = keys_new + colnames[int(col)] + '__'
		keys_new = keys_new[0:-2]
		stability_score_new[keys_new] = round(stability_score[keys],3)
	
	stability_score_save = pd.DataFrame.from_dict(stability_score_new, orient='index')  ### make a dataframe from a dictionary
	stability_score_save = stability_score_save.sort_values(by = 0,ascending=False)
	stability_score_save.to_csv(args.save_name + '_stability_score.txt', header=True, index=True,sep='\t')
	
	# save the weights for each iteration
	all_rf_weights_save = pd.DataFrame.from_dict(all_rf_weights, orient='columns')
	all_rf_weights_save = all_rf_weights_save.iloc[:,1:]
	all_rf_weights_save.index = df.columns.tolist()[1:]
	all_rf_weights_save = all_rf_weights_save.sort_values(by = 'rf_weight%s'%args.K,ascending=False)
	all_rf_weights_save.to_csv(args.save_name + '_weights.txt', header=True, index=True,sep='\t')
	
	# save the results of each iteration
	if args.save_model == "Y":
		all_K_iter_rf_data_save = pd.DataFrame.from_dict(all_K_iter_rf_data, orient='columns')
		import pickle
		f = open(args.save_name + '_data.pkl',"wb")
		pickle.dump(all_K_iter_rf_data_save,f)
		f.close()

	# save predictions of the each iteration
	Prediction = pd.DataFrame(index=y_test.index)
	if args.save_predictions == "Y":
		for k in range(1,int(args.K)+1):
			predictions = pd.DataFrame.from_dict(all_K_iter_rf_data['rf_iter%s'%k]['rf_validation_metrics']['prediction'])
			predictions.index = y_test.index
			Prediction['rf_iter%s'%k] = predictions
	Prediction.to_csv(args.save_name + '_predictions.txt', header=True, index=True,sep='\t')
	
	# save feature importance of the each iteration
	Importance = pd.DataFrame(index=X_test.columns)
	if args.save_importance == "Y":
		for k in range(1,int(args.K)+1):
			importance = pd.DataFrame.from_dict(all_K_iter_rf_data['rf_iter%s'%k]['feature_importances'])
			importance.index = X_test.columns
			importance_std = pd.DataFrame.from_dict(all_K_iter_rf_data['rf_iter%s'%k]['feature_importances_std'])
			importance_std.index = X_test.columns
			feature_importances_rank_idx = pd.DataFrame.from_dict(all_K_iter_rf_data['rf_iter%s'%k]['feature_importances_rank_idx'])
			feature_importances_rank_idx.index = X_test.columns
			Importance['imp%s'%k] = importance
			Importance['imp%s_std'%k] = importance_std
			Importance['imp%s_rank'%k] = feature_importances_rank_idx
	Importance.to_csv(args.save_name + '_importance.txt', header=True, index=True,sep='\t')
	
	
if __name__ == '__main__':
	main()


