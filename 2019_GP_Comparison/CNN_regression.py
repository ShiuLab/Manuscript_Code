from __future__ import absolute_import, division, print_function, unicode_literals
import os
import numpy as np
import pandas as pd
import timeit
import argparse
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow.keras import datasets, layers, models, optimizers
from tensorflow.keras.callbacks import EarlyStopping
import tensorflow.keras.backend as K
from datetime import datetime
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV
from keras.wrappers.scikit_learn import KerasClassifier
from sklearn import preprocessing, metrics
from sklearn.preprocessing import OneHotEncoder, LabelEncoder, MinMaxScaler
from sklearn.utils import class_weight
import ANN_Functions as ANN

start_time = timeit.default_timer()

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn


def main():
	########################
	### Parse Input Args ###
	########################
	parser = argparse.ArgumentParser(
		description='Predicting Traits Using Convolutional Neural Networks. \
			See README.md for more information about the pipeline usage.\n'
			'Written by: Emily Bolger\nModified by: Christina Azodi',
		epilog='https://github.com/ShiuLab')

	## Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-x', help='Feature data. Format options available: '
		' 1) Matrix with one row for every instance. 2) Directory with one '
		'matrix for every instance with file name matching instance name '
		'provided in -y', required=True)
	req_group.add_argument('-y', help='Matrix with label to predict',
		required=True)
	req_group.add_argument('-test', help='List of instances to use as test set',
		required=True)
	req_group.add_argument('-save', help='Name for Output File', required=True)


	## Optional 
	inp_group = parser.add_argument_group(title='OPTIONAL INPUT')
	inp_group.add_argument('-feat', help="List of column names in -x to use",
		default=False)
	inp_group.add_argument('-x_sort', help='Method to sort feature (-x) data '
		'by column. Options: (False, alpha, file_with_order, cluster)',
		 default=False)
	inp_group.add_argument('-shape', help='Dimension of -x (e.g. for input '
		'with 4 rows and 6 columns: -shape 4,6)', default='default')
	inp_group.add_argument('-onehot', help='T/F to convert 1xX data into matrix'
		' by one-hot encoding', default='F')
	inp_group.add_argument('-onehot_order', help='Order for 1-hot y axis. For '
		'example: -onehot_order 1,0,-1)', default=False)
	inp_group.add_argument('-y_name', help='Col name to predict', default='Y')
	inp_group.add_argument('-sep', help='Deliminator for X & Y', default='\t')


	# How to run CNN
	inp_group.add_argument('-run', help='T/F to run the final model. If F, will'
		'only run the grid search if -gs T', default='f')
	inp_group.add_argument('-n', help='Num replicates of model... Different'
		'train/validation split each replicate', type=int, default=10)
	inp_group.add_argument('-n_jobs', '-p', help='Number of processors for '
		'parallel computing (max for HPCC = 14)', type=int, default=1)
	inp_group.add_argument('-cv', help='Number of cross validation folds',
		type=int, default=5)

	# Parameter Selection
	inp_group.add_argument('-params', help='How to select parameters. Options: '
		'grid search (gs), default, from XXX_GridSearch.txt (provide path).',
		default='default')
	inp_group.add_argument('-gs_reps', help='Number of combinations of '
		'parameters to test in the grid search', type=int, default=100)

	# Default CNN parameters
	inp_group.add_argument('-cnn_type', help='CNN architecture. Options: '
		'(simple, DeepGS)', default="simple")
	inp_group.add_argument('-filters', help='Number of kernels/filters in each '
		'CNN layer', type=int, default=32)
	inp_group.add_argument('-kernel_l', help='Length of kernel (height '
		'defaults to the full height of the dataset)', type=int, default=16)
	inp_group.add_argument('-stride_len', help='Stride of Convolution kernels '
		'(width defaults to 1)', type=int, default=1)
	inp_group.add_argument('-activation', help='Activation function in all but '
		'last dense layer, which is set to linear', type=str, default='relu')
	inp_group.add_argument('-pool_size', help='Size of max pooling layer filter '
		'(first number only, second defaults to 1)', type=int, default=8)
	inp_group.add_argument('-optimizer', help='Optimization function to use)',
		type=str, default='Adam')
	inp_group.add_argument('-dropout', help='Value for Dropout Rate',
		type=float, default=0.5)
	inp_group.add_argument('-l2', help='Value for L2 regularization',
		type=float, default=0.2)
	inp_group.add_argument('-learn_rate', help='Value for Learning Rate',
		type=float, default=0.01)
	inp_group.add_argument('-clip_value', help='Clip Value', type=float,
		default=0.5)
	inp_group.add_argument('-patience', help='Patience for Early Stopping',
		type=int, default=10)
	inp_group.add_argument('-min_delta', help='Minimum Delta Value for Early '
		'Stopping', type=float, default=0)
	inp_group.add_argument('-num_epochs', help='Max number of Epochs',
		type=int, default=1000)
	inp_group.add_argument('-n_channels', help='Num channels', type=int,
		default=1)

	# Argument parsing
	args = parser.parse_args()

	if args.shape == 'default':
		tmp = pd.read_csv(args.x, sep=args.sep, index_col=0)
		shape_r, shape_c = 1, tmp.shape[1]
	else:
		shape_r, shape_c = args.shape.strip().split(',')
		shape_r = int(shape_r)
		shape_c = int(shape_c)


	########################
	### Parse Input Data ###
	########################
	print("\n***** Loading Data ******\n")

	# Step 1: Read in x file, if feat file given only keep those features
	if os.path.isfile(args.x):
		x = pd.read_csv(args.x, sep=args.sep, index_col=0)
		x.index = x.index.astype('str')
		instance_order = list(x.index.values)
		with open(args.test) as test_file:
			test_instances = test_file.read().splitlines()
			test_instances = [str(i) for i in test_instances]
			train_val_instances = list(set(instance_order) - set(test_instances))
			test_index = [x.index.get_loc(i) for i in test_instances]
			train_val_index = [x.index.get_loc(i) for i in train_val_instances]
		if args.feat:
			with open(args.feat) as f:
				features = f.read().strip().splitlines()
			x = x.loc[:, features]
	elif os.path.isdir(args.x):
		x = ANN.fun.Image2Features(args.x, shape_r, shape_c)
	n_instances = x.shape[0]
	n_feats = x.shape[1]
	print("Total number of instances: %i" % n_instances)
	print("Number of features used: %i" % n_feats)

	# Step 2: Sort x data
	if args.x_sort == 'alpha':
		print('Sorting feature data by column alpha numerically...')
		x = x.reindex(sorted(x.columns), axis=1)
	elif args.x_sort == 'cluster':
		print('Sorting feature data by column using clustering...')
		print('\n\nNOT IMPLEMENTED YET... PROGRESSING WITHOUT SORTING...\n\n')
	else:
		if not args.x_sort:
			print('Using -x in the order provided in -x or in -feat')
		else:
			with open(args.x_sort) as order:
				order_list = order.read().strip().splitlines()
			x = x.loc[:, order_list]
	print('\nSnapshot of input feature data:')
	print(x.head())

	# Step 3: One-hot-encode X if required
	if args.onehot.lower() in ['t', 'true']:
		x_1hot_list = []
		x = x.round(0)
		labels = pd.unique(x.values.ravel())
		ohe = preprocessing.OneHotEncoder(categories='auto', sparse=False)
		for i in range(len(x)):
			x_row = np.array(x.iloc[i, ]).reshape(n_feats, -1)
			oh_matrix = ohe.fit_transform(x_row)
			if oh_matrix.shape[1] < len(labels):
				labels_present = pd.unique(x_row.ravel())
				missing = list(set(labels) - set(labels_present))
				print("Instance in row %i is has no '%s', so adding by hand..."
					% (i, missing))
				x_row = np.append(x_row, np.array([missing]), axis=0)
				oh_matrix = ohe.fit_transform(x_row)
				oh_matrix = oh_matrix[:-1, :]
			x_1hot_list.append(oh_matrix)
		x = np.swapaxes(np.array(x_1hot_list), 1, 2)

	data_height = x.shape[1]
	data_width = x.shape[2]
	x = x.reshape((n_instances, data_height, data_width, args.n_channels))
	print("\nShape of feature data used for training/validation/testing:")
	print(x.shape)

	print("\nSnapshot of feature data for first instance in data set:")
	print(x[0, :, :, 0])

	# Step 4: Read in Y data and make sure sorted as in -x
	y = pd.read_csv(args.y, sep=args.sep, index_col=0)
	y.index = y.index.astype('str')
	y = y[[args.y_name]]

	print("\nShape of Label Data:")
	print(y.shape)

	# Step 5: Remove testing data
	x_test = x[test_index, :, :, :]
	x_train = x[train_val_index, :, :, :]
	y_test = y.ix[test_instances]
	y_train = y.ix[train_val_instances]


	################################
	### Define CNN architectures ###
	################################

	def tfp_pearson(y_true, y_pred):
		return tfp.stats.correlation(y_pred, y_true, event_axis=None)

	def make_cnn_model(cnn_type=args.cnn_type, learn_rate=args.learn_rate,
		filters=args.filters, pool_size=args.pool_size,
		kernel_l=args.kernel_l, kernel_h=data_height,
		activation=args.activation, optimizer=args.optimizer, units=1):

		if optimizer.lower() == 'adam':
			opt = tf.keras.optimizers.Adam(lr=learn_rate,
				clipvalue=args.clip_value)
		elif optimizer.lower() == 'nadam':
			opt = tf.keras.optimizers.Nadam(lr=learn_rate,
				clipvalue=args.clip_value)
		elif optimizer.lower() == 'rmsprop':
			opt = tf.keras.optimizers.RMSprop(lr=learn_rate,
				clipvalue=args.clip_value)

		if cnn_type.lower() == 'simple':
			K.clear_session()
			model = models.Sequential()
			model.add(layers.Conv2D(filters=filters,
				kernel_size=tuple([kernel_h, kernel_l]),
				kernel_regularizer=tf.keras.regularizers.l2(args.l2),
				strides=tuple([args.stride_len, 1]),
				activation=activation,
				kernel_initializer='glorot_normal',
				input_shape=(data_height, data_width, args.n_channels)))
			model.add(layers.MaxPooling2D(pool_size=tuple([1, pool_size])))
			model.add(layers.Flatten())
			model.add(layers.Dropout(args.dropout))
			model.add(layers.Dense(24, activation=activation))
			model.add(layers.BatchNormalization())
			model.add(layers.Dense(units=units, activation='linear'))
			model.compile(optimizer=opt, loss='mean_squared_error')

		elif cnn_type.lower() == 'deepgs':
			K.clear_session()
			model = models.Sequential()
			model.add(layers.Conv2D(filters=filters,
				kernel_size=tuple([kernel_h, kernel_l]),
				strides=tuple([args.stride_len, 1]),
				activation=activation,
				kernel_initializer='glorot_normal',
				input_shape=(data_height, data_width, args.n_channels)))
			model.add(layers.MaxPooling2D(pool_size=tuple([1, pool_size])))
			model.add(layers.Conv2D(filters=filters,
				kernel_size=tuple([1, kernel_l]),
				strides=tuple([args.stride_len, 1]),
				activation=activation))
			model.add(layers.MaxPooling2D(pool_size=tuple([1, pool_size])))
			model.add(layers.Dropout(args.dropout))
			model.add(layers.Flatten())
			model.add(layers.Dense(units=24, activation=activation))
			model.add(layers.BatchNormalization())
			model.add(layers.Dropout(args.dropout))
			model.add(layers.Dense(units=units, activation='linear'))
			model.compile(optimizer=opt, loss='mean_squared_error')
		return model

	####################
	### Grid Search  ###
	####################

	if args.params.lower() == 'gs':
		print('\n***** Starting Random Search with %i reps using %i testing '
			'instances and %i fold cross-validation *****\n' % (
			args.gs_reps, x_train.shape[0], args.cv))
		scoring = {'neg_mse': 'neg_mean_squared_error', 'exp_var': 'explained_variance'}
		param_grid = dict(learn_rate=[1, 0.1, 0.01, 0.001, 0.0001, 0.00001],
			filters=[8, 16, 32],
			kernel_l=[8, 16, 32],
			pool_size=[4, 8, 16],
			activation=["relu", "selu", "elu"],
			optimizer=['RMSprop', 'Adam', 'nadam'],
			cnn_type=['simple', 'deepgs'])
		model = KerasClassifier(build_fn=make_cnn_model,
			batch_size=100,
			epochs=50,
			verbose=1)
		rand_search = RandomizedSearchCV(estimator=model,
			param_distributions=param_grid,
			cv=args.cv,
			n_iter=args.gs_reps,
			n_jobs=args.n_jobs,
			verbose=1,
			scoring=scoring,
			refit='neg_mse')
		gs_result = rand_search.fit(x_train, y_train)
		gs_result_df = pd.DataFrame.from_dict(gs_result.cv_results_)

		print("Saving Grid Search Results....")
		print(gs_result_df.head())
		with open(args.save + "_GridSearch.txt", 'a') as out_gs:
			gs_result_df.to_csv(out_gs, header=out_gs.tell() == 0, sep='\t')

	########################
	### Run Final Models ###
	########################

	if args.run.lower() in ['t', 'true']:

		# Step 1: Define the parameters from the Grid Search or use default
		if args.params.lower() != 'default':
			if args.params.lower() != 'gs':
				gs_result_df = pd.read_csv(args.params, sep='\t')
				gs_result_df.fillna(0, inplace=True)

			gs_mean = gs_result_df.groupby(['param_filters', 'param_optimizer',
				'param_learn_rate', 'param_kernel_l', 'param_pool_size',
				'param_cnn_type', 'param_activation']).agg({'mean_test_score':
				'mean', 'std_test_score': 'mean'}).reset_index()

			gs_mean = gs_mean.sort_values(by='mean_test_score', ascending=False)
			print('\nSnapshot of grid search results:')
			print(gs_mean.head())

			args.cnn_type = gs_mean['param_cnn_type'].iloc[0]
			args.pool_size = int(gs_mean['param_pool_size'].iloc[0])
			args.learn_rate = float(gs_mean['param_learn_rate'].iloc[0])
			args.kernel_l = int(gs_mean['param_kernel_l'].iloc[0])
			args.filters = int(gs_mean['param_filters'].iloc[0])
			args.activation = gs_mean['param_activation'].iloc[0]
			args.optimizer = gs_mean['param_optimizer'].iloc[0]

		print('\n***** Running CNN models ******')
		print('CNN Architecture: %s\nOptimizer: %s\nActivation function:'
			' %s\nLearning Rate: %f\nNumber of kernels: '
			'%i\nKernel shape: [%i, %i]\nPooling Size: [%i, 1]\n' % (
				args.cnn_type, args.optimizer, args.activation, args.learn_rate,
				args.filters, args.kernel_l, data_height, args.pool_size))

		for i in range(args.n):
			print('Rep %i of %i' % (i, args.n))
			run = True

			while run:
				# Step 2: Creating CNN model using Tensorflow
				model = make_cnn_model(cnn_type=args.cnn_type,
					learn_rate=args.learn_rate,
					optimizer=args.optimizer,
					filters=args.filters,
					pool_size=args.pool_size,
					kernel_l=args.kernel_l,
					kernel_h=data_height,
					activation=args.activation,
					units=1)
				#print(model.summary())

				# Step 3: Split training into training2 and validation
				x_train2, x_val, y_train2, y_val = train_test_split(x_train,
					y_train, test_size=0.2)
				print('Train on %i, validate on %i, test on %i' % (
					x_train2.shape[0], x_val.shape[0], x_test.shape[0]))

				# Step 4: Define optimizer and early stopping criteria & train
				model.compile(optimizer=args.optimizer,
					loss='mean_squared_error', metrics=[tfp_pearson])

				earlystop_callback = EarlyStopping(monitor='val_loss',
					mode='min',
					min_delta=args.min_delta,
					patience=args.patience,
					restore_best_weights=True,
					verbose=1)

				model.fit(x_train2, y_train2,
					batch_size=100,
					epochs=args.num_epochs,
					verbose=1,
					callbacks=[earlystop_callback],
					validation_data=(x_val, y_val))

				# Step 5: Apply best model to train, val, test, and report results
				train_mse, train_pcc = model.evaluate(x_train2, y_train2)
				val_mse, val_pcc = model.evaluate(x_val, y_val)
				if val_pcc > 0:
					run = False
				else:
					print('\nPCC was negative on valid data.. retraining...')

			test_mse, test_pcc = model.evaluate(x_test, y_test)
			if np.isnan(test_pcc):
				# Still don't know why this happens, but this fixes it...
				print('Recalculating PCC using Numpy...')
				pred = model.predict(x_test).tolist()
				pred2 = [i for sublist in pred for i in sublist]
				test_pcc = np.corrcoef(pred2, y_test[args.y_name].values)[0, 1]

			print('PCC: train, val, and test: %3f, %3f, %3f' % (train_pcc,
				val_pcc, test_pcc))

			if not os.path.isfile('RESULTS.txt'):
				out = open('RESULTS.txt', 'w')
				out.write('ID\tX\tY\ty_name\ttest_set\t'
					'CNN_Type\tLearn_Rate\tMin_Delta\tPatience\tActivation\t'
					'Optimizer\tKernel_num\tKernel_len\tPooling_Size\tDropout'
					'\tTrain_mse\tTrain_PCC\tVal_mse\tVal_PCC\tTest_mse\t'
					'Test_PCC\n')
				out.close()

			out = open('RESULTS.txt', "a")
			out.write('%s\t%s\t%s\t%s\t%s\t'
				'%s\t%f\t%f\t%i\t%s\t%s\t'
				'%i\t%i\t%i\t%f\t%f\t'
				'%f\t%f\t%f\t%f\t%f\n' % (
				args.save, args.x, args.y, args.y_name, args.test, args.cnn_type,
				args.learn_rate, args.min_delta, args.patience, args.activation,
				args.optimizer, args.filters, args.kernel_l, args.pool_size,
				args.dropout, train_mse, train_pcc, val_mse, val_pcc, test_mse,
				test_pcc))
			out.close()

		print('\nDone!')

if __name__ == '__main__':
	main()
