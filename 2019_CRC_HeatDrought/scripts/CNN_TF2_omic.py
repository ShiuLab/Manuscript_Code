'''
module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load CUDA/10.0.130 cuDNN/7.5.0.56-CUDA-10.0.130
module load Python/3.6.4
source /mnt/home/azodichr/tf2env/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2
cd ~/01_CombinedStress/Prasch_HD/08_10clusters/05_OmicIntegration/04_CNN/
python ~/GitHub/Combined_Stress_Response/TF_CNN.py -x test_x.npy -y test_y.npy

'''

from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import pandas as pd
import tensorflow as tf
import timeit, os
import argparse
import copy
from tensorflow.keras import datasets, layers, models, optimizers
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.model_selection import RandomizedSearchCV
from keras.wrappers.scikit_learn import KerasClassifier
from keras.constraints import max_norm
import tensorflow.keras.backend as K
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.utils import resample
from sklearn import preprocessing, metrics
from sklearn.utils import class_weight
from sklearn.model_selection import ParameterSampler
start_time = timeit.default_timer()
from tf_explain.callbacks.activations_visualization import ActivationsVisualizationCallback

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn


def main():

	########################
	### Parse Input Args ###
	########################
	parser = argparse.ArgumentParser(
		description='CNN classification code implemented using TensorFlow v2.0',
		epilog='https://github.com/azodichr')

	parser.add_argument('-x', help='Feature numpy dataset', required=True)
	parser.add_argument('-y', help='Class/Y numpy dataset', required=True)
	parser.add_argument('-run', help='T/F to run final models', default='t')
	parser.add_argument('-splits', help='Values for train/val/test',
		default='70,10,20')
	parser.add_argument('-y_name', help='Phenotype Trait')
	parser.add_argument('-f', help='Function: gs, run, full', default='full')
	parser.add_argument('-save', help='Name for Output File', default='test')
	parser.add_argument('-balance', help='t/f to downsample so balance classes',
		default='t')
	parser.add_argument('-n_channels', help='Number of channels', default=1,
		type=int)
	parser.add_argument('-cv', help='Number of cross validation folds',
		type=int, default=5)
	parser.add_argument('-n_jobs', '-p', help='Number of processors for '
		'parallel computing (max for HPCC = 14)', type=int, default=1)
	parser.add_argument('-save_model', help='T/F if want to save final models',
		type=str, default='f')
	parser.add_argument('-tag', help='Identifier String to add to RESULTS',
		type=str, default='cnn')
	parser.add_argument('-save_detailed', help='T/F Save detailed model performance',
		type=str, default='f')
	parser.add_argument('-original_df', help='DF fed into input_converter.py',
		type=str, default='')
	parser.add_argument('-imp_m', help='T/F to calculate importance of each motif',
		type=str, default='f')
	parser.add_argument('-imp_k', help='T/F to calculate importance of each kernel',
		type=str, default='f')

	# Default Hyperparameters
	parser.add_argument('-params', help='Output from -f gs (i.e. '
		'SAVE_GridSearch.txt)', default='default')
	parser.add_argument('-actfun', help='Activation function. (relu, sigmoid)',
		default='relu')
	parser.add_argument('-learn_rate', help='Learning Rate', default=0.01,
		type=float)
	parser.add_argument('-dropout', help='Dropout rate', default=0.25,
		type=float)
	parser.add_argument('-l2', help='Shrinkage parameter for L2 regularization',
	 default=0.25, type=float)
	parser.add_argument('-filters', help='Number of Kernels/filters',
		default=8, type=int)
	parser.add_argument('-optimizer', help='Optimization function to use)',
		type=str, default='Adam')
	parser.add_argument('-dense', help='Number of nodes in dense layer',
		type=int, default=16)
	parser.add_argument('-activation', help='Activation function in all but '
		'last dense layer, which is set to linear', type=str, default='relu')
	parser.add_argument('-n_reps', '-n', help='Number of replicates (unique '
		'validation set/starting weights for each)', default=100, type=int)
	parser.add_argument('-clip_value', help='Clip Value', type=float,
		default=0.5)
	parser.add_argument('-patience', help='Patience for Early Stopping',
		type=int, default=5)
	parser.add_argument('-min_delta', help='Minimum Delta Value for Early '
		'Stopping', type=float, default=0)

	# Grid Search reps/space
	parser.add_argument('-gs_reps', '-gs_n', help='Number of Grid Search Reps'
		'(will append results if SAVE_GridSearch.csv exists)', type=int,
		default=10)
	parser.add_argument('-actfun_gs', help='Activation functions for Grid '
		'Search', nargs='*', default=['relu', 'selu', 'elu'])
	parser.add_argument('-dropout_gs', help='Dropout rates for Grid Search',
		nargs='*', type=float, default=[0.0, 0.1, 0.25])
	parser.add_argument('-l2_gs', help='Shrinkage parameters for L2 for Grid '
		'Search', nargs='*', type=float, default=[0.01, 0.1, 0.25])
	parser.add_argument('-lrate_gs', help='Learning Rate', nargs='*', 
		type=float, default=[0.1, 0.01, 0.001, 0.0001])
	parser.add_argument('-kernels_gs', help='Number of Kernels for Grid Search',
	 default=[4, 8, 16, 24], type=int)

	args = parser.parse_args()
	k_height = 'tmp'
	args.k_len = 'tmp'

	def downsample(x, y):
		unique, counts = np.unique(y_all, return_counts=True)
		smaller_index = list(counts).index(min(counts))
		bigger_index = list(counts).index(max(counts))

		i_smaller = np.where(y_all == unique[smaller_index])[0]
		i_bigger = np.where(y_all == unique[bigger_index])[0]
		downsample_n = len(i_smaller)
		i_bigger_downsampled = np.random.choice(i_bigger, size=downsample_n,
			replace=False)

		i_keep = list(i_smaller) + list(i_bigger_downsampled)
		y = y_all[i_keep]
		x = x_all[i_keep]

		return x, y

	def make_cnn_model(learn_rate=args.learn_rate, filters=args.filters,
		dropout=args.dropout, dense=args.dense,
		l2=args.l2, activation=args.activation, optimizer=args.optimizer,
		units=1):

		if optimizer.lower() == 'adam':
			opt = tf.keras.optimizers.Adam(lr=learn_rate,
				clipvalue=args.clip_value)
		elif optimizer.lower() == 'nadam':
			opt = tf.keras.optimizers.Nadam(lr=learn_rate,
				clipvalue=args.clip_value)
		elif optimizer.lower() == 'rmsprop':
			opt = tf.keras.optimizers.RMSprop(lr=learn_rate,
				clipvalue=args.clip_value)
		elif optimizer.lower() == 'sgdm':
			opt = tf.keras.optimizers.SGD(lr=learn_rate, decay=1e-6,
				clipvalue=args.clip_value, momentum=0.9, nesterov=True)

		conv2d_layer = layers.Conv2D(filters=filters,
			kernel_size=tuple([k_height, 1]),
			kernel_regularizer=tf.keras.regularizers.l2(l2),
			activation=activation,
			kernel_initializer='glorot_normal',
			input_shape=(n_rows, n_columns, args.n_channels),
			name='conv2d_layer')
		K.clear_session()
		model = models.Sequential()
		model.add(conv2d_layer)
		model.add(layers.Flatten())
		model.add(layers.Dense(dense, activation=activation))
		model.add(layers.Dropout(dropout))
		model.add(layers.Dense(units=1, activation='sigmoid'))
		model.compile(optimizer=opt, loss='binary_crossentropy')

		return model, conv2d_layer


	##########################
	### Data preprocessing ###
	##########################
	x_all = np.load(args.x)
	y_all = np.load(args.y)
	x_all = x_all.reshape(x_all.shape + (args.n_channels,))

	if args.balance.lower() in ['t', 'true']:
		x, y = downsample(x_all, y_all)

		print('Y shape (down-sampled): %s' % str(y.shape))
		print('X shape (down-sampled): %s' % str(x.shape))
	else:
		y = y_all
		x = x_all

	print("\nSnapshot of feature data for first instance in data set:")
	print(x[0, :, 0:5, 0])
	n = y.shape[0]
	n_rows = x.shape[1]
	n_columns = x.shape[2]

	k_height = x.shape[1]
	args.k_len = 1
	print('Kernel dimensions: ', k_height, args.k_len)

	###################
	### Grid Search ###
	###################

	if args.params.lower() == 'gs':
		print('\n***** Starting Random Search with %i reps using %i testing '
			'instances and %i fold cross-validation *****\n' % (
			args.gs_reps, x.shape[0], args.cv))
		scoring = {'acc': 'accuracy', 'f1': 'f1'}
		param_grid = dict(learn_rate=[0.1, 0.01, 0.001],
			filters=[8, 16],
			dense=[8, 16, 32],
			l2=[0.1, 0.25],  #, 0.5],
			dropout=[0.1, 0.25],  #, 0.5],
			activation=["relu"],  #, 'selu', 'elu'],
			optimizer=['RMSprop', 'Adam', 'nadam'])
		model, conv2d_layer = KerasClassifier(build_fn=make_cnn_model,
			batch_size=100,
			epochs=30,
			verbose=0)
		rand_search = RandomizedSearchCV(estimator=model,
			param_distributions=param_grid,
			cv=args.cv,
			n_iter=args.gs_reps,
			n_jobs=args.n_jobs,
			scoring=scoring,
			refit='acc',
			verbose=0)
		gs_result = rand_search.fit(x, y)
		gs_result_df = pd.DataFrame.from_dict(gs_result.cv_results_)

		print("Saving Grid Search Results....")
		print(gs_result_df.head())
		with open(args.save + "_GridSearch.txt", 'a') as out_gs:
			gs_result_df.to_csv(out_gs, header=out_gs.tell() == 0, sep='\t')

	print('\n\n Grid Search results saved to: %s_GridSearch.txt\n' % args.save)

	################
	### Run final model 
	################

	if args.run.lower() in ['t', 'true']:
		print('####### Running Final Model(s) ###########')

		# Step 1: Define the parameters from the Grid Search or use default
		if args.params.lower() != 'default':
			if args.params.lower() != 'gs':
				gs_result_df = pd.read_csv(args.params, sep='\t')
				gs_result_df.fillna(0, inplace=True)

			gs_mean = gs_result_df.groupby(['param_filters', 'param_optimizer',
				'param_learn_rate', 'param_dropout', 'param_l2',
				'param_dense', 'param_activation']).agg({'mean_test_acc':
				'mean', 'std_test_acc': 'mean', 'mean_fit_time': 'count'
				}).reset_index()

			print('Parameter Search Coverage: \nMin: %i\nMean: %3f\nMax:%i' %
				(gs_mean['mean_fit_time'].min(),
					gs_mean['mean_fit_time'].mean(),
					gs_mean['mean_fit_time'].max()))

			if gs_mean['mean_fit_time'].min() == 1:
				print('Dropping parameter combinations with < 2 replicates...')
				gs_mean = gs_mean[gs_mean['mean_fit_time'] >= 2]

			gs_mean = gs_mean.sort_values(by='mean_test_acc', ascending=False)
			print('\nSnapshot of grid search results:')
			print(gs_mean.head())

			args.learn_rate = float(gs_mean['param_learn_rate'].iloc[0])
			args.l2 = float(gs_mean['param_l2'].iloc[0])
			args.dropout = float(gs_mean['param_dropout'].iloc[0])
			args.filters = int(gs_mean['param_filters'].iloc[0])
			args.dense = int(gs_mean['param_dense'].iloc[0])
			args.activation = gs_mean['param_activation'].iloc[0]
			args.optimizer = gs_mean['param_optimizer'].iloc[0]

		print('\n***** Running CNN models ******')
		print('Optimizer: %s\nActivation function:'
			' %s\nLearning Rate: %4f\nNumber of kernels: '
			'%i\nL2: %4f\nDropout: %4f\nDense nodes: %s\n' % (
				args.optimizer, args.activation, args.learn_rate,
				args.filters, args.l2, args.dropout, args.dense))

		final_results = pd.DataFrame()
		motif_imps = pd.DataFrame()
		kern_imp = []

		for n in range(args.n_reps):
			print("\nReplicate %i/%i" % (n, args.n_reps))
			x, y = downsample(x_all, y_all)
			print(x.shape)

			model, conv2d_layer = make_cnn_model(learn_rate=args.learn_rate,
				optimizer='sgdm',
				filters=args.filters,
				dense=args.dense,
				l2=args.l2,
				dropout=args.dropout,
				activation=args.activation)
			#print(model.summary())

			# Step 3: Split training into training2 and validation
			x_train, x_test, y_train, y_test = train_test_split(x,
				y, stratify=y, test_size=0.1)
			x_train, x_val, y_train, y_val = train_test_split(x_train,
				y_train, stratify=y_train, test_size=0.111)
			print('Train on %i, validate on %i, test on %i' % (
				x_train.shape[0], x_val.shape[0], x_test.shape[0]))

			# Step 4: Define optimizer and early stopping criteria & train
			model.compile(optimizer=args.optimizer,
				loss='binary_crossentropy', metrics=['accuracy'])

			earlystop_callback = EarlyStopping(monitor='val_loss',
				mode='min',
				min_delta=args.min_delta,
				patience=args.patience,
				restore_best_weights=True,
				verbose=0)

			model.fit(x_train, y_train,
				batch_size=50,
				epochs=1000,
				verbose=0,
				callbacks=[earlystop_callback],
				validation_data=(x_val, y_val))

			train_loss, train_acc = model.evaluate(x_train, y_train)
			val_loss, val_acc = model.evaluate(x_val, y_val)
			test_loss, test_acc = model.evaluate(x_test, y_test)

			val_yhat = model.predict(x_val)
			max_f1 = 0
			best_thresh = 0
			for thr in np.arange(0.01, 1, 0.01):
				thr_pred = val_yhat.copy()
				thr_pred[thr_pred >= thr] = 1
				thr_pred[thr_pred < thr] = 0
				if sum(thr_pred) > 1:  # Eliminates cases where all predictions are negative and the f1 and auROC are undefined
					f1 = f1_score(y_val, thr_pred, pos_label=1)  # Returns F1 for positive class
					if f1 >= max_f1:
						max_f1 = f1
						best_thresh = thr
			print('Threshold for F1 measure: %3f' % best_thresh)

			# Calculate AUC-ROC and F-measure from train, val, and test.
			yhat_train = model.predict(x_train)
			train_auroc = roc_auc_score(y_train, yhat_train)
			yhat_train[yhat_train >= best_thresh] = 1
			yhat_train[yhat_train < best_thresh] = 0
			train_f1 = f1_score(y_train, yhat_train, pos_label=1)

			yhat_val = model.predict(x_val)
			val_auroc = roc_auc_score(y_val, yhat_val)
			yhat_val[yhat_val >= best_thresh] = 1
			yhat_val[yhat_val < best_thresh] = 0
			val_f1 = f1_score(y_val, yhat_val, pos_label=1)

			yhat_test = model.predict(x_test)
			test_auroc = roc_auc_score(y_test, yhat_test)
			yhat_test[yhat_test >= best_thresh] = 1
			yhat_test[yhat_test < best_thresh] = 0
			test_f1 = f1_score(y_test, yhat_test, pos_label=1)

			if args.save_model.lower() in ['t', 'true']:
				model.save(args.save + '_model_' + str(n) + '.h5')

			final_results = final_results.append({'ID': args.save, 'Tag': args.tag,
				'Rep': n, 'X_file': args.x, 'Y_file': args.y,
				'ActFun': args.activation, 'dropout': args.dropout, 
				'L2': args.l2, 'LearnRate': args.learn_rate,
				'Optimizer': args.optimizer, 'n_Kernels': args.filters,
				'F1_threshold': best_thresh, 'n_Dense': args.dense,
				'Acc_train': train_acc, 'Loss_train': train_loss,
				'auROC_train': train_auroc, 'F1_train': train_f1,
				'Acc_val': val_acc, 'Loss_val': val_loss,
				'auROC_val': val_auroc, 'F1_val': val_f1,
				'Acc_test': test_acc, 'Loss_test': test_loss,
				'auROC_test': test_auroc, 'F1_test': test_f1}, ignore_index=True)

			##########################
			## Model Interpretation ##
			##########################

			if (args.imp_m.lower() in ['t', 'true'] or
				args.imp_k.lower() in ['t', 'true']):
				# Step 1: Read in x data meta data
				key = pd.read_csv(args.original_df, sep='\t', index_col=0,)
				key_index_list = key.columns.str.split('_', expand=True).values
				key.columns = pd.MultiIndex.from_tuples([(x[1], x[0]) for x in key_index_list])
				key = key.sort_index(axis=1)
				motifs = key.columns.levels[0].values
				omic_stack = list(key[list(key.columns.levels[0])[0]])
				omic_stack.append('PA')

				# Calculate Motif importance (zero-out-each-feature)
				if args.imp_m.lower() in ['t', 'true']:
					motif_imp = np.empty((0, 2))
					model_mot_imp = model
					for mx in range(0, x_test.shape[2] - 1):
						x_test_tmp = np.copy(x_test)
						x_test_tmp[:, ..., mx, :] = 0
						yhat_m_imp = model_mot_imp.predict(x_test_tmp)
						auroc_m_imp = roc_auc_score(y_test, yhat_m_imp)
						imp_m_auc = test_auroc - auroc_m_imp
						motif_imp = np.vstack((motif_imp, np.array([motifs[mx],
							imp_m_auc])))
					motif_imp = pd.DataFrame(motif_imp, columns=['motif',
						'auROC_test_decrease'])
					if n == 0:
						motif_imps = motif_imp
					else:
						motif_imps = pd.merge(motif_imps, motif_imp, on='motif')

				# Calculate Kernel Importance (zero-out-weights)
				if args.imp_k.lower() in ['t', 'true']:
					all_weights = model.get_weights()
					all_weights_2 = all_weights.copy()
					print('Performing Leave-One-Kernel-Out importance analysis...')
					for kx in range(0, args.filters):
						orig_weights = all_weights[0][:, :, 0, kx].copy()
						orig_weights = orig_weights.tolist()
						orig_weights = [i for l in orig_weights for i in l]
						conv2d_drop = copy.deepcopy(all_weights)
						conv2d_drop[0][:, :, 0, kx] = 0.0
						print(conv2d_drop[0][1, :, 0, 0:10])
						model_LOKO = tf.keras.models.clone_model(model)
						model_LOKO.set_weights(weights=conv2d_drop)
						yhat_k_imp = model_LOKO.predict(x_test)
						auroc_k_imp = roc_auc_score(y_test, yhat_k_imp)
						imp_k_auc = test_auroc - auroc_k_imp
						old = roc_auc_score(y_test, model.predict(x_test))
						print(old, imp_k_auc)
						kern_imp.append([n, imp_k_auc, orig_weights])

		if args.imp_m.lower() in ['t', 'true']:
			print('Snapshor ot motif importance scores...')
			motif_imps = motif_imps.set_index('motif')
			motif_imps = motif_imps.apply(pd.to_numeric, errors='coerce')
			motif_imps['mean_imp'] = motif_imps.mean(axis=1)
			motif_imps = motif_imps.sort_values('mean_imp', 0, ascending=False)
			print(motif_imps['mean_imp'].head())
			motif_imps['mean_imp'].to_csv(args.save + "_Motif_imp", sep="\t", index=True)

		if args.imp_k.lower() in ['t', 'true']:
			print('\nSnapshot of kernel importance scores:')
			kern_imp = pd.DataFrame(kern_imp,
				columns=['rep', 'auROC_test_decrease', 'kernel'])
			print(kern_imp.head())
			kern_imp.to_csv(args.save + "_Kernel_imp", sep="\t", index=True)

		final_results.to_csv(args.save + "_results.txt", header=True, sep='\t')

		# Save summary of results to RESULTS.txt
		calc_cols = ['F1_threshold', 'Acc_train', 'Acc_val', 'Acc_test',
		'Loss_train', 'Loss_val', 'Loss_test', 'auROC_train', 'auROC_val',
		'auROC_test', 'F1_train', 'F1_val', 'F1_test']
		final_results = final_results.drop(['Rep'], axis=1)
		std = final_results[calc_cols].std(axis=0, skipna=True)
		std = std.add_suffix('_std')
		mean = final_results[calc_cols].mean(axis=0, skipna=True)
		mean = mean.add_suffix('_mean')
		str_cols = final_results.drop(calc_cols, axis=1).iloc[0]
		str_cols = str_cols.append(pd.Series([args.n_reps], index=['Reps']))
		summary = pd.concat([str_cols, mean, std])

		#summary.set_index('index', inplace=True)
		print('\n### Summary of results on test set ###')
		print(summary.filter(like='test_mean', axis=0) )
		with open("RESULTS.txt", 'a') as f:
			summary.to_frame().transpose().to_csv(f, header=f.tell() == 0, sep='\t')

	print('Done!')


if __name__ == '__main__':
	main()
