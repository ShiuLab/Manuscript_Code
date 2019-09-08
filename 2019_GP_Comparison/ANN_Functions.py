from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys, os
import subprocess as sp
import numpy as np
import pandas as pd
import tensorflow as tf
from datetime import datetime
import timeit
import math

class fun(object):
	def __init__(self, filename):
		self.tokenList = open(filename, 'r')


	def train_valid_test_split(df, ho, y_name, val_perc):
	  """ Splitting data into training, validation (using val_perc) and testing (using holdout)"""

	  with open(ho) as ho_file:
	      ho_instances = ho_file.read().splitlines()
	      num_ho = len(ho_instances)
	  try:
	      test = df.loc[ho_instances, :]
	      train = df.drop(ho_instances)
	  except:
	      ho_instances = [int(x) for x in ho_instances]
	      test = df.loc[ho_instances, :]
	      train = df.drop(ho_instances)

	  val_set_index = np.random.rand(len(train)) < val_perc
	  valid = train[val_set_index]
	  train = train[~val_set_index]

	  X_train = train.drop(y_name, axis=1).values
	  X_valid = valid.drop(y_name, axis=1).values
	  X_test = test.drop(y_name, axis=1).values
	  Y_train = train.loc[:, y_name].values
	  Y_valid = valid.loc[:, y_name].values
	  Y_test = test.loc[:, y_name].values
	  X = df.drop(y_name, axis=1).values
	  Y = df.loc[:, y_name].values

	  return X, Y, X_train, X_valid, X_test, Y_train, Y_valid, Y_test

	def initialize_starting_weights(WEIGHTS, n_input, n_classes, archit, layer_number, df, mu, sigma, X_train, Y_train):
		"""
		Initializes starting weights and biases, with weights determined by:
			- random	Standard normal distribution 
			- xavier 	Useful when # features is large (Glorot and Yoshua (2010))
			- file 		Use weights from input file with noise infusion (options: mu, sigma)
		"""
		weights = {}
		biases = {}

		if WEIGHTS == 'random':
			weights['h1'] = tf.Variable(tf.random_normal([n_input, archit[0]]))
			biases['b1'] = tf.Variable(tf.random_normal([archit[0]]))
			for l in range(1,layer_number):
				w_name = 'h' + str(l+1)
				b_name = 'b' + str(l+1)
				weights[w_name] = tf.Variable(tf.random_normal([archit[l-1], archit[l]]))
				biases[b_name] = tf.Variable(tf.random_normal([archit[l]]))
			weights['out'] = tf.Variable(tf.random_normal([archit[-1], n_classes]))
			biases['out'] = tf.Variable(tf.random_normal([n_classes]))
		
		else:
			initializer = tf.contrib.layers.xavier_initializer()

			if WEIGHTS == 'xavier':
				weights['h1'] = tf.Variable(initializer([n_input, archit[0]]))

			else: # If weights is a file with features and weights:
				# Use Xavier weights for non-assigned nodes
				xav_init = tf.Variable(initializer([n_input, archit[0]]))
				init = tf.global_variables_initializer()
				with tf.Session() as sess:
					sess.run(init)
					xav = sess.run(xav_init)
				xav_mean = np.absolute(xav).mean().mean()
				xav_sd = xav.std()
				
				if WEIGHTS.lower() == 'rf':
					w_file = fun.Weights_RF(X_train, Y_train)
				
				elif WEIGHTS.lower() in ['rrblup', 'rrb']:
					w_file = fun.Weights_rrBLUP(X_train, Y_train)

				elif WEIGHTS.lower() in ['bl', 'bayesa', 'bayesb', 'brr']:
					w_file = fun.Weights_BGLR(WEIGHTS, X_train, Y_train)

				else:
					# Read in input weights
					w_file = pd.read_csv(WEIGHTS, sep=',', index_col = None, header=None, names=['ID', 'weight'])
					w_file = w_file[w_file['ID'].isin(list(df))] # Drop weights not included in x
					w_file.ID = w_file.ID.astype('category') # Convert to category 
					w_file.ID.cat.set_categories(list(df), inplace=True) # Set sorter for categories as X order
					w_file = w_file.sort_values(['ID'])
				w_mean = w_file['weight'].abs().mean().mean()

				# Scale input weights so that mean matches Xavier weights mean & add noise
				scale_by = xav_mean / w_mean
				print('Scaling input weights by: %.5f, with noise sd: %.5f' % (scale_by, xav_sd))
				weight_df = pd.concat([w_file['weight']]*archit[0], axis=1, ignore_index=True)
				noise = np.random.normal(mu, xav_sd, [n_input, archit[0]])
				weights_noise = (weight_df * scale_by) + noise
				weights_noise = weights_noise.astype(np.float32)
				
				# Set 75% of nodes to be from Xavier and 25% from input weights
				#per75 = int(archit[0]*0.75) # 25% seeded
				per75 = int(archit[0]*0.5) # 50% seeded
				#per75 = int(archit[0]*0.25) # 75% seeded
				#per75 = int(archit[0]*0.02) # 100% seeded
				weights_noise.iloc[:,0:per75] = xav[:,0:per75]
				weights['h1'] = tf.Variable(weights_noise)
			
			
			biases['b1'] = tf.Variable(tf.random_normal([archit[0]]))
			for l in range(1,layer_number):
				w_name = 'h' + str(l+1)
				b_name = 'b' + str(l+1)
				weights[w_name] = tf.Variable(initializer([archit[l-1], archit[l]]))
				biases[b_name] = tf.Variable(tf.random_normal([archit[l]]))
			weights['out'] = tf.Variable(tf.random_normal([archit[-1], n_classes]))
			biases['out'] = tf.Variable(tf.random_normal([n_classes]))


		return weights, biases

	def define_architecture(arc):
		# Define Architecture
		try:
			hidden_units = arc.strip().split(',')
			archit = list(map(int, hidden_units))
		except:
			archit = [arc]
		layer_number = len(archit)

		return archit, layer_number
	
	def define_loss(loss_type, nn_y, pred, l2, weights):
		"""Define loss function, taking into account L2 penalty if needed"""
		
		if loss_type == 'mse':
	 		loss = tf.losses.mean_squared_error(nn_y, pred)
	 		if l2 != 0.0:
	 			try:
	 				regularizer = tf.nn.l2_loss(weights['h1']) + tf.nn.l2_loss(weights['h2'])
	 			except:
	 				regularizer = tf.nn.l2_loss(weights['h1'])
	 			loss = tf.reduce_mean(loss + l2 * regularizer)
	 		else:
	 			loss = tf.reduce_mean(loss)
		else:
	 		print('That loss function is not implemented...')
	 		quit()

		return loss

	def multilayer_perceptron(x, weights, biases, layer_number, activation_function, dropout):
		"""
		Generate MLP model of any shape with options for activation function type and dropout levels
		currently: sigmoid, relu, elu
		"""
		layer = x
		for l in range(1,layer_number+1):
			weight_name = 'h' + str(l)
			bias_name = 'b' + str(l)
			layer = tf.add(tf.matmul(layer, weights[weight_name]), biases[bias_name])
			if activation_function.lower() == 'sigmoid':
				layer = tf.nn.sigmoid(layer)
			elif activation_function.lower() == 'relu':
				layer = tf.nn.relu(layer)
			elif activation_function.lower() == 'elu':
				layer = tf.nn.elu(layer)
			else:
				print("Given activation function is not supported")
				quit()
			if dropout != 0:
				dropout_rate = 1.0 - dropout
				drop_out = tf.nn.dropout(layer, dropout_rate)
		out_layer = tf.matmul(layer, weights['out']) + biases['out']
		return out_layer

	def save_trained_weights(SAVE, tvars, tvars_vals, archit, feat_list):
		for var, val in zip(tvars, tvars_vals):
			if var.name == 'Variable:0':
				col_names = ['Node_1_%s' % s for s in range(1,archit[0]+1)] 
				weights_l1 = pd.DataFrame(val, index=feat_list, columns=col_names)
				print('\n\n\nSnapshot of weights in first hidden layer:')
				print(weights_l1.head())

		if len(archit) == 1:
			for var, val in zip(tvars, tvars_vals):
				if var.name == 'Variable_2:0':
					weights_final = pd.DataFrame(val, index=list(weights_l1), columns=['Final'])

		elif len(archit) == 2:
			for var, val in zip(tvars, tvars_vals):
				if var.name == 'Variable_2:0':
					col_names = ['Node_2_%s' % s for s in range(1,archit[1]+1)] 
					weights_l2 = pd.DataFrame(val, index=list(weights_l1), columns=col_names )
					weights_l2.to_csv(SAVE+'_Weights_HL2.csv', index=True)
				elif var.name == 'Variable_4:0':
					weights_final = pd.DataFrame(val, index=list(weights_l2), columns=['Final'])

		elif len(archit) == 3:
			for var, val in zip(tvars, tvars_vals):
				if var.name == 'Variable_2:0':
					col_names = ['Node_2_%s' % s for s in range(1,archit[1]+1)] 
					weights_l2 = pd.DataFrame(val, index=list(weights_l1), columns=col_names )
					weights_l2.to_csv(SAVE+'_Weights_HL2.csv', index=True)
				elif var.name == 'Variable_4:0':
					col_names = ['Node_3_%s' % s for s in range(1,archit[2]+1)] 
					weights_l3 = pd.DataFrame(val, index=list(weights_l2), columns=col_names )
					weights_l3.to_csv(SAVE+'_Weights_HL3.csv', index=True)
				elif var.name == 'Variable_6:0':
					weights_final = pd.DataFrame(val, index=list(weights_l3), columns=['Final'])

		weights_l1.to_csv(SAVE+'_Weights_HL1.csv', index=True)
		weights_final.to_csv(SAVE+'_Weights_Fin.csv', index=True)
	
	def conv2d(data_in, W, conv_r, conv_c):
		return tf.nn.conv2d(data_in, W, strides=[1,1,1,1], padding='SAME') # strides=[1, 1, 1, 1]

	def maxpool2d(x):
  		return tf.nn.max_pool(x, ksize=[1, 2, 2, 1], strides=[1, 2, 2, 1], padding='SAME')

	def convolutional_neural_network(x, conv_r, conv_c, shape_r, shape_c, dropout, actfun):
		shape_r_2pool = int(math.ceil(shape_r/4))
		shape_c_2pool = int(math.ceil(shape_c/4))
		weights = {
		'W_conv1': tf.Variable(tf.random_normal([conv_r, conv_c, 1, 32])),
		  'W_conv2': tf.Variable(tf.random_normal([conv_r, conv_c, 32, 64])),
		  'W_fc': tf.Variable(tf.random_normal([shape_r_2pool*shape_c_2pool*64, 1024])),
		  'out': tf.Variable(tf.random_normal([1024, 1]))}
		biases = {
		'b_conv1': tf.Variable(tf.random_normal([32])),
		  'b_conv2': tf.Variable(tf.random_normal([64])),
		  'b_fc': tf.Variable(tf.random_normal([1024])),
		  'out': tf.Variable(tf.random_normal([1]))}

		# Reshape input to a 4D tensor 
		x = tf.reshape(x, shape=[-1, shape_r, shape_c, 1])

		# Convolution Layer 1 - Max Pooling 1 - Convlution Layer 2 - Max Pooling 2
		if actfun.lower() == 'relu':
			conv1 = tf.nn.relu(fun.conv2d(x, weights['W_conv1'], conv_r, conv_c) + biases['b_conv1'])
			conv1 = fun.maxpool2d(conv1)
			conv2 = tf.nn.relu(fun.conv2d(conv1, weights['W_conv2'], conv_r, conv_c) + biases['b_conv2'])
			conv2 = fun.maxpool2d(conv2)
		elif actfun.lower() == 'sigmoid':
			conv1 = tf.nn.sigmoid(fun.conv2d(x, weights['W_conv1'], conv_r, conv_c) + biases['b_conv1'])
			conv1 = fun.maxpool2d(conv1)
			conv2 = tf.nn.sigmoid(fun.conv2d(conv1, weights['W_conv2'], conv_r, conv_c) + biases['b_conv2'])
			conv2 = fun.maxpool2d(conv2)
		else:
			print('Activation function not accepted. Use: sigmoid or relu')
			quit()

		# Fully connected layer
		fc = tf.reshape(conv2, [-1, shape_r_2pool * shape_c_2pool * 64]) # Reshape conv2 output to fit fully connected layer
		if actfun.lower() == 'relu':
			fc = tf.nn.relu(tf.matmul(fc, weights['W_fc']) + biases['b_fc'])
		if actfun.lower() == 'sigmoid':
			fc = tf.nn.sigmoid(tf.matmul(fc, weights['W_fc']) + biases['b_fc'])
		if dropout > 0:
			print('Applying Dropout (dropout) with dropout rate = %f' % dropout)
			fc = tf.nn.dropout(fc, 1-dropout) # Keep rate = 1 - dropout rate
		
		output = tf.matmul(fc, weights['out']) + biases['out']

		return output 

	def Image2Features(X_file, shape_r, shape_c):
		"""
		Read in images from directory (X_file), scale to fit shape, then flatten into one X file.
		"""
		from PIL import Image

		if os.path.isfile(X_file + 'X_processed.csv'):
			print('Pulling already processed image data from %sX_processed.csv' % X_file)
			x = pd.read_csv(X_file + 'X_processed.csv' , sep=',', index_col = 0)

		else:
			print('Converting images to pixle feature file...')
			col_names = ["p_" + str(i+1) for i in range(shape_r*shape_c)]
			x = pd.DataFrame(columns=col_names)
			size = shape_r, shape_c

			for f in os.listdir(X_file):
				if f.startswith("."):
					pass
				else:
					try:
						print(f)
						img = Image.open(X_file + f).convert('L') # Converts to 8-bit grayscale
						img.thumbnail(size)
						img_data = list(img.getdata())
						x.loc[f] = img_data
					except:
						print('%s would not convert to pixels...' % f)
			x.to_csv(X_file + 'X_processed.csv', sep=',')

		print(x.head())
		print(x.shape)
		return x 

	def Weights_RF(X_train, Y_train):
		""" Run Random Forest to get coefficients to use as starting weights """
		from sklearn.ensemble import RandomForestRegressor
		from math import sqrt
		print('Calculating seeded starting weights using RF')
		feat_sel_forest = RandomForestRegressor(max_features= round(sqrt(X_train.shape[1])), n_estimators=500, n_jobs=1)
		feat_sel_forest = feat_sel_forest.fit(X_train, Y_train)
		w_file = pd.DataFrame(data=feat_sel_forest.feature_importances_, columns=['weight'])
		return(w_file)

	def Weights_rrBLUP(X_train, Y_train):
		""" Run rrBLUP to get coefficients to use as starting weights """
		from random import randint
		cwd = os.getcwd()
		tmp = str(randint(1000000, 9999999))
		np.savetxt('x_file_'+ tmp, X_train, delimiter=',')
		np.savetxt('y_file_'+ tmp, Y_train, delimiter=',')

		# Write temp Rscript  
		tmpR=open("rrB_%s.R" % tmp,"w")
		tmpR.write("setwd('%s')\n" % cwd)
		tmpR.write('library(rrBLUP)\n')
		tmpR.write('library(data.table)\n')
		tmpR.write("X <- fread('%s', sep=',', header=F)\n" % ('x_file_'+ tmp))
		tmpR.write("Y <- fread('%s', sep=',', header=F)\n" % ('y_file_'+ tmp))
		tmpR.write("mod <- mixed.solve(Y$V1, Z=X, K=NULL, SE=FALSE, return.Hinv=FALSE)\n")
		tmpR.write("weight <- mod$u\n")
		tmpR.write("weight_df <- as.data.frame(weight)\n")
		tmpR.write("write.table(weight_df, file='%s', sep=',', row.names=TRUE, quote=FALSE)\n" % ('rrBScores_'+tmp+'.txt'))
		tmpR.close()

		process = sp.Popen('module load R && export R_LIBS_USER=~/R/library && Rscript rrB_%s.R' % tmp, shell=True)
		process.wait()
		w_file = pd.read_csv('rrBScores_'+tmp+'.txt', sep = ',')
		os.system("rm rrBScores_%s.txt rrB_%s.R x_file_%s y_file_%s" % (tmp, tmp, tmp, tmp))
		
		return(w_file)


	def Weights_BGLR(WEIGHTS, X_train, Y_train):
		""" Run BL, BayesA, BayesB, or BRR using the BGLR package implemented in R to get seeded starting weights """
		from random import randint
		cwd = os.getcwd()
		tmp = str(randint(1000000, 9999999))
		np.savetxt('x_file_'+ tmp, X_train, delimiter=',')
		np.savetxt('y_file_'+ tmp, Y_train, delimiter=',')

		# Write temp Rscript  
		tmpR=open("BGLR_%s.R" % tmp,"w")
		tmpR.write("setwd('%s')\n" % cwd)
		tmpR.write('library(BGLR)\n')
		tmpR.write('library(data.table)\n')
		tmpR.write("X <- fread('%s', sep=',', header=F)\n" % ('x_file_'+ tmp))
		tmpR.write("Y <- fread('%s', sep=',', header=F)\n" % ('y_file_'+ tmp))
		tmpR.write("ETA=list(list(X=X,model='%s'))\n" % WEIGHTS)
		tmpR.write("fm=BGLR(y=Y$V1,ETA=ETA,verbose=FALSE, nIter=12000,burnIn=2000)\n")
		tmpR.write("weight <- fm$ETA[[1]]$b\n")
		tmpR.write("weight_df <- as.data.frame(weight)\n")
		tmpR.write("write.table(weight_df, file='%s', sep=',', row.names=TRUE, quote=FALSE)\n" % ('BGLRScores_'+tmp+'.txt'))
		tmpR.close()

		process = sp.Popen('module load R && export R_LIBS_USER=~/R/library && Rscript BGLR_%s.R' % tmp, shell=True)
		process.wait()
		w_file = pd.read_csv('BGLRScores_'+tmp+'.txt', sep = ',')
		os.system("rm BGLRScores_%s.txt BGLR_%s.R x_file_%s y_file_%s" % (tmp, tmp, tmp, tmp))
		
		return(w_file)



