"""
Generate feature selection subsets for a given Random Forest model
"""

import sys
import argparse
import joblib
import datatable as dt
import pandas as pd
import numpy as np


def parse_arguments():
	"""Argument parser"""
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-m", "--model", type=str, help="path to Random Forest (RF) model file",
		default=""
	)
	parser.add_argument(
		"-i", "--imp", type=str, help="path to RF feature importance file",
		default=""
	)
	parser.add_argument(
		"-start", type=int, help="Number of features to start with",
		required=True
	)
	parser.add_argument(
		"-stop", type=int, help="Number of features to end with",
		required=True
	)
	parser.add_argument(
		"-o", "--save", type=str, help="path and prefix to save feature subset as",
		required=True
	)
	parser.add_argument(
		"-step", type=int, help="Step size of features to select", default=None
	)
	parser.add_argument(
		"-base", type=int, help="Base of exponential number of features to select, e.g., base^start to base^stop",
		default=None
	)
	parser.add_argument(
		"-x", "--X_file", type=str, help="Path to input feature table of model (for sklearn models using < v1.0",
		required=False
	)
	parser.add_argument(
		"-d", "--drop", type=str, help="Whether to drop features with zero gini importance before feature selection (y/n)",
		required=False, default="n")
	return parser.parse_args()

if __name__ == "__main__":
	# Arguments help
	if len(sys.argv)==1:
		print("python 0_feature_selection.py [-m optional] [-i optional] [-start] [-stop] [-step] [-o] [-x optional]")
		sys.exit()
	
	# Arguments
	args = parse_arguments()
	model_file = args.model
	imp_file = args.imp
	start = args.start
	stop = args.stop
	step = args.step
	base = args.base
	save = args.save
	drop = args.drop

	# Read in feature importances
	if model_file != "": # From the saved model
		mod = joblib.load(open(model_file, "rb"))
		try:
			# .feature_names_in_ only works for models built using scikit-learn version > 1.0
			imp = mod.feature_importances_
			imp = pd.DataFrame(imp, index=mod.feature_names_in_)
		except:
			X_file = args.X_file
			X = dt.fread(X_file).to_pandas()
			imp = mod.feature_importances_
			imp = pd.DataFrame(imp, index=X.columns[1:])

	if imp_file != "": # From the feature importance file
		imp = pd.read_csv(imp_file, sep="\t", index_col=0)
		imp = pd.DataFrame(imp["mean_imp"])
	
	# Sort and remove unimportant features
	try:
		imp.sort_values(by=0, ascending=False, inplace=True) # sort features
	except:
		imp.sort_values(by="mean_imp", ascending=False, inplace=True) # sort features
	
	if drop=="y":
		og_featnum = len(imp)
		imp = imp.loc[imp.iloc[:,0]!=0.0,:] # drop non-important features
		print(f"After removing unimportant features, only {len(imp)} out of {og_featnum} remained.")

	# Select subsets of features
	if stop > len(imp):
		print(f"Setting stop argument to {len(imp)}")
		stop = len(imp)

	if step != None:
		for i in range(start, stop+step, step):
			if i == start:
				continue
			else:
				subset = imp.iloc[start:i,:].index.values
				print(f"saving {start} to {i} features to file {save}_top_{len(subset)}")
				np.savetxt(f"{save}_top_{len(subset)}", subset, fmt="%s")
	
	if base != None:
		for i in range(start, stop+1):
			if i == start-1:
				continue
			else:
				subset = imp.iloc[start-1:base**i,:].index.values
				print(f"saving {start-1} to {base**i} features to file {save}_top_{len(subset)}")
				np.savetxt(f"{save}_top_{len(subset)}", subset, fmt="%s")