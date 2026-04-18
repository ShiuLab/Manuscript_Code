###########################################################
# Sum the SHAP feature interaction values across instances
# 
#
# Written by: Peipei Wang
# Modified by: Kenia Segura Abá
###########################################################

import sys, os, argparse, warnings, gc
import datatable as dt
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
from scipy import sparse

def warn(*args, **kwargs):
	pass
warnings.warn = warn


def calc_total_shap(path, file, idx, n):
	'''Calculate the total SHAP interaction values for each feature pair.
	The on-diagonal effects correspond to the difference between the SHAP value
	and the sum of the off-diagonal SHAP interaction values for a feature.'''
	shap = sparse.load_npz(path + file)
	shap = pd.DataFrame(shap.toarray(), index=idx, columns=idx)
	phi_ij = np.triu(shap, k=1) # upper triangle
	phi_ji = np.triu(shap.T, k=1) # lower triangle
	tot_phi = phi_ij + phi_ji # total interaction effects
	tot_phi = pd.DataFrame(tot_phi, index=idx, columns=idx)
	rem = pd.Series(np.diag(shap), index=idx) # remaining effects
	phi = tot_phi.where(np.triu(tot_phi, k=1).astype(bool)).stack()
	rem_name = file.split('.')[0].split('_')[-1] # instance name
	phi_name = file.split('.')[0].split('_')[-1] # instance name
	rem = {rem_name: rem}
	phi = {phi_name: phi}
	return rem, phi


def combine_across_instances(rem, phi, SHAP_rem, SHAP_phi):
	'''Combine local interaction effects and remaining effects across instances.'''
	SHAP_rem.update(rem)
	SHAP_phi.update(phi)
	return SHAP_rem, SHAP_phi


def sum_across_instances(SHAP_rem, SHAP_phi):
	'''Sum the local interaction effects and remaining effects across instances.'''
	SHAP_rem = pd.DataFrame(SHAP_rem)
	SHAP_phi = pd.DataFrame(SHAP_phi)
	SHAP_rem["sum_across_instances"] = SHAP_rem.sum(axis=1)
	SHAP_phi["Interaction"] = SHAP_phi.sum(axis=1)
	return SHAP_rem, SHAP_phi


def main():
	parser = argparse.ArgumentParser(
		description='Sum the SHAP feature interaction values across instances.')
	
	### Input arguments ###
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-path', help='path saving the interaction score files', required=True)
	req_group.add_argument('-save', help='name to save', required=True)
	# Optional
	opt_group = parser.add_argument_group(title='OPTIONAL INPUT')
	opt_group.add_argument('-y', help='name of label', default='')
	opt_group.add_argument('-dtype', help='name of feature type', default='')
	
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()
	
	if not args.path.endswith('/'):
		path = args.path + '/'
	else:
		path = args.path
	
	n = 0
	SHAP_rem = {}
	SHAP_phi = {}
	if (args.y == '' and args.dtype == ''):
		# Get row and column feature names (.npz files don't have such info)
		for file in glob(path + 'shap_interaction_scores_*.txt'):
			if not file.endswith('_summed.txt'):
				idx = dt.fread(file).to_pandas()
				idx = idx.iloc[:,0]
				idx.name = 'ID'
		for file in tqdm(os.listdir(path)):
			if file.endswith('.npz') and file.startswith('shap_interaction_scores'):
				# Calculate local total shap interaction effects and remaining effects
				rem, phi = calc_total_shap(path, file, idx, n)
				SHAP_rem, SHAP_phi = combine_across_instances(rem, phi, SHAP_rem, SHAP_phi)
				n += 1
				del rem, phi
				gc.collect()
	else: 
		# Get row and column feature names (.npz files don't have such info)
		for file in glob(path + 'shap_interaction_scores_*.txt'):
			if (args.y in file and args.dtype in file and not file.endswith('_summed.txt')):
				idx = dt.fread(file).to_pandas()
				idx = idx.iloc[:,0]
				idx.name = 'ID'
		for file in tqdm(os.listdir(path)):
			if (file.endswith('.npz') and file.startswith('shap_interaction_scores')) \
				and (args.y in file and args.dtype in file):
				# Calculate local total shap interaction effects and remaining effects
				rem, phi = calc_total_shap(path, file, idx, n)
				SHAP_rem, SHAP_phi = combine_across_instances(rem, phi, SHAP_rem, SHAP_phi)
				n += 1
				del rem, phi
				gc.collect()
		
	# Sum effects across instances
	SHAP_rem, SHAP_phi = sum_across_instances(SHAP_rem, SHAP_phi)
	
	# Save results
	dt.Frame(SHAP_rem.reset_index()).to_csv(args.save + '_remaining_summed.txt', sep='\t', header=True)
	SHAP_phi.rename_axis(index=['Feature1', 'Feature2'], inplace=True)
	dt.Frame(SHAP_phi.reset_index()).to_csv(args.save + '_interaction_summed.txt', sep='\t', header=True)


if __name__ == '__main__':
	main()
