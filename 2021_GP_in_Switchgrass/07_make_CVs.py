import sys, os
import pandas as pd
import numpy as np
import math
import random
import sys,os,argparse


def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for making the CVs files')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-file', help='the pheno matrix', required=True)
	req_group.add_argument('-cv', help='the fold number for cross-validation', required=True)
	req_group.add_argument('-number', help='how many times you want to repeat your cross-validation scheme', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	df = pd.read_csv(args.file, sep=',', header =0, index_col = 0)
	cv = int(args.cv)
	number = int(args.number)

	cvs = pd.DataFrame(index = df.index)

	n_lines = len(df)
	n_reps = int((n_lines/cv) + 1) #math.ceil
	print(n_lines)

	for i in range(1,number+1):
		name = 'cv_' + str(i)
		mix = np.repeat(range(1,6), n_reps)
		np.random.shuffle(mix)
		cvs[name] = mix[0:n_lines]

	cvs.to_csv('CVFs.csv', sep=',')

if __name__ == '__main__':
	main()
