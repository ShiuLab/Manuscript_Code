'''
input1: classification file of variation
input2: file with ploidy of each individual
input3: genotype matrix
input4: indel or SNP
'''

import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for extracting the biallelic SNPs or indels')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-classification', help='classification file of variation', required=True)
	req_group.add_argument('-file', help='genotype matrix', required=True)
	req_group.add_argument('-type', help='indel or SNP', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	type = open(args.classification,'r').readlines()
	inp = open(args.file,'r').readlines()

	D = {}	
	for inl in type[1:]:
		tem = inl.split('\t')
		chr = tem[0]
		pos = tem[1]
		type = tem[4]
		allelic = int(tem[5])
		if type == sys.argv[4] and allelic == 2:
			D[chr + '-' + pos] = 1

	out = open(args.file + '_biallelic_%s.txt'%args.type,'w')
	out.write(inp[0])
	for inl in inp[1:]:
		tem = inl.split('\t')
		chr = tem[0]
		pos = tem[1]
		if chr + '-' + pos in D:
			out.write(inl)

	out.close()

if __name__ == '__main__':
	main()

