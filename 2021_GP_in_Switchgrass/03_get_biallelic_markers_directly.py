'''
classify the variation into SNP, indel, or SNP/indel; biallelic or non-biallelic; in genic or intergenic, three_UTR or five_UTR region, exonic or intronic, splicing regions
input1: genotype matrix, which has been filtered with MAF and missing data
'''
import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def judge(ref,alt):
	T = {}
	for snp in alt.split(','):
		if len(ref) == len(snp) and '*' not in snp:
			T['SNP'] = 1
		else:
			T['indel'] = 1
	return '/'.join(sorted(T.keys()))

def main():
	parser = argparse.ArgumentParser(description='This code is for getting the biallelic SNPs or indels directly from the genotype matrix')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-file', help='genotype matrix, which has been filtered with MAF and missing data', required=True)
	req_group.add_argument('-type', help='indel or SNP', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	file = args.file
	inp = open(file,'r').readlines()
	out = open(file + '_biallelic_%s.txt'%args.type,'w')
	out.write(inp[0])
	for inl in inp[1:]:
		tem = inl.split('\t')
		chr = tem[0]
		pos = int(tem[1])
		ref = tem[2]
		alt = tem[3].strip()
		type = judge(ref,alt)
		allelic = len(alt.split(',')) + 1
		if allelic == 2:
			if type == args.type:
				out.write(inl)

	out.close()


if __name__ == '__main__':
	main()










