'''
input1: biallelic variation matrix
input2: imputed biallelic variation matrix
'''
import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for converting the imputed genotype matrix to the format previously used')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-matrix', help='biallelic variation matrix', required=True)
	req_group.add_argument('-imputed_matrix', help='imputed biallelic variation matrix', required=True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	
	
	inp = open(args.matrix,'r').readlines()
	impute = open(args.imputed_matrix,'r').readlines()
	out = open(args.matrix + '_imputed_geno.csv','w')
	out2 = open(args.matrix + '_imputed.txt','w')
	D = {}
	i = 1
	while i < len(inp):
		inl = inp[i]
		tem = inl.split('\t')
		D[i] = [tem[0]+ '_' +tem[1],tem[2],tem[3]] ### D[i] = [chr_pos,ref,alt]
		i += 1
		
	title = 'ID'
	i = 1
	while i <= len(D):
		title = title + ',' + D[i][0]
		i += 1

	out.write(title + '\n')

	out2.write(inp[0])
	R= {}

	x = 0
	while x < len(impute):
		inl = impute[x]
		if inl.startswith('# id'):
			ind = inl.split('# id ')[1].strip()
			R[ind] = {}
			res = ind
			x += 1
			a1 = impute[x]
			x += 1
			a2 = impute[x]
			allele1 = a1.split()
			allele2 = a2.split()
			y = 0
			while y < len(allele1):
				pos = D[y+1][0]
				ref = D[y+1][1]
				alt = D[y+1][2]
				if len(ref) > len(alt):
					ref = 'R'
					alt = 'D'
				if len(ref) < len(alt):
					ref = 'R'
					alt = 'I'
				if allele1[y] == allele2[y] and allele1[y]==ref:
					res = res + ',1'
				if allele1[y] == allele2[y] and allele1[y]==alt:
					res = res + ',-1'
				if allele1[y] != allele2[y]:
					res = res + ',0'
				R[ind][D[y+1][0]] = '%s/%s'%(allele1[y],allele2[y])  ### R[ind][chr_pos] = 'A/T'
				y += 1
			out.write(res + '\n')
		x += 1

	out.close()

	IND = inp[0].strip().split('\t')[4:]
	for inl in inp[1:]:
		tem = inl.split('\t')
		chr = tem[0]
		pos = tem[1]
		ref = tem[2]
		alt = tem[3]
		res = '%s\t%s\t%s\t%s'%(chr,pos,ref,alt)
		for ind in IND:
			res = res + '\t' + R[ind]['%s_%s'%(chr,pos)]
		out2.write(res + '\n')
		
	out2.close()
	
if __name__ == '__main__':
	main()
	
	
