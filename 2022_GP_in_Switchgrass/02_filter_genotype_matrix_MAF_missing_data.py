import sys,os,argparse
import pandas as pd
import numpy as np

def warn(*args, **kwargs):
	pass
import warnings
warnings.warn = warn

def main():
	parser = argparse.ArgumentParser(description='This code is for filtering the genotype matrix with the criteria that: missing data < 20%, MAF(the second most common allele) > 0.05')
	# Required
	req_group = parser.add_argument_group(title='REQUIRED INPUT')
	req_group.add_argument('-file', help='genotype matrix file', required=True)
	
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(0)
	args = parser.parse_args()	

	file = args.file
	inp = open(file,'r')
	out = open(file + '_filtered','w')
	inl = inp.readline()
	while inl:
		tem = inl.strip().split('\t')
		if inl.startswith('Chr\t'):
			title = 'Chr\tPos\tRef\tAlt'
			for n in tem[4:]:
				title = title + '\t' + n
			out.write(title + '\n')
			out.flush()
		else:
			keep = 0
			var = {}
			for n in tem[4:]:
				if n != './././.' and n != './.':
					keep += 1
					### here, take tetraploid as the same as octoploid
					if len(n.split('/')) == 2:
						n = '%s/%s/%s/%s'%(n.split('/')[0],n.split('/')[0],n.split('/')[1],n.split('/')[1])
					### count the frequency of each allele
					for allele in n.split('/'):
						if allele not in var:
							var[allele] = 1
						else:
							var[allele] += 1
			if len(var) > 1:
				L = []
				for allele in var:
					L.append(var[allele])
				L.sort()
				### if missing data < 20%, MAF(the second most common allele) > 0.05
				if float(keep)/(len(tem)-4) >= 0.8 and float(L[-2])/(keep*4) > 0.05:
					out.write(inl)
					out.flush()
		inl = inp.readline()

	out.close()

if __name__ == '__main__':
	main()
