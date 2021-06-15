import sys,os
import pandas as pd
'''
input1: classification of SNPs
input2: SNP matrix
'''
type = open(sys.argv[1],'r').readlines()
inp = open(sys.argv[2],'r').readlines()
D = {}	
for inl in type[1:]:
	tem = inl.split('\t')
	chr = tem[0]
	pos = tem[1]
	type = tem[4]
	allelic = int(tem[5])
	if allelic > 2:
		D[chr + '-' + pos] = 1

out = open(sys.argv[2].split('.txt')[0] + '_multi-allelic.txt','w')
out.write(inp[0])
for inl in inp[1:]:
	tem = inl.split('\t')
	chr = tem[0]
	pos = tem[1]
	if chr + '-' + pos in D:
		out.write(inl)

out.close()
		
