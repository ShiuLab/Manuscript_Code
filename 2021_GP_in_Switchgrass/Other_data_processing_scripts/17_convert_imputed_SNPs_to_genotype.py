import sys,os
'''
input1: biallelic-SNPs for only tetraploids
input2: imputed SNPs
'''
inp = open(sys.argv[1],'r').readlines()
impute = open(sys.argv[2],'r').readlines()
out = open(sys.argv[1] + '_imputed_geno.csv','w')
out2 = open(sys.argv[1] + '_imputed.txt','w')
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
	
	
	
