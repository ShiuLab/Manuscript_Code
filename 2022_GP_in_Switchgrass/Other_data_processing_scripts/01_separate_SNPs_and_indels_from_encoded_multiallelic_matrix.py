import sys,os
import copy
inp = open(sys.argv[1],'r').readlines()
snp = open(sys.argv[1] + '_SNP','w')
indel = open(sys.argv[1] + '_indel','w')
snp.write(inp[0])
indel.write(inp[0])
D = {}
for k in range(1,len(inp)):
	inl = inp[k]
	tem = inl.strip().split('\t')
	if '%s__%s'%(tem[0],tem[1]) not in D:
		D['%s__%s'%(tem[0],tem[1])] = 1
		if len(tem[3])==1 and len(tem[2])==1 and tem[3] != '*':
			snp.write(inl)
		elif len(tem[3]) == len(tem[2]) and len(tem[3]) != 1:
			for i in range(0,len(tem[3])):
				if tem[3][i] != tem[2][i]:
					new = copy.deepcopy(tem)
					new[1] = '%s_%s'%(int(tem[1].split('_')[0]) + i - 1,tem[3][i]) 
					new[2] = tem[2][i]
					new[3] = tem[3][i]
					for j in range(4,len(tem)):
						new[j] = '%s/%s'%(tem[j].split('/')[0][i],tem[j].split('/')[1][i])
					snp.write('\t'.join(new) + '\n')
		else:
			indel.write(inl)

snp.close()
indel.close()