import sys,os
inp = open(sys.argv[1],'r')
out = open(sys.argv[1].split('.txt')[0] + '_One_hot_encoding.txt','w')
inl = inp.readline()
out.write(inl)
inl = inp.readline()
while inl:
	tem = inl.strip().split('\t')
	pos = tem[1]
	ref = tem[2]
	alt = tem[3].split(',')
	for alternative in alt:
		tem_new = inl.strip().split('\t')
		tem_new[1] = pos+'_'+alternative
		tem_new[3] = alternative
		for i in range(4,len(tem_new)):
			snp = tem_new[i]
			alle1 = snp.split('/')[0]
			alle2 = snp.split('/')[1]
			if (alle1 == alternative and alle2 != alternative) or (alle1 != alternative and alle2 == alternative):
				tem_new[i] = '%s/%s'%(ref,alternative)
			if alle1 != alternative and alle2 != alternative and alle1 != '.':
				tem_new[i] = '%s/%s'%(ref,ref)
		res = '\t'.join(tem_new)
		out.write(res + '\n')
	inl = inp.readline()
	
out.close()

