import sys,os
inp = open(sys.argv[1],'r')
out = open(sys.argv[1].split('_')[2].split('.')[0] + '_methylated.txt','w')
out.write('Loc\t%s\n'%sys.argv[1].split('_')[2].split('.')[0])
inl = inp.readline()
while inl:
	if not inl.startswith('chrom'):
		tem = inl.strip().split('\t')
		if tem[6]=='1':
			out.write('Chr%s_%s_%s_%s\t%s\n'%(tem[0],tem[1],tem[3],tem[2],tem[6]))
			out.flush()
	inl = inp.readline()

out.close()
