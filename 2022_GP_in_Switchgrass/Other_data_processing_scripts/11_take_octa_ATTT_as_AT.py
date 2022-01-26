import sys,os
import numpy
inp = open(sys.argv[1],'r').readlines()
out = open(sys.argv[1] + '_converted_ATTT_to_AT.txt','w')
out.write(inp[0])
for inl in inp[1:]:
	tem = inl.strip().split('\t')
	ref = tem[2]
	alt = tem[3]
	for i in range(4,len(tem)):
		snp = tem[i].split('/')
		# if len(numpy.unique(snp)) == 1 and snp[0]==ref:
			# tem[i] = '%s/%s'%(ref,ref)
		# if len(numpy.unique(snp)) == 1 and snp[0]==alt:
			# tem[i] = '%s/%s'%(alt,alt)
		if len(numpy.unique(snp)) == 1:
			tem[i] = '%s/%s'%(numpy.unique(snp)[0],numpy.unique(snp)[0])
		# if len(numpy.unique(snp)) == 2:
			# tem[i] = '%s/%s'%(ref,alt)
		if len(numpy.unique(snp)) == 2:
			tem[i] = '%s/%s'%(numpy.unique(snp)[0],numpy.unique(snp)[1])
		if len(numpy.unique(snp)) > 2:
			print(snp)
		if tem[i] == './././.':
			tem[i] = './.'
	out.write('\t'.join(tem) + '\n')
	
out.close()