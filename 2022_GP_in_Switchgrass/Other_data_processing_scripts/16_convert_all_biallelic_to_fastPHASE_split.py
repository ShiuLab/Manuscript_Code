import sys,os
inp = open(sys.argv[1],'r').readlines()
number = int(sys.argv[2])
startnum = int(sys.argv[3])
for n in range(startnum,number//5000 + 1):
	out = open(sys.argv[1].split('.txt')[0] + '_fastPHASE_%s_%s.txt'%((n-1)*5000+1, n*5000),'w')
	out.write('%s\n'%(len(inp[0].split('\t')) - 4))
	out.write('%s\n'%5000)
	tem = inp[0].strip().split('\t')
	for i in range(4,len(inp[0].split('\t'))):
		ind = tem[i]
		out.write('# id %s\n'%ind)
		out1 = ''
		out2 = ''
		for inl in inp[((n-1)*5000+1):(n*5000+1)]:
			ref = inl.split('\t')[2]
			alt = inl.split('\t')[3]
			if len(alt) > len(ref):
				indel = 'I'
			elif len(alt) < len(ref):
				indel = 'D'
			snp = inl.strip().split('\t')[i]
			snp = snp.replace('.','?')
			if len(alt) > len(ref):
				snp = snp.replace(alt,'I')
				snp = snp.replace(ref,'R')
			if len(alt) < len(ref):
				snp = snp.replace(ref,'R')
				snp = snp.replace(alt,'D')
			allele = snp.split('/')
			if out1 == '':
				out1 = allele[0]
				out2 = allele[1]
			else:
				out1 = out1 + ' ' + allele[0]
				out2 = out2 + ' ' + allele[1]
		out.write(out1 + '\n')
		out.write(out2 + '\n')

	out.close()

if 1==1:	
	out = open(sys.argv[1].split('.txt')[0] + '_fastPHASE_%s_%s.txt'%((number//5000) * 5000 + 1, number),'w')
	out.write('%s\n'%(len(inp[0].split('\t')) - 4))
	out.write('%s\n'%(number - (number//5000) * 5000))
	tem = inp[0].strip().split('\t')
	for i in range(4,len(inp[0].split('\t'))):
		ind = tem[i]
		out.write('# id %s\n'%ind)
		out1 = ''
		out2 = ''
		for inl in inp[((number//5000) * 5000 + 1):(number + 1)]:
			ref = inl.split('\t')[2]
			alt = inl.split('\t')[3]
			if len(alt) > len(ref):
				indel = 'I'
			elif len(alt) < len(ref):
				indel = 'D'
			snp = inl.strip().split('\t')[i]
			snp = snp.replace('.','?')
			if len(alt) > len(ref):
				snp = snp.replace(alt,'I')
				snp = snp.replace(ref,'R')
			if len(alt) < len(ref):
				snp = snp.replace(ref,'R')
				snp = snp.replace(alt,'D')
			allele = snp.split('/')
			if out1 == '':
				out1 = allele[0]
				out2 = allele[1]
			else:
				out1 = out1 + ' ' + allele[0]
				out2 = out2 + ' ' + allele[1]
		out.write(out1 + '\n')
		out.write(out2 + '\n')

	out.close()
				
		