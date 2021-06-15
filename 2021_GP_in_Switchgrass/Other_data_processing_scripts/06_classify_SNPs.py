import sys,os
gff = open('/mnt/home/peipeiw/Documents/Genome_selection/Genotype_matrix/Pvirgatum_516_v5.1.gene.gff3','r').readlines()
G = {}
for inl in gff:
	if not inl.startswith('#'):
		tem = inl.split('\t')
		chr = tem[0]
		type = tem[2]
		name = tem[8].split('ID=')[1].split(';')[0].replace('.v5.1','')
		if type == 'gene':
			G[name] = {}
		if type == 'mRNA':
			gene = '.'.join(name.split('.')[0:2])
			G[gene][name] = {}
		if type == 'CDS':
			gene = '.'.join(name.split('.')[0:2])
			mRNA = '.'.join(name.split('.')[0:3])
			if 'CDS' not in G[gene][mRNA]:
				G[gene][mRNA]['CDS'] = []
			G[gene][mRNA]['CDS'].append(name)
		if type == 'five_prime_UTR':
			gene = '.'.join(name.split('.')[0:2])
			mRNA = '.'.join(name.split('.')[0:3])
			if 'five_prime_UTR' not in G[gene][mRNA]:
				G[gene][mRNA]['five_prime_UTR'] = []
			G[gene][mRNA]['five_prime_UTR'].append(name)
		if type == 'three_prime_UTR':
			gene = '.'.join(name.split('.')[0:2])
			mRNA = '.'.join(name.split('.')[0:3])
			if 'three_prime_UTR' not in G[gene][mRNA]:
				G[gene][mRNA]['three_prime_UTR'] = []
			G[gene][mRNA]['three_prime_UTR'].append(name)
			
GFF = {}
for inl in gff:
	if not inl.startswith('#'):
		tem = inl.split('\t')
		chr = tem[0]
		type = tem[2]
		left = int(tem[3])
		right = int(tem[4])
		dir = tem[6]
		name = tem[8].split('ID=')[1].split(';')[0].replace('.v5.1','')
		if type not in GFF:
			GFF[type] = {}
		if chr not in GFF[type]:
			GFF[type][chr] = {}
		if type == 'gene':
			if left not in GFF[type][chr]:
				GFF[type][chr][left] = [dir,right,name]
			else:
				print('%s\t%s\t%s'%(chr,left,name))
				GFF[type][chr][left].append([dir,right,name])
		if type == 'CDS' or type =='five_prime_UTR' or type == 'three_prime_UTR':
			GFF[type][chr][name] = [dir,right,left]
		
def judge(ref,alt):
	T = {}
	for snp in alt.split(','):
		if len(ref) == len(snp) and '*' not in snp:
			T['SNP'] = 1
		else:
			T['indel'] = 1
	return '/'.join(sorted(T.keys()))

def location(chr,pos):
	genic = 0
	if chr in GFF['gene']:
		for left in GFF['gene'][chr]:
			if len(GFF['gene'][chr][left]) >= 3: ### the left is unique to a gene
				if left <= pos and GFF['gene'][chr][left][1] >= pos:
					genic = 1
					gene_name = GFF['gene'][chr][left][2]
					res = []
					for mRNA_name in G[gene_name]: ### if in exonic and intronic regions in different transcripts, then location = 'splicing'
						exonic = 0
						five_UTR = 0
						three_UTR = 0
						for cds_name in G[gene_name][mRNA_name]['CDS']:
							if GFF['CDS'][chr][cds_name][2] <= pos and GFF['CDS'][chr][cds_name][1] >= pos:
								exonic = 1
						if exonic == 0:
							if exonic == 0:
								try:
									for three_name in G[gene_name][mRNA_name]['three_prime_UTR']:
										if GFF['three_prime_UTR'][chr][three_name][2] <= pos and GFF['three_prime_UTR'][chr][three_name][1] >= pos:
											three_UTR = 1
									if three_UTR == 0:
										for five_name in G[gene_name][mRNA_name]['five_prime_UTR']:
											try:
												if GFF['five_prime_UTR'][chr][five_name][2] <= pos and GFF['five_prime_UTR'][chr][five_name][1] >= pos:
													five_UTR = 1
											except:
												print('No five UTR for %s'%gene_name)
								except:
									print('No three UTR for %s'%gene_name)
						res.append([exonic,five_UTR,three_UTR])
					mul = 1
					sum = 0
					sum3 = 0
					sum5 = 0
					for result in res:
						mul = mul*result[0]
						sum = sum + result[0]
						sum3 = sum3 + result[2]
						sum5 = sum5 + result[1]
					if mul == 1:
						location = 'exonic'
					if mul == 0 and sum >= 1:
						location = 'splicing'
					if sum == 0:
						if sum3 >= 1 and sum5 == 0:
							location = 'three_UTR'
						if sum5 >= 1 and sum3 == 0:
							location = 'five_UTR'
						if sum3 >= 1 and sum5 >= 1:
							print('something is wrong with %s_%s'%(chr,pos))
						if sum3 == 0 and sum5 == 0:
							location = 'intronic'
					break
			if len(GFF['gene'][chr][left]) > 3: ### the left is not unique to a gene
				for j in range(3,len(GFF['gene'][chr][left])):
					if left < pos and GFF['gene'][chr][left][j][1] > pos:
						genic = 1
						gene_name = GFF['gene'][chr][left][j][2]
						res = []
						for mRNA_name in G[gene_name]: ### if in exonic and intronic regions in different transcripts, then location = 'splicing'
							exonic = 0
							five_UTR = 0
							three_UTR = 0
							for cds_name in G[gene_name][mRNA_name]['CDS']:
								if GFF['CDS'][chr][cds_name][2] <= pos and GFF['CDS'][chr][cds_name][1] >= pos:
									exonic = 1
							if exonic == 0:
								try:
									for three_name in G[gene_name][mRNA_name]['three_prime_UTR']:
										if GFF['three_prime_UTR'][chr][three_name][2] <= pos and GFF['three_prime_UTR'][chr][three_name][1] >= pos:
											three_UTR = 1
									if three_UTR == 0:
										for five_name in G[gene_name][mRNA_name]['five_prime_UTR']:
											try:
												if GFF['five_prime_UTR'][chr][five_name][2] <= pos and GFF['five_prime_UTR'][chr][five_name][1] >= pos:
													five_UTR = 1
											except:
												print('No five UTR for %s'%gene_name)
								except:
									print('No three UTR for %s'%gene_name)
							res.append([exonic,five_UTR,three_UTR])
						mul = 1
						sum = 0
						sum3 = 0
						sum5 = 0
						for result in res:
							mul = mul*result[0]
							sum = sum + result[0]
							sum3 = sum3 + result[2]
							sum5 = sum5 + result[1]
						if mul == 1:
							location = 'exonic'
						if mul == 0 and sum >= 1:
							location = 'splicing'
						if sum == 0:
							if sum3 >= 1 and sum5 == 0:
								location = 'three_UTR'
							if sum5 >= 1 and sum3 == 0:
								location = 'five_UTR'
							if sum3 >= 1 and sum5 >= 1:
								print('something is wrong with %s_%s'%(chr,pos))
							if sum3 == 0 and sum5 == 0:
								location = 'intronic'
						break
		if genic == 0:
			location = 'intergenic'
	else:
		location = 'intergenic'
	return(location)	
	
	
inp = open(sys.argv[1],'r').readlines()
out = open(sys.argv[1] + '_classification.txt','w')
out.write('Chr\tPos\tRef\tAlt\tType\tallelic\tLocation\n')
for inl in inp[1:]:
	tem = inl.split('\t')
	chr = tem[0]
	pos = int(tem[1])
	ref = tem[2]
	alt = tem[3].strip()
	type = judge(ref,alt)
	allelic = len(alt.split(',')) + 1
	loc = location(chr,pos)
	out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr,pos,ref,alt,type,allelic,loc))

out.close()











