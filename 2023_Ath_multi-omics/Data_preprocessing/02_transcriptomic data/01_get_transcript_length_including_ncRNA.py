import sys,os
import pprint

### get the TAIR10 GFF3 ready before run this script
inp = open('TAIR10_GFF3_genes.gff','r').readlines()
D = {}
D_mRNA = {}
D_ncRNA = {}
for inl in inp:
	if not inl.startswith('#') and inl.split('\t')[2] in ["mRNA",'CDS','five_prime_UTR','three_prime_UTR','exon','ncRNA','mRNA_TE_gene','tRNA','snoRNA','miRNA','pseudogenic_transcript','pseudogenic_exon',"snRNA",'rRNA']:
		if inl.split('\t')[2] == "mRNA":
			tem = inl.split('\t')
			gene = tem[8].strip().split('Parent=')[1].split(';')[0]
			if gene not in D:
				D[gene] = {}
				D_mRNA[gene] = 1
		elif inl.split('\t')[2] in ["ncRNA","mRNA_TE_gene","tRNA","snoRNA","miRNA","pseudogenic_transcript","snRNA",'rRNA']:
			tem = inl.split('\t')
			gene = tem[8].strip().split('Parent=')[1].split(';')[0]
			if gene not in D:
				D[gene] = {}
				D_ncRNA[gene] = 1
		else:
			tem = inl.split('\t')
			mRNA = tem[8].strip().split('Parent=')[1].split(';')[0].split(',')[0]## name of mRNA
			if mRNA.split('.')[0] in D_mRNA and tem[2] != 'exon':
				if mRNA not in D[gene]:
					D[gene][mRNA] = 0
				D[gene][mRNA] += abs(int(tem[3]) - int(tem[4])) + 1
			if mRNA.split('.')[0] in D_ncRNA:
				if mRNA not in D[gene]:
					D[gene][mRNA] = 0
				D[gene][mRNA] += abs(int(tem[3]) - int(tem[4])) + 1
				
out = open('TAIR10_longest_mRNA_length_including_ncRNA.txt','w')
out.write('mRNA_ID\tLength\n')
for gene in D:
	length = 0
	for mRNA in D[gene]:
		length = max(length,D[gene][mRNA])
	out.write('%s\t%s\n'%(mRNA, length))
	
out.close()
