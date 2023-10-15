import sys
inp = open(sys.argv[1],"r").readlines() # the SNP matrix
out = open(sys.argv[2],"w") # output the hmp format for kinship calculation

def list2string(s):
	newstr = "" 
	for ele in s:
		if ele == '0':
			ele = 'AA'
		if ele == '1':
			ele = 'TT'
		newstr += ele + "\t" 
	return newstr

i = 0
for line in inp:
	li = line.strip().split(',')
	if i == 0:
		rs = "rs#"
		alleles = "alleles"
		chrom = "chrom"
		pos = "pos"
		strand = "strand"
		assembly = "assembly#"
		center = 'center'
		protLSID = 'protLSID'
		assayLSID = "assayLSID"
		panel = "panel"
		QCcode = "QCcode"
		genotype = list2string(li[1:])	
		result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(rs,alleles,chrom,pos,strand,assembly,center,protLSID,assayLSID,panel,QCcode,genotype)
		out.write(result)
		i = i+1
	else:
		rs = "Ath_"+str(i)
		alleles = "A/T"
		chrom = li[0].split('_')[0]
		pos = li[0].split('_')[1]
		strand = "+"
		assembly = "NA"
		center = 'NA'
		protLSID = 'NA'
		assayLSID = "NA"
		panel = "NA"
		QCcode = "NA"
		genotype = list2string(li[1:])
		result = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(rs,alleles,chrom,pos,strand,assembly,center,protLSID,assayLSID,panel,QCcode,genotype)
		out.write(result)
		i = i+1

out.close()