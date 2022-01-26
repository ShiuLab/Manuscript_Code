import sys,os
import pandas as pd
gff = open('Pvirgatum_516_v5.1.gene.gff3','r').readlines()
chrlen = open('Switchgrass_chr_length.txt','r').readlines()
L = {}
for inl in chrlen:
	L[inl.split('\t')[0]] = int(inl.strip().split('\t')[1])

GFF = {}
for inl in gff:
	if not inl.startswith('#'):
		if inl.split('\t')[2] == 'gene':
			tem = inl.split('\t')
			chr = tem[0]
			left = int(tem[3])
			right = int(tem[4])
			dir = tem[6]
			gene = tem[8].strip().split('Name=')[1].split(';')[0]
			if chr not in GFF:
				GFF[chr] = {}
			if left not in GFF[chr]:
				GFF[chr][left] = [gene,right,dir]

Gene = {}
for chr in GFF:
	lefts = sorted(GFF[chr].keys())
	for i in range(0,len(lefts)):
		gene = GFF[chr][lefts[i]][0]
		left = lefts[i]
		right = GFF[chr][lefts[i]][1]
		dir = GFF[chr][lefts[i]][2]
		if gene not in Gene:
			Gene[gene] = {}
		Gene[gene]['chr'] = chr
		Gene[gene]['dir'] = dir
		if dir == '+':
			Gene[gene]['GBL'] = left
			Gene[gene]['GBR'] = right
			Gene[gene]['USR'] = left -1
			Gene[gene]['DSL'] = right +  1
			if i == 0:
				Gene[gene]['USL'] = 1
			else:
				Gene[gene]['USL'] = GFF[chr][lefts[i-1]][1] + 1
			if i == len(lefts) -1:
				Gene[gene]['DSR'] = L[chr]
			else:
				Gene[gene]['DSR'] = lefts[i+1] - 1
		if dir == '-':
			Gene[gene]['GBL'] = right
			Gene[gene]['GBR'] = left
			Gene[gene]['DSL'] = left -1
			Gene[gene]['USR'] = right +  1
			if i == 0:
				Gene[gene]['DSR'] = 1
			else:
				Gene[gene]['DSR'] = GFF[chr][lefts[i-1]][1] + 1
			if i == len(lefts) -1:
				Gene[gene]['USL'] = L[chr]
			else:
				Gene[gene]['USL'] = lefts[i+1] - 1

Gene_loc = pd.DataFrame.from_dict(Gene, orient='index')
GBS = pd.read_csv('All_filtered_SNPs_indels_20181203_GBS_4col.txt_classification.txt',header=0,index_col=None,sep='\t')
#GBS = GBS.loc[GBS['Type']=='SNP',:]
#GBS = GBS.loc[GBS['allelic']==2,:]

GBS_markers = open('Markers_distribution_all_markers.txt','w')
for gene in Gene_loc.index.tolist():
	tem = GBS.loc[GBS['Chr']==Gene_loc.loc[gene,'chr']]
	GB = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'GBL']) * (tem['Pos'] - Gene_loc.loc[gene,'GBR']) <= 0,]
	UP = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'USL']) * (tem['Pos'] - Gene_loc.loc[gene,'USR']) <= 0,]
	Down = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'DSL']) * (tem['Pos'] - Gene_loc.loc[gene,'DSR']) <= 0,]
	if GB.shape[0] > 0:
		for pos in GB['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('GBS\t%s\t%s\t+\tGene_body\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'GBL'],Gene_loc.loc[gene,'GBR'],pos,(pos - Gene_loc.loc[gene,'GBL'] + 1)/float(Gene_loc.loc[gene,'GBR'] - Gene_loc.loc[gene,'GBL'] + 1)*3500 ))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('GBS\t%s\t%s\t-\tGene_body\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'GBL'],Gene_loc.loc[gene,'GBR'],pos,(Gene_loc.loc[gene,'GBL'] - pos + 1)/float(Gene_loc.loc[gene,'GBL'] - Gene_loc.loc[gene,'GBR'] + 1)*3500 ))
	if UP.shape[0] > 0:
		for pos in UP['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('GBS\t%s\t%s\t+\tUpstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos,pos - Gene_loc.loc[gene,'USR']))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('GBS\t%s\t%s\t-\tUpstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos, Gene_loc.loc[gene,'USR'] - pos))
	if Down.shape[0] > 0:
		for pos in Down['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('GBS\t%s\t%s\t+\tDownstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos,3500 + pos - Gene_loc.loc[gene,'DSL']))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('GBS\t%s\t%s\t-\tDownstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos, 3500 + Gene_loc.loc[gene,'DSL'] - pos))

GBS = pd.read_csv('All_filtered_SNPs_indels_20181203_exome_capture_4col.txt_classification.txt',header=0,index_col=None,sep='\t')
#GBS = GBS.loc[GBS['Type']=='SNP',:]
#GBS = GBS.loc[GBS['allelic']==2,:]
for gene in Gene_loc.index.tolist():
	tem = GBS.loc[GBS['Chr']==Gene_loc.loc[gene,'chr']]
	GB = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'GBL']) * (tem['Pos'] - Gene_loc.loc[gene,'GBR']) <= 0,]
	UP = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'USL']) * (tem['Pos'] - Gene_loc.loc[gene,'USR']) <= 0,]
	Down = tem.loc[(tem['Pos'] - Gene_loc.loc[gene,'DSL']) * (tem['Pos'] - Gene_loc.loc[gene,'DSR']) <= 0,]
	if GB.shape[0] > 0:
		for pos in GB['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('Exome\t%s\t%s\t+\tGene_body\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'GBL'],Gene_loc.loc[gene,'GBR'],pos,(pos - Gene_loc.loc[gene,'GBL'] + 1)/float(Gene_loc.loc[gene,'GBR'] - Gene_loc.loc[gene,'GBL'] + 1)*3500 ))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('Exome\t%s\t%s\t-\tGene_body\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'GBL'],Gene_loc.loc[gene,'GBR'],pos,(Gene_loc.loc[gene,'GBL'] - pos + 1)/float(Gene_loc.loc[gene,'GBL'] - Gene_loc.loc[gene,'GBR'] + 1)*3500 ))
	if UP.shape[0] > 0:
		for pos in UP['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('Exome\t%s\t%s\t+\tUpstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos,pos - Gene_loc.loc[gene,'USR']))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('Exome\t%s\t%s\t-\tUpstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos, Gene_loc.loc[gene,'USR'] - pos))
	if Down.shape[0] > 0:
		for pos in Down['Pos']:
			if Gene_loc.loc[gene,'dir']=='+':
				GBS_markers.write('Exome\t%s\t%s\t+\tDownstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos,3500 + pos - Gene_loc.loc[gene,'DSL']))
			if Gene_loc.loc[gene,'dir']=='-':
				GBS_markers.write('Exome\t%s\t%s\t-\tDownstream\t%s\t%s\t%s\t%s\n'%(Gene_loc.loc[gene,'chr'],gene,Gene_loc.loc[gene,'USL'],Gene_loc.loc[gene,'USR'],pos, 3500 + Gene_loc.loc[gene,'DSL'] - pos))

GBS_markers.close()

df = pd.read_csv('Markers_distribution_all_markers.txt',header=None,index_col=None,sep='\t')
df['Marker'] = df.apply(lambda row: '%s_%s'%(row[1],row[7]), axis = 1)
df['distance'] = abs(df[8])
df.loc[df[4]=='Downstream','distance'] = abs(df.loc[df[4]=='Downstream',8] - 3500)
df['Type'] = 0
df.loc[df[4]=='Upstream','Type'] = 1
df.loc[df[4]=='Downstream','Type'] = 1
GBS = df[df[0]=='GBS']
Exome = df[df[0]=='Exome']
GBS = GBS.sort_values(['Type','distance'], ascending=True).drop_duplicates('Marker',keep="first")
Exome = Exome.sort_values(['Type','distance'], ascending=True).drop_duplicates('Marker',keep="first")
df2 = pd.concat([GBS,Exome],axis=0)
df2.to_csv('Markers_distribution_all_markers_unique.txt',header=False,index=False,sep='\t')