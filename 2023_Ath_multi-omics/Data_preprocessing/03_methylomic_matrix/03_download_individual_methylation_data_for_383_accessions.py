import sys,os
import os.path

# the input file can be found in the folder "Datasets"
inp = open('Link_for_methylation_383_accessions.txt','r')

for inl in inp:
	try:
		os.system('wget %s'%inl.split('\t')[5])
	except:
		try:
            os.system('wget %s'%inl.split('\t')[6])
		except:
            try:
                os.system('wget %s'%inl.split('\t')[7])
            except:
                os.system('wget %s'%inl.split('\t')[8])
           
# for the ones not downloaded
for inl in inp:
	tem = inl.strip().split('\t')
	for link in tem[5:9]:
		if not os.path.isfile(link.split('/')[-1]): 
			os.system('wget %s'%link)
		
