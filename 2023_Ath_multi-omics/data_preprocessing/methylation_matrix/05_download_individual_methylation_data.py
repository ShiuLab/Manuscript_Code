import sys,os
import os.path
inp = open('Link_for_methylation_383_accessions.txt','r')
'''
for inl in inp:
	try:
		os.system('wget %s'%inl.split('\t')[5])
	except:
		os.system('wget %s'%inl.split('\t')[6])

for inl in inp:
	try:
		os.system('wget %s'%inl.split('\t')[7])
	except:
		os.system('wget %s'%inl.split('\t')[8])
'''
for inl in inp:
	tem = inl.strip().split('\t')
	for link in tem[5:9]:
		if not os.path.isfile(link.split('/')[-1]): 
			os.system('wget %s'%link)
		
