import sys,os
import pickle
files = sys.argv[1]
inp = open(files,'r') # *proportion.txt or *_meted.txt files from 06_write_slurm_jobs_for_methylation_proportion_and_save_methylated_site.py
D = {}
inl = inp.readline()
inl = inp.readline()
while inl:
	D[inl.split('\t')[0]] = float(inl.split('\t')[1].split('/')[0])/int(inl.split('\t')[1].split('/')[1])
	inl = inp.readline()
	
f = open("%s.pkl"%files,"wb")
pickle.dump(D,f)
f.close()
