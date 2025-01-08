import sys,os
import pickle
files = sys.argv[1]
inp = open(files,'r') 
D = {}
inl = inp.readline()
inl = inp.readline()
while inl:
	D[inl.split('\t')[0]] = float(inl.split('\t')[1].split('/')[0])/int(inl.split('\t')[1].split('/')[1])
	inl = inp.readline()
	
f = open("%s.pkl"%files,"wb")
pickle.dump(D,f)
f.close()
