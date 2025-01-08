import sys,os

S = {}
seq = open('TAIR10_chr1.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chr2.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chr3.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chr4.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chr5.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chrC.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

seq = open('TAIR10_chrM.fas','r').readlines()
chr = ''
for inl in seq:
	if inl.startswith('>'):
		chr = inl.split()[0][1:]
		seq = ''
	else:
		seq = seq + inl.strip()

S[chr] = seq

import pickle
f = open("Chromosome_seq_TAIR10.pkl","wb")
pickle.dump(S,f)
f.close()

# load dictionary
import pickle
with open('Chromosome_seq_TAIR10.pkl', 'rb') as f1:
	S = pickle.load(f1)

# make a new column from existing columns

R = {'A':'T','T':'A','C':'G','G':'C'}

inp = open('Methylation_genome_wide_383_accessions_overlapping_with_SNP.csv','r').readlines()[1:]
out =  open('Methylation_genome_wide_383_accessions_overlapping_with_SNP_ref_seq.txt','w')
for inl in inp:
	chr = inl.split(',')[0].split('_')[0]
	pos = int(inl.split(',')[0].split('_')[1])
	dir = inl.strip().split(',')[1][-1]
	if dir == '+':
		nuc = S[chr][pos]
	else:
		nuc = R[S[chr][pos]]
	out.write('%s,%s\n'%(inl.strip(),nuc))

out.close()

# only save the one need to be saved
ref = open('Methylation_genome_wide_383_accessions_overlapping_with_SNP_ref_seq.txt','r').readlines()
Ref = {}
for inl in ref:
	tem = inl.strip().split(',')
	if tem[0] not in Ref:
		Ref[tem[0]] = {}
	if tem[1] not in Ref[tem[0]]:
		Ref[tem[0]][tem[1]] = tem[2]

f = open("Methylation_genome_wide_383_accessions_overlapping_with_SNP_ref_seq.pkl","wb")
pickle.dump(Ref,f)
f.close()

import pickle
with open('Methylation_genome_wide_383_accessions_overlapping_with_SNP_ref_seq.pkl', 'rb') as f1:
	Ref = pickle.load(f1)

SNP = open('/mnt/home/peipeiw/Documents/Ath_GS/datasets/SNP_binary_matrix_383_accessions.csv','r')
S = {}
L = {}
inl = SNP.readline()
tem = inl.strip().split(',')
for i in range(1,len(tem)):
	L[i] = tem[i]

Result = {}
inl = SNP.readline()
while inl:
	tem = inl.strip().split(',')
	if tem[0] in Ref:
		for met in Ref[tem[0]]:
			if Ref[tem[0]][met] == 'C':
				for i in range(1,len(tem)):
					if tem[i] == '1':
						if met not in Result:
							Result[met] = []
						Result[met].append(L[i])
			if Ref[tem[0]][met] != 'C':
				for i in range(1,len(tem)):
					if tem[i] == '0':
						if met not in Result:
							Result[met] = []
						Result[met].append(L[i])
	inl = SNP.readline()

f = open("Methylation_genome_wide_383_accessions_overlapping_with_SNP_targeted.pkl","wb")
pickle.dump(Result,f)
f.close()


out = open('Methylation_genome_wide_383_accessions_overlapping_with_SNP_targeted','w')
for met in Result:
	out.write('%s\t%s\n'%(met,','.join(Result[met])))

out.close()









