import sys,os
import pandas as pd
n = 0
for files in os.listdir('./'):
	if files.endswith('filtered'):
		if n == 0:
			df = pd.read_csv(files, sep='\t', index_col = 0, header = 0)
			n += 1
		else:
			df2 = pd.read_csv(files, sep='\t', index_col = 0, header = 0)
			df = pd.concat([df,df2],axis=0)

df.to_csv('All_filtered_SNPs_indels_20181203_exome_capture.txt', index=True, header=True,sep="\t")


inp = open('All_filtered_SNPs_indels_20181203_exome_capture.txt','r')
out = open('All_filtered_SNPs_indels_20181203_exome_capture_4col.txt','w')
inl = inp.readline()
while inl:
	out.write('\t'.join(inl.split('\t')[0:4]) + '\n')
	inl = inp.readline()
	
out.close()