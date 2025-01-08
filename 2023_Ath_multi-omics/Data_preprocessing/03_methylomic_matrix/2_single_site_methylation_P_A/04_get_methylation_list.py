import sys,os

out = open('Methylation_sites_listgenome_wide_383_accessions.txt','w')
D = {}
for file in os.listdir('./'):
    if file.endswith('_methylated.txt'):
        lines = open(file,'r').readlines()
        for inl in lines[1:]:
            if inl.split('\t')[0] not in D:
                out.write('%s\n'%inl.split('\t')[0])
                D[inl.split('\t')[0]] = 1

out.close()

