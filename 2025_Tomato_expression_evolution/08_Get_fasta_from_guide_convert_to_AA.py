'''
This is a script to take FASTA for different groups of syntenic homelogs and convert them to AA sequnces
Author: LT
'''
from Bio.Seq import Seq
import argparse
import os
#required
parser = argparse.ArgumentParser(description='This code looks at alignmnt guide document and create AA fasta files for given sets of gene sequnces')
parser.add_argument('-f', '--cds',required=True, help="fasa file containing all the coding sequnces in the genome in question")
parser.add_argument('-a', '--aln', required=True, help="alignmnt guide file generated from previous code")
parser.add_argument('-m', '--modifying',required=False, default='n', help="specify if you are using to modify OGs after running re_pick_OGs_after_ML_run.py. Default is n.")
parser.add_argument('-o', '--out', required=True, help="output directory")
args = parser.parse_args()

cds_fas = open(args.cds,"r").readlines()
Cds_dict = {}
alignmnt_guide = open(args.aln,"r").readlines()

for fas in cds_fas:
    if fas.startswith(">"):
        gene_id = fas.replace(">","").lstrip().rstrip()
        if gene_id not in Cds_dict.keys():
            Cds_dict[gene_id] = []
    
    Cds_dict[gene_id] = fas.lstrip().rstrip()
if args.modifying == 'n':    
    for line in alignmnt_guide:
        if not line.startswith("in_group"):
            #print(line)
            ingroups = line.split("\t")[0].lstrip().rstrip().replace('"','').split(",")
            outgrups = line.split("\t")[1].lstrip().rstrip().replace('"','').split(",")
            group = line.split("\t")[2].lstrip().rstrip()
            gene_list = ingroups
            gene_list.extend(outgrups)
            out = open(os.path.join(args.out,f'Group{group}_{outgrups[0]}_{outgrups[1]}.fas'),"w")
            for gene in gene_list:
                gene_id_fas = f'>{gene}'
                seq = Cds_dict[gene]
                dna_seq = Seq(seq)
                aa_seq = dna_seq.translate() 
                #print(aa_seq)
                out.write(f'{gene_id_fas}\n{str(aa_seq)}\n')
            out.close()
elif args.modifying == 'y':
        for line in alignmnt_guide:
            if not line.startswith("in_group"):
                #print(line)
                ingroups = line.split("\t")[0].lstrip().rstrip().replace('"','').split(",")
                outgrups = line.split("\t")[1].lstrip().rstrip().replace('"','').split(",")
                group = line.split("\t")[2].lstrip().rstrip()
                OG_status = line.split("\t")[3].lstrip().rstrip()
                if OG_status == "not_abide":
                    gene_list = ingroups
                    gene_list.extend(outgrups)
                    out = open(os.path.join(args.out,f'Group{group}_{outgrups[0]}_{outgrups[1]}.fas'),"w")
                    for gene in gene_list:
                        gene_id_fas = f'>{gene}'
                        seq = Cds_dict[gene]
                        dna_seq = Seq(seq)
                        aa_seq = dna_seq.translate() 
                        #print(aa_seq)
                        out.write(f'{gene_id_fas}\n{str(aa_seq)}\n')
                    out.close()
        