'''
This is a script to take FASTA made from previous code (get_fasta_from_guide_convert_to_AA.py or get_fasta_from_guid.py) and align them using MAFFT and then construct boostrapped phylogenies using RAxML
Author: LT
'''
import os
import numpy as np
import random
import argparse
#required
parser = argparse.ArgumentParser(description='This code writes and submit a shell script to align cds using MAFFT and construct boostrpped phylogenies using RAxML')
parser.add_argument('-p', '--path',required=True, help="Path of the directory containing alignmnt files")
parser.add_argument('-pf', '--prefix', required=True, help="prefix to identify the alignment file")
#not requred
parser.add_argument('-t', '--time', type=int, required=False, default=1,help="Time in hours for slurm submission")
args = parser.parse_args()

wd = args.path
suffix = args.prefix
time = args.time

os.chdir(wd)
for file in os.listdir(wd):
    if (file.startswith(suffix)) and (file.endswith(".fas")) and (not file.endswith("_aligned.fas")):
        #print(f'{file} found')
        og = ",".join(file.replace(".fas","").split("_")[1:])
        group = file.split("_")[0].lstrip().rstrip()
        out_aln_nm = file.replace(".fas","_aligned.fas")
        out_tree_nm = file.replace(".fas","_RAxML_1000.out")
        cmd = open(f'{group}.sh',"w")
        cmd.write(f"#!/bin/bash\n########## BATCH Lines for Resource Request ##########\n#SBATCH --time={time}:00:00\n#SBATCH --nodes=1\n#SBATCH --ntasks=2\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=80G\n#SBATCH --job-name {group}\n########## Command Lines to Run #########\nmodule purge\nmodule load GCC/10.2.0 OpenMPI/4.0.5 MAFFT\ncd ./\n")
        cmd.write(f"mafft --genafpair --maxiterate 1000 {file} > {out_aln_nm}\n")
        cmd.write(f"module purge\nmodule load icc/2019.1.144-GCC-8.2.0-2.31.1 impi/2018.4.274 RAxML/8.2.12-hybrid-avx2\n ")
        cmd.write(f"raxmlHPC -n {out_tree_nm} -f a -x {np.random.choice(range(1,20000),1,replace=False)[0]} -T 4 -p {np.random.choice(range(1,20000),1,replace=False)[0]} -# 1000 -m PROTGAMMAAUTO –autoprot=bic -s {out_aln_nm} -o {og}\n")
        cmd.close()
        os.system(f"sbatch {group}.sh")
     