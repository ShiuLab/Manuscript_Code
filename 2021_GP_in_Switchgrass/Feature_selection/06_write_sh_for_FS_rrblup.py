import sys,os
path = sys.argv[1]
trait = sys.argv[2]
start = int(sys.argv[3])
stop = int(sys.argv[4])
step = int(sys.argv[5])
out = open('rrBLUP_FS.sh','w')
out.write('#!/bin/sh --login\n#SBATCH --time=20:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=20G\n#SBATCH --job-name rrBLUP_FS.sh\n#SBATCH -e rrBLUP_FS.sh.e\n#SBATCH -o rrBLUP_FS.sh.o\ncd %s\nmodule load R\nmodule load Python/3.6.4\n'%(path))
for i in range(start,stop,step):
	out.write('Rscript /mnt/home/peipeiw/Documents/Genome_selection/CV/13_rrBLUP_training_test_split_fread.r geno.csv pheno.csv Markers_top%s.txt %s Test.txt 5 10 CVFs.csv Top%s_markers_%s_geno\n'%(i,trait,i,trait))

out.close()