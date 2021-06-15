import sys,os
path = sys.argv[1]
trait = sys.argv[2]
start = int(sys.argv[3])
stop = int(sys.argv[4])
step = int(sys.argv[5])
out = open('rrBLUP_pca_%s.sh'%trait,'w')
out.write('#!/bin/sh --login\n#SBATCH --time=4:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=20G\n#SBATCH --job-name rrBLUP_pca_%s.sh\n#SBATCH -e rrBLUP_pca_%s.sh.e\n#SBATCH -o rrBLUP_pca_%s.sh.o\ncd %s\nmodule load R\n'%(path,trait,trait,trait))
for i in range(start,stop,step):
	out.write('Rscript /mnt/home/peipeiw/Documents/Genome_selection/CV/14_rrBLUP_pca_for_subset_markers.r geno_%s.csv pheno.csv all %s Test.txt 5 10 CVFs.csv Random_%s_markers_%s_pca\n'%(i,trait,i,trait))

out.close()
