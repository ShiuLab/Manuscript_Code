# write jobs for target methylated sites for each accession
import sys,os
path = sys.argv[1]
for files in os.listdir('./'):
	if files.endswith('.tsv'):
		accession = files.split('_')[2].split('.')[0]
		out = open('%s.sh'%accession,'w')
		out.write('#!/bin/sh --login\n#SBATCH --time=1:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=80G\n#SBATCH --job-name %s.sh\n#SBATCH -e %s.sh.e\n#SBATCH -o %s.sh.o\ncd %s\n'%(accession,accession,accession,path))
		out.write('awk \'{print "Chr"$1"_"$2"_"$4"_"$3"\\t"$5"/"$6"\\t"$5}\' < %s > %s_proportion.txt\n'%(files,accession))
		out.write('awk \'$3>0\' < %s_proportion > %s_meted.txt'%(accession, accession))
		out.close()