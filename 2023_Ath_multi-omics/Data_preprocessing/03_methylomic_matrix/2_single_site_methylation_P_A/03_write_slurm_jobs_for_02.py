import sys,os
path = sys.argv[1]
for files in os.listdir('./'):
	if files.endswith('.tsv'):
		accession = files.split('_')[2].split('.')[0]
		out = open('%s.sh'%accession,'w')
		out.write('#!/bin/sh --login\n#SBATCH --time=1:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=20G\n#SBATCH --job-name %s.sh\n#SBATCH -e %s.sh.e\n#SBATCH -o %s.sh.o\ncd %s\n'%(accession,accession,accession,path))
		out.write('python 02_only_save_methylated_sites.py %s'%files)
		out.close()
