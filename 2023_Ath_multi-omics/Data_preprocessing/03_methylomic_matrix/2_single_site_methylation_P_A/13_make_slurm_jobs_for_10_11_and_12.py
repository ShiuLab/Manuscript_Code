import sys,os
path = sys.argv[1]
for i in range(1,101):
	files = 'Methylation_genome_wide_383_accessions.csv_%s'%i
	out = open('Job16_%s.sh'%i,'w')
	out.write('#!/bin/sh --login\n#SBATCH --time=4:00:00\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=80G\n#SBATCH --job-name Job16_%s.sh\n#SBATCH -e Job16_%s.sh.e\n#SBATCH -o Job16_%s.sh.o\ncd %s\n'%(i,i,i,path))
	out.write('module load Python/3.6.4\n')
	out.write('python 10_fill_0_with_NaN.py %s\n'%(files))
	out.write('python 11_fill_NA_back_with_0_for_shared_sites.py %s\n'%(files))
	out.write('python 12_check_methylation_data_NaN_proportions_and_drop_separate_files.py %s_transposed_filled_with_NaN_considering_same_site.csv\n'%(files))
	out.close()