import sys,os
save_name = sys.argv[1]
path_pca = sys.argv[2]
trait = sys.argv[3]
out = open('FS_for_%s.txt'%save_name,'w')
out.write('R2\tType\tNumber_of_markers\n')
for files in os.listdir(path_pca):
	if files.startswith('R2_cv_results_Random_') and files.endswith('_pca.csv') and trait in files:
		inp = open(path_pca + files,'r').readlines()
		for inl in inp[1:]:
			out.write('%s\tpca\t%s\n'%(inl.strip(),files.split('_Random_')[1].split('_')[0]))

for files in os.listdir('./'):
	if files.startswith('R2_cv_results_') and files.endswith('_geno.csv'):
		inp = open(files,'r').readlines()
		for inl in inp[1:]:
			out.write('%s\tgeno\t%s\n'%(inl.strip(),files.split('_Top')[1].split('_')[0]))

out.close()

out = open('FS_for_%s_on_test.txt'%save_name,'w')
out.write('R2\tType\tNumber_of_markers\n')
for files in os.listdir(path_pca):
	if files.startswith('R2_cv_results_Random_') and files.endswith('_pca.csv') and trait in files:
		inp = open(path_pca + files,'r').readlines()
		for inl in inp[1:]:
			out.write('%s\tpca\t%s\n'%(inl.strip(),files.split('_Random_')[1].split('_')[0]))

for files in os.listdir('./'):
	if files.startswith('R2_test_results_') and files.endswith('_geno.csv'):
		inp = open(files,'r').readlines()
		for inl in inp[1:]:
			out.write('%s\tgeno\t%s\n'%(inl.strip(),files.split('_Top')[1].split('_')[0]))

out.close()