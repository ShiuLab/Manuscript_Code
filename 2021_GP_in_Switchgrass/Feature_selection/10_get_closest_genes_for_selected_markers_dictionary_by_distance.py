import sys,os
import pickle
with open('/mnt/home/peipeiw/Documents/Genome_selection/Other_markers/Intermediate_files/Corresponding_marker_locus_for_multiallelic.pkl', 'rb') as f1:
	Corr = pickle.load(f1)


file = sys.argv[1]
distance = int(sys.argv[2])
marker = open(file,'r').readlines()
inp = open('/mnt/home/peipeiw/Documents/Genome_selection/Distribution_of_markers/Markers_distribution_all_markers_unique.txt','r').readlines()
D = {}
for inl in inp:
	tem = inl.strip().split('\t')
	D[tem[9]] = [tem[2],tem[4],float(tem[8])]

R = {}
#out = open('Genes_closest_to_%s_%s'%(file,distance),'w')
out = open('Genes_closest_to_%s_%s_info.txt'%(file,distance),'w')
for inl in marker:
	tem = inl.strip().split('_')
	if 'scaffold' not in inl:
		m = '_'.join(tem[3:5])
		m2 = '_'.join(tem[3:])
		if m in D:
			if D[m][1] == 'Downstream':
				if D[m][2]-3500 <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Upstream':
				if abs(D[m][2]) <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Gene_body':
				out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
		else:
			m = '_'.join(Corr[m2].split('_')[0:2])
			if D[m][1] == 'Downstream':
				if D[m][2]-3500 <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Upstream':
				if abs(D[m][2]) <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Gene_body':
				out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
	else:
		m = '_'.join(tem[3:6])
		m2 = '_'.join(tem[3:])
		if m in D:
			if D[m][1] == 'Downstream':
				if D[m][2]-3500 <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Upstream':
				if abs(D[m][2]) <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Gene_body':
				out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
		else:
			m = '_'.join(Corr[m2].split('_')[0:3])
			if D[m][1] == 'Downstream':
				if D[m][2]-3500 <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Upstream':
				if abs(D[m][2]) <= distance:
					out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
			if D[m][1] == 'Gene_body':
				out.write('%s\t%s\t%s\t%s\n'%(D[m][0],D[m][1],D[m][2],m2))
	out.flush()

out.close()
