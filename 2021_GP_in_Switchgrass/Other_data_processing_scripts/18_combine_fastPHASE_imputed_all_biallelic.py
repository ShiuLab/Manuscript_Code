import sys,os
number = int(sys.argv[1])
out = open('Octa_exome_biallelic_SNPs_imputed_combined.txt','w')
D = {}
for n in range(1,number//5000+1):
	impute = open('Octa_exome_biallelic_SNP_%s_%s_hapguess_switch.out'%((n-1)*5000+1, n*5000),'r').readlines()
	x = 1
	while x < len(impute):
		inl = impute[x]
		if inl.startswith('# id'):
			ind = inl.split('# id ')[1].strip()
			x += 1
			a1 = impute[x].strip()
			x += 1
			a2 = impute[x].strip()
			if ind not in D:
				D[ind] = [a1,a2]
			else:
				D[ind][0] = D[ind][0] + ' ' + a1
				D[ind][1] = D[ind][1] + ' ' + a2
		x += 1
			
if 1==1:
	impute = open('Octa_exome_biallelic_SNP_%s_%s_hapguess_switch.out'%((number//5000) * 5000 + 1, number),'r').readlines()
	x = 1
	while x < len(impute):
		inl = impute[x]
		if inl.startswith('# id'):
			ind = inl.split('# id ')[1].strip()
			x += 1
			a1 = impute[x].strip()
			x += 1
			a2 = impute[x].strip()
			D[ind][0] = D[ind][0] + ' ' + a1
			D[ind][1] = D[ind][1] + ' ' + a2
		x += 1

for ind in D:
	out.write('# id %s\n'%ind)
	out.write(D[ind][0] + '\n')
	out.write(D[ind][1] + '\n')

out.close()
