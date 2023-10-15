import sys,os
out = open('Job_for_90per.sh','w')
for i in range(2,101):
	out.write("cut -d ',' -f2- Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_90per_imputed_MAF > Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_90per_imputed_MAF_cut\n"%(i,i))

out.write("paste Methylation_genome_wide_383_accessions.csv_1_transposed_filled_with_NaN_considering_same_site.csv_90per_imputed_MAF *_90per_imputed_MAF_cut -d ',' > Methylation_genome_wide_383_accessions_90per_imputed_MAF.csv")
out.close()


out = open('Job_for_75per.sh','w')
for i in range(2,101):
	out.write("cut -d ',' -f2- Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_75per_imputed_MAF > Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_75per_imputed_MAF_cut\n"%(i,i))

out.write("paste Methylation_genome_wide_383_accessions.csv_1_transposed_filled_with_NaN_considering_same_site.csv_75per_imputed_MAF *_75per_imputed_MAF_cut -d ',' > Methylation_genome_wide_383_accessions_75per_imputed_MAF.csv")
out.close()



out = open('Job_for_50per.sh','w')
for i in range(2,101):
	out.write("cut -d ',' -f2- Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_50per_imputed_MAF > Methylation_genome_wide_383_accessions.csv_%s_transposed_filled_with_NaN_considering_same_site.csv_50per_imputed_MAF_cut\n"%(i,i))

out.write("paste Methylation_genome_wide_383_accessions.csv_1_transposed_filled_with_NaN_considering_same_site.csv_50per_imputed_MAF *_50per_imputed_MAF_cut -d ',' > Methylation_genome_wide_383_accessions_50per_imputed_MAF.csv")
out.close()
