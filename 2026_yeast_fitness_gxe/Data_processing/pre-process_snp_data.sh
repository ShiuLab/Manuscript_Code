# Description: Steps to pre-process SNP .gvcf file and generate a SNP matrix file

# 1. Filter the phenotype matrix to keep only diploid individuals
# see pre-process_pheno_data.py

# 2. Filter SNPs to keep only biallelic markers with MAF > 0.05
vcftools --vcf 1011Matrix.gvcf --min-alleles 2 --max-alleles 2 \
	--keep diploid_750_IDs.txt --max-missing 0.2 --maf 0.05 \
	--out 1011Matrix_diploid750_biallelic_maf0.05 --recode --remove-indels
# out: 1011Matrix_biallelic_maf0.05.recode.vcf

# 3. Convert filtered .vcf into fastPHASE format
vcftools --vcf 1011Matrix_diploid750_biallelic_maf0.05.recode.vcf --plink \
	--out 1011Matrix_diploid750_biallelic_maf0.05_11032022.plink

module load PLINK/1.9b_4.1-x86_64
plink --file 1011Matrix_diploid750_biallelic_maf0.05.plink --recode12 fastphase \
	--out 1011Matrix_diploid750_biallelic_maf0.05.fastphase # convert to fastPHASE format

plink --file 1011Matrix_diploid750_biallelic_maf0.05.plink --recode12 fastphase \
	--out 1011Matrix_diploid750_biallelic_maf0.05.fastphase --freq # get allele frequencies

# 4. Impute missing genotypes using fastPHASE
./fastPHASE -T10 -H-4 -o1011Matrix_diploid750_biallelic_maf0.05 1011Matrix_diploid750_biallelic_maf0.05.fastphase.chr-0.recode.phase.inp

# 5. Convert the imputed genotypes file into a CSV with the -1,0,1 encoding
# -1: homozygous major allele; 0: heterozygous; 1: homozygous minor allele
path=Data/Peter_2018/0_raw_data
file=1011Matrix_diploid750_biallelic_maf0.05_genotypes.out
map=1011Matrix_diploid750_biallelic_maf0.05.fastphase.frq
recoded=1011Matrix_diploid750_biallelic_maf0.05.fastphase.chr-0.recode.phase.inp
save=1011Matrix_diploid750_biallelic_maf0.05_genotypes.csv # S2 File SNP matrix
python fastPHASE_to_csv.py -path $path -file $file -map $map -recoded $recoded -save $save

# 6. Recode the imputed genotypes file into 0,1,2 format with PLINK1.9
path=Data/Peter_2018/0_raw_data
file=1011Matrix_diploid750_biallelic_maf0.05_genotypes.csv
frq=1011Matrix_diploid750_biallelic_maf0.05.fastphase.frq
vcf=1011Matrix_diploid750_biallelic_maf0.05.recode.vcf
save=geno_012.csv # S12 File SNP matrix (does not include S288C genotypes, yet)
python fastPHASE_to_012csv.py -path $path -file $file -frq $frq -vcf $vcf -save $save
