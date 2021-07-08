# Genomic_prediction_in_Switchgrass

## Pre-processing
Note that when using python script in HPCC, please module load Python/3.6.4

### Step 1. convert vcf file to genetic matrix

	python 01_conver_genotype_gvcf_to_genotype_matrix.py -file your_vcf
 
### Step 2. filter the genotype matrix

	python 02_filter_genotype_matrix_MAF_missing_data.py -file genotype_matrix
 
Note that if you only want to get the biallelic SNPs or indels, rather than the classification of markers to genic, intergenic, etc, please skip step 3 and 4, and try the script below:

	python 03_get_biallelic_markers_directly.py -file 1011Matrix_genotype_matrix.txt_filtered -type SNP

### Step 3. classify the variation into SNP, indel, or SNP/indel; biallelic or non-biallelic; in genic or intergenic, three_UTR or five_UTR region, exonic or intronic, splicing regions

Note that, this script is for switchgrass specifically. For other species, be careful about the gff format and the how the gene names are encoded**

	python 03_classify_variations.py -file genotype_matrix_filtered -gff gff_file

### Step 4. extract the biallelic SNPs or indels

	python 04_get_bi-allelic_SNP_or_indel.py -classification marker_classification -file genotype_matrix_filtered -type SNP_or_indel

 
### Step 5. convert the genotype matrix to the fastPHASE format

	python 05_convert_genotype_matrix_to_fastPHASE_format.py -file 1011Matrix_genotype_matrix.txt_filtered_biallelic_SNP.txt
 
### Step 6. impute the missing data using fastPHASE, please download and install the [fastPHASE](http://scheet.org/software.html) first

	./fastPHASE -T10 -oName_for_output 1011Matrix_genotype_matrix.txt_filtered_biallelic_SNP.txt_fastPHASE.txt

### Step 7. convert the imputed genotype matrix back to the format used previously

	python 06_convert_imputed_biallelic_variation_to_genotype.py -matrix 1011Matrix_genotype_matrix.txt_filtered_biallelic_SNP.txt -imputed_matrix Name_for_output_hapguess_switch.out
 
 
## Now you can build the genomic prediction models, using the rrBLUP

Note that make sure you have the geno and pheno matrices beforehand.

### Step 8. make CVs file, which will be used for the cross-validation scheme, here 5-fold cross-validation scheme is repeated 10 times

	python 07_make_CVs.py -file pheno.csv -cv 5 -number 10

### Step 9. get the population structure, which is defined as the top 5 pricinple components from the genetic markers

	Rscript 08_getPCs.r geno.csv pheno.csv
 
Note that if you have too large geno matrix, try the following two steps to get the top five PCs

### Step 9_1. remove constant columns in the geno matrix

	python 08_1_remove_constant_columns.py geno.csv
 
### Step 9_2. get the top 5 PCs

	Rscript 08_2_PCA_after_removing.r geno_non_constant.csv

The logical for the following script is that: for each cross-validation fold, using the training fold to build a model, then apply the model to the validation fold. So now you have the predicted values for individuals in the  validation fold. After run for each of the CV fold, you would have the predicted values for all your individuals. Finally, the r2 was calculated using the true and predicted values of all your individuals. This will be repeated n times as you set and n r2 values will be reported.

### Step 10. genomic prediction using the genetic markers or population structure within a cross-validation scheme

Note: 09_rrBLUP_fread_PCs.r builds a simple linear regression model using the PCs, which has almost the same results of models built using mixed.solve from rrBLUP

	Rscript 09_rrBLUP_fread_predict_values.r geno.csv pheno.csv all all 5 10 CVFs.csv exome_geno
	
	Rscript 09_rrBLUP_fread_predict_values.r PCA5_geno.csv pheno.csv all all 5 10 CVFs.csv exome_pca
	
or

	Rscript 09_rrBLUP_fread_PCs.r PCA5_geno.csv pheno.csv all all 5 10 CVFs.csv exome_pca

Steps 11-14 are how I did the feature selection. Markers were selected using the training set, and in the end, models built using selected markers were applied on the test set. Since I did 5-fold CV, thus 1/6 of individuals will be used as test set, and the remaining 5/6 individuals will be used in the 5-fold cv.

### Step 11. hold out individuals for test set, do the stratified sampling

	python 10_holdout_test_stratified.py pheno.csv target_trait 6

### Step 12. if you want to do feature selection, then you should build models without the test set. So first, get the matrix for the training set, make the CVs file using the training individuals, then build models using the training matrices and output the coef of markers

	Rscript 11_split_geno_pheno.r geno.csv pheno.csv Test.txt
	
	python 07_make_CVs.py -file pheno_training.csv -cv 5 -number 10
	
	Rscript 09_rrBLUP.r geno_training.csv pheno_training.csv all target_trait 5 10 CVFs.csv exome_geno

### Step 13. select the number of markers based on the abs coef

	python 12_select_markers_according_to_abs_coef.py  -coef coef_file -start 250 -stop 5250 -step 250

### Step 14. genomic prediction using the genetic markers or population structure within a cross-validation scheme

The logical for the following script is that: first X% of all the individuals will be held out as test set, which will never be used in the model training process, the remaining 1-X% will be used to train the model, using exactly the same approach as step 10. For each cv fold, the model was also applied to the test set, and after run for all the cv folds, the average r2 across n cv folds will be used for the test set. 

	Rscript 13_rrBLUP_training_test_split_fread_predict_values.r geno.csv pheno.csv selected_markers target_trait Test.txt 5 10 CVFs.csv selected_markers_geno

### Step 15. get the prediction using the top 5 PCs for ramdomly selected markers. Make sure the title of the first column in your geno matrix is "ID".

	Rscript 14_random_select_subset.r geno_file start stop step total_number
	
	Rscript 15_rrBLUP_pca_for_subset_markers.r geno_250.csv pheno.csv selected_markers target_trait Test.txt 5 10 CVFs.csv Random_250_markers_pca
 
