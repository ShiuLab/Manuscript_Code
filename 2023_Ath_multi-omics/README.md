# multi-omics integration project

Codes and datasets for our manuscript "Prediction of plant complex traits via integration of multi-omics data".

# **1. Datasets**

## 1.1 SNP matrix

The SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz

## 1.2 transcriptomic data

The transcriptomic data was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz
which are read count files for 727 accessions

## 1.3 methylomic data

The gene body methylation data was download from http://signal-genet.salk.edu/1001.php , which contains 1107 methylomes, with the ID name.
The tsv files for individual accessions were download from NCBI with GEO accession ID: GSE43857

## 1.4 phenotypic data

Flowering time at 10℃ (1163 accessions) and 16℃ (1123 accessions) were downloaded from https://arapheno.1001genomes.org/study/12/
Other phenotypic data were download from https://arapheno.1001genomes.org/study/38/

## 1.5 benchmark flowering time genes

Benchmark flowering time genes were downloaded from the FLOR-ID database http://www.phytosystems.ulg.ac.be/florid/

# **2. Data preprocessing**

## **2.1 SNP matrix**

Related scripts for SNP matrix can be found in the folder Data_preprocessing/01_SNP_matrix

>Get the SNP matrix for all accessions, and also output the list of all SNPs

```  
python 01_extract_the_h5py_binary_SNP_matrix.py
```
	
>Save the SNP matrix for 383 accession; convert the common allele to 1, and the alternative allele to -1; get rid of the rare variants (MAF < 5%)	
  
```
python 02_get_snp_matrix_for_383_accessions.py
```
	
## **2.2 transcriptomic data**

Related scripts can be found in the folder Data_preprocessing/02_transcriptomic data 
  
>To calculate the TPM (transcripts per kilobase million), we need to get the transcript length for genes at first. Get the GFF3 file ready before run this script. Here we used the TAIR10 GFF3 file, which is downloaded from the TAIR database 

```  
python 01_get_transcript_length_including_ncRNA.py
```

>Calculate TPM using the function “calculateTPM” from the R package “scater”

>Put the downloaded file GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv and the transcript length file TAIR10_longest_mRNA_length_including_ncRNA.txt in your work directory
  
```
Rscript 02_calculate_TPM.r
```

>Log the TPM and get the logged TPM matrix for 383 accessions
  
```
python 03_log_TPM_and_get_transcriptomic_matrix_for_383_accessions.py
```	
	
## **2.3 methylomic data**

### **2.3.1 For gene-body methylation**

Related scripts can be found in the folder Data_preprocessing/03_methylomic_matrix/1_gene_body_methylation

>Parse the downloaded gene-body methylation matrix, and save the corresponding matrix for 383 accessions

```
python 01_get_methylation_for_383_accessions.py
```

>Impute missing data
 
```
python 02_knn_imputation.py input_file Test_for_383_accessions.txt
```

### **2.3.2 For single-site methylation-based formats**
	
>Download individual methylation data for 383 accessions 

```
python 01_download_individual_methylation_data_for_383_accessions.py
```
	
#### 2.3.2.1 For presence/absence of methylation

Related scripts can be found in the folder Data_preprocessing/03_methylomic_matrix/2_single_site_methylation_P_A
	
>Only save methylated sites, which are sites with the sixth column at 1 in the downloaded *.tsv file. 

```
python 02_only_save_methylated_sites.py inputFile
```

>>To make slurm jobs for the above script, you can use the script below, your_work_dir is the directory containing all the downloaded *.tsv files

```
python 03_write_slurm_jobs_for_03.py inputFile
```

>Get the methylated single site list for all 383 accessions, and order the list

```
python 04_get_methylation_list.py
```
	
```
Rscript 05_order_methylation_sites.r
```

>Combine all files with the single-site presence/absence state of methylation, and the output file is Methylation_genome_wide_383_accessions.csv

```
python 06_combine_all_methylation_files.py
```

>To distinguish the un-methylated sites (value should be 0) from sites with missing information (values should be NA to be imputed in the future), we first converted all 0 values to NAs, and then filled the NAs using information from SNP matrix. 

>>Get targeted methylated sites (all single sites which were methylated in at least one accession) for each accession. The input files are the downloaded *.tsv files, and the file Methylation_sites_listgenome_wide_383_accessions.txt is output from python 04_get_methylation_list.py

```
awk \'{print "Chr"$1"_"$2"_"$4"_"$3}\' < inputFile > inputFile_selected_columns
```

```
awk \'NR==FNR { lines[$0]=1; next } $0 in lines\' inputFile_selected_columns Methylation_sites_listgenome_wide_383_accessions.txt > inputFile_met.txt
```

>>To make slurm jobs for the above two awk commands, you can use the script below. Your_work_dir is the directory containing all the downloaded *.tsv files

```
python 07_get_targeted_sites.py your_work_dir
```

>Get overlapping sites between the SNPs and the methylation sites

```
Rscript 08_get_overlapping_site_between_M_and_G.r
```

>Get the SNP information for the overlapping sites betweem SNPs and methylation sites

```
python 09_get_ref_seq_for_overlapping_methylation_sites.py inputFile
```


>>To facilitate the job running, the Methylation_genome_wide_383_accessions.csv file was split into 100 small files

```
awk 'NR==1{header=$0; count=1; print header > "Methylation_genome_wide_383_accessions.csv_" count; next } !( (NR-1) % 173766){		count++; print header > "Methylation_genome_wide_383_accessions.csv_" count; }  {print $0 > "Methylation_genome_wide_383_accessions.csv_" count	 }' Methylation_genome_wide_383_accessions.csv
```

>Fill 0s with NAs, and then fill the NAs back with 0s for overlapping site between SNPs and methylation sites. The inputFile is one of the small files output from the above awk command.

```
python 10_fill_0_with_NaN.py inputFile
```

```
python 11_fill_NA_back_with_0_for_shared_sites.py inputFile
```

>Count the number of accessions which have NAs for each methylation site, and the inputFile is one of the small files

```
python 12_check_methylation_data_NaN_proportions_and_drop_separate_files.py inputFile
```

>>To make slurm jobs for the above three python scripts, you can use the script below. your_work_dir is the directory containing all the small files

```
python 13_make_slurm_jobs_for_10_11_and_12.py your_work_dir
```

>Concat all count files together

```
cat *count > Methylation_genome_wide_618_accessions_considering_same_site_NaN_count.txt
```

>Find out the number of accessions with missing methylation data for sites at 90th, 75th and 50th percentiles. This is done in R

```
dat <- read.table('Methylation_genome_wide_618_accessions_considering_same_site_NaN_count.txt',head=F,sep='\t')
dat <- cbind(dat,618-dat[,2])
quantile(dat[,3], c(.50, .75, .90))
```

>Distinguish methylation sites that are missing in <=2 (90th percentile), <=8 (75th percentile), <=38 (50th percentile) accessions

```
python 14_keep_50_75_90_percentile_MAF.py 
```

>Impute the 50th percentile files; fit on the training dataset, then transform it on the test dataset; filter with MAF (minor allele frequency); extract the 75th percentile and 90th percentile informations from these 50th percentile files. The file Test_for_383_accessions.txt can be found in the folder /Datasets

```
python 15_knn_imputation_for_0_1_training_fit_on_test.py 50per_file Test_for_383_accessions.txt
```

>>Make slurm jobs to concat the final imputed methylation files, using "cut" and "paste" shell commands

```
python 16_write_jobs_for_merging_methylation_data_using_paste.py
```

>Split the methylation metrix into CG, CHG, CHH-type methylation files

```
python 17_split_methylation_to_CG_CHH_CHG.py inputFile
```

#### 2.3.2.2 For methylation proportion matrix

Related scripts can be found in the folder Data_preprocessing/03_methylomic_matrix/3_single_site_methylation_Prop

>Calculate the methylation proportion for each C site, and save the methylated sites. Write the slurm jobs to parse the raw single-site methylation files separately. your_work_dir is the directory containing all the small files

```
python 01_write_slurm_jobs_for_methylation_proportion_and_save_methylated_site.py your_work_dir
```
	
>>Note: for individual files, you can run awk command lines, and the input file is the downloaded *.tsv file 

```
awk \'{print "Chr"$1"_"$2"_"$4"_"$3"\\t"$5"/"$6"\\t"$5}\' < inputFile > inputFile_proportion
```
```	
awk \'$3>0\' < inputFile_proportion > inputFile_meted.txt
```	

>Make dictionary for single-site methylation proportion, and the input file is *_meted.txt

```	
python 02_make_dic_for_methylation_proportion.py input_file 
```	

>*Repeat what have been done for methylation P/A analysis*. Or if you want to keep the same methylation sites for methy_prop and methy_P/A, you can just fill the corresponding values in methy_P/A matrix using methylation proportion before imputation. The inputFile is the *_50per file.

```
python 03_fill_single_base_methylation_as_proportion.py inputFile
```

>Impute the methylation proportion matrix

```
python 04_knn_imputation_for_proportion_training_fit_on_test.py inputFile
```

### 2.3.2.3 For binning of methylation presence/absence or proportion values

Related scripts can be found in the folder Data_preprocessing/04_methylation_binning

### 2.3.2.4 For clustering of binned methylation presence/absence or proportion values

Related scripts can be found in the folder Data_preprocessing/05_methylation_profile_clustering

# **3. get the Kinship and correlation matrix for omics data**

Related scripts can be found in the folder Data_preprocessing/06_correlation_matrix_for_omics_data

## 3.1 Get the kinship

>Convert the SNP matrix to hmp format

```
python 01_genotype2hapmap.py SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05.csv SNP_383_accessions.hmp.txt
```

>To get the kinship matrix, please follow the illustration of the software tassel, which can be found in https://www.maizegenetics.net/tassel

```
tassel/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile SNP_383_accessions.hmp.txt -outputFile SNP_383_accessions_order.hmp.txt -fileType Hapmap -Xmx100g -Xms2g
```
```
tassel/tassel-5-standalone/run_pipeline.pl -importGuess SNP_383_accessions_order.hmp.txt -KinshipPlugin -method Centered_IBS -endPlugin -export SNP_383_accessions_kinship.txt -exportType SqrMatrix -Xmx100g -Xms2g
```

## 3.2 Get the correlation matrix for other multi-omics data

```  
Rscript 02_get_correlation_matrix_for_omics_data.r
```

## 3.3 Remove the potential confounding effects of K on mCor. You may want to replace input files in this script

```
Rscript 03_removing_confounding_effects_of_K_from_mCor.r
```

# **4. Genomic prediction using machine learning algorithms**

Example data can be found in the folder [/Example_data_for_model_building](https://github.com/ShiuLab/Manuscript_Code/tree/master/2023_Ath_multi-omics/Example_data_for_model_building).

For the RandomForest regression model building, we used the script from Azodi et al., 2020. Plant Cell paper. Here is the link to the script [ML_regression.py](https://github.com/ShiuLab/Manuscript_Code/blob/master/2019_expression_GP/scripts/ML_regression.py). Example run as below:

```
python ML_regression.py -df SNP_383_accessions_kinship.csv_training.csv -df2 Phenotype_value_383_common_accessions_2017_Grimm.csv -y_name Length -sep ',' -alg RF -gs T -cv_set CVFs.csv -save RF_kinship_Length -test Test.txt -n_jobs 4 -n 10 -cv_num 5
```

For the rrBLUP model, we used the script from our another project. Here is the link to the script [13_rrBLUP_training_test_split_fread_predict_values.r](https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_GP_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r). Example run as below:

```
Rscript 13_rrBLUP_training_test_split_fread_predict_values.r SNP_383_accessions_kinship.csv_training.csv Phenotype_value_383_common_accessions_2017_Grimm.csv all Length Test.txt 5 10 CVFs.csv Kinship
```

# **5. SHAP values and feature interaction values for RF models**

The scrips and example data can be found in the folder [/SHAP](https://github.com/ShiuLab/Manuscript_Code/tree/master/2023_Ath_multi-omics/SHAP)

>Note: the feature interaction value calculation is pretty slow for datasets with a large number of features. In our study, we only took the benchmark flowering time gene-related features to calculate the feature interactions.

Example run for the script SHAP_training_only_saving_interaction_figures_for_given_feature_list.py. Example input data can be found in the folder [/SHAP/Example_data](https://github.com/ShiuLab/Manuscript_Code/tree/master/2023_Ath_multi-omics/SHAP/Example_data)

```
python SHAP_training_only_saving_interaction_figures_for_given_feature_list.py -df Matrix_top_GTM_features_for_426_flowering_time_genes.csv -df2 Phenotype_value_383_common_accessions_2017_Grimm.csv -sep "," -y_name FT10_mean -test Test.txt -save SHAP_426_benchmark_top_GTM -model RF_top_GTM_426_benchmark_genes_model.pkl -n_jobs 16 -top 20 -interaction y -interaction_score y -feature_list Feature_list_selected_for_426_flowering_genes.txt
```

To summarize and order the feature interaction values, please run the script Sum_SHAP_interaction.py. The path contains all the interaction value files output from the above script. Example interaction files can be found in the folder [/SHAP/Example_data/Interaction_files](https://github.com/ShiuLab/Manuscript_Code/tree/master/2023_Ath_multi-omics/SHAP/Example_data/Interaction_files), where only interaction values among 30 features in two individuals were kept.

```
python Sum_SHAP_interaction.py -path ./ -save Summarized_and_ordered_interactions
```