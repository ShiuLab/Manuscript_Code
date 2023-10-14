# multi-omics integration project
Codes for our manuscript "Prediction of plant complex traits via integration of multi-omics data"

## **1. Datasets**

### 1.1 SNP matrix

The SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz

### 1.2 transcriptomic data
The transcriptomic data was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz
which are read count files for 727 accessions

### 1.3 methylomic data
The gene body methylation data was download from http://signal-genet.salk.edu/1001.php , which contains 1107 methylomes, with the ID name.
The tsv files for individual accessions were download from NCBI with GEO accession ID: GSE43857

### 1.4 phenotypic data
Flowering time at 10℃ (1163 accessions) and 16℃ (1123 accessions) were downloaded from https://arapheno.1001genomes.org/study/12/
Other phenotypic data were download from https://arapheno.1001genomes.org/study/38/

### 1.5 benchmark flowering time genes

Benchmark flowering time genes were downloaded from the FLOR-ID database http://www.phytosystems.ulg.ac.be/florid/

## **2. Data_preprocessing**

### **2.1 SNP matrix**

Related scripts for SNP matrix can be found in the folder Data_preprocessing\01_SNP_matrix

  *  Get the SNP matrix

```  
python 01_extract_the_h5py_binary_SNP_matrix.py
```
	
  *  Save the SNP matrix for 383 accession, convert common allele to 1, and not common allele to -1, get rid of the rare variants (MAF < 5%)	
  
```
python 02_get_snp_matrix_for_383_accessions.py
```
	
### **2.2 transcriptomic data**

Related scripts can be found in the folder Data_preprocessing\02_transcriptomic data 
  
  *  To calculate the TPM, we need to get the transcript length for genes at first. Get the GFF3 file ready before run this script. Here we used the TAIR10 GFF3 file 

```  
python 01_get_transcript_length_including_ncRNA.py
```

  *  Calculate TPM using the function “calculateTPM” from the R package “scater”
  *  Put the downloaded file GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv and the transcript length file TAIR10_longest_mRNA_length_including_ncRNA.txt in your work directory
  
```
Rscript 02_calculate_TPM.r
```

  *  Log the TPM and get the logged TPM matrix for 383 accessions
  
```
python 03_log_TPM_and_get_transcriptomic_matrix_for_383_accessions.py
```	
	
### **2.3 methylomic data**

#### **2.3.1 For gene-body methylation**

Related scripts can be found in the folder Data_preprocessing\03_methylomic_matrix

  *  Parse the downloaded gene-body methylation matrix, and save the corresponding matrix for 383 accessions

```
python 01_get_methylation_for_383_accessions.py
```

  *  Impute missing data
 
```
python 02_knn_imputation.py
```

#### **2.3.2 For single-site methylation-based formats**
	
  *  Download individual methylation data for 383 accessions 

```
python 03_download_individual_methylation_data_for_383_accessions.py
```
	
##### 2.3.2.1 For presence/absence of methylation
	
  *  For presence/absence of methylation, save the methylated sites 

```
python 05_write_slurm_jobs_for_04.py inputFile
```
	
  *  write slurm jobs for all downloaded methylation files
  
```
python 06_write_slurm_jobs_for_03.py yourWorkPath
```

  *  combine all files with single-site presence/absence of methylation information

```
python 07_combine_all_methylation_files.py
```

```
python 08_make_dic_for_methylation_PA_proportion.py input_file
```


#### 2.3.2.2 For methylation proportion

  *  Calculate the methylation proportion for each C site, and save the methylated sites. Write the slurm jobs to parse the raw single-site methylation files separately.

```
python 05_write_slurm_jobs_for_methylation_proportion_and_save_methylated_site.py
```
	
>>>>Note: for individual files, you can run awk command lines

```
>>>>awk \'{print "Chr"$1"_"$2"_"$4"_"$3"\\t"$5"/"$6"\\t"$5}\' < inputFile > inputFile_proportion
```
```	
>>>>awk \'$3>0\' < inputFile_proportion > inputFile_meted.txt
```	


## **3. get the Kinship and correlation matrix for omics data**

Related scripts can be found in the folder Data_preprocessing\04_correlation_matrix_for_omics_data

### 3.1 get the kinship

  *  Convert the SNP matrix to hmp format

```
python 01_genotype2hapmap.py SNP_binary_matrix_383_accessions_drop_all_zero_MAF_larger_than_0.05.csv SNP_383_accessions.hmp.txt
```

  *  To get the kinship matrix, please follow the illustration of the software tassel, which can be found in https://www.maizegenetics.net/tassel

```
tassel/tassel-5-standalone/run_pipeline.pl -SortGenotypeFilePlugin -inputFile SNP_383_accessions.hmp.txt -outputFile SNP_383_accessions_order.hmp.txt -fileType Hapmap -Xmx100g -Xms2g
```
```
tassel/tassel-5-standalone/run_pipeline.pl -importGuess SNP_383_accessions_order.hmp.txt -KinshipPlugin -method Centered_IBS -endPlugin -export SNP_383_accessions_kinship.txt -exportType SqrMatrix -Xmx100g -Xms2g
```

  *  Get correlation matrix for other multi-omics data

```  
Rscript 02_get_correlation_matrix_for_omics_data.r
```

  *  Remove the potential confounding effects of K on mCor. You may want to replace input files in this script

```
Rscript 04_removing_confounding_effects_of_K_from_mCor.r
```

## **4. Genomic prediction using machine learning algorithms**



## **5. SHAP values and feature interaction values**

  * Note: the feature interaction value calculation is pretty small for datasets with a large number of features. In our study, we only took the benchmark flowering time gene-related features to calculate feature interactions.