# multi-omics
Codes for multi-omics integration project 

## 1. Datasets

## 1.1 SNP matrix
The SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz

### 1.2 transcriptomic data
The transcriptomic data was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz
which are read count files for 727 accessions

### 1.3 methylomic data
The gene body methylation data was download from http://signal.salk.edu/1001.php (which is currently not accessible), which contains 1107 methylomes, with the ID name.
The tsv files for individual accessions were download from NCBI with GEO accession ID: GSE43857

### 1.4 phenotypic data
Flowering time at 10℃ (1163 accessions) and 16℃ (1123 accessions) were downloaded from https://arapheno.1001genomes.org/study/12/
Other phenotypic data were download from https://arapheno.1001genomes.org/study/38/

### 1.5 benchmark flowering time genes
Benchmark flowering time genes were downloaded from the FLOR-ID database http://www.phytosystems.ulg.ac.be/florid/

## 2. Data_preprocessing

### 2.1 SNP matrix
The scripts for SNP matrix can be found in the folder Data_preprocessing\SNP_matrix

### 2.2 transcriptomic data

  * To calculate the TPM, we need to get the transcript length for genes at first. Get the GFF3 file ready before run this script. Here we used the TAIR10 GFF3 file 
  
	`python 01_get_transcript_length_including_ncRNA.py`

### 2.3 methylomic data
