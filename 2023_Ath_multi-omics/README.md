# multi-omics
Codes for multi-omics integration project 

## Datasets
## 1. SNP matrix
The SNP matrix was download from https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX.tar.gz

### 2. transcriptomic data
The transcriptomic data was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz
which are read count files for 727 accessions

### 3. methylomic data
The gene body methylation data was download from http://signal.salk.edu/1001.php (which is currently not accessible), which contains 1107 methylomes, with the ID name.
The tsv files for individual accessions were download from NCBI with GEO accession ID: GSE43857

### 4. phenotypic data
Flowering time at 10℃ (1163 accessions) and 16℃ (1123 accessions) were downloaded from https://arapheno.1001genomes.org/study/12/
Other phenotypic data were download from https://arapheno.1001genomes.org/study/38/

### 5. benchmark flowering time genes
Benchmark flowering time genes were downloaded from the FLOR-ID database http://www.phytosystems.ulg.ac.be/florid/