#install.packages('MatrixEQTL')
#source("Matrix_eQTL_R/Matrix_eQTL_engine.r");

args = commandArgs(trailingOnly=TRUE)

library(MatrixEQTL)
library(data.table)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');

## Input
SNP_file_name = args[1]
expression_file_name = args[2]
save = args[3]
transcript_loc_name = "/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/00_RawData/v3_v4_xref.txt"

# Output
output_file_name_cis = paste("eqtl_cis_", save, ".csv", sep="")
output_file_name_tra = paste("eqtl_trans_", save, ".csv", sep="")

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6 # Used by Zan et al - Genetic regulation of transcriptional variation in natural A. thaliana accessions. 2013

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Format SNP and gene position dataframes.
snpspos <- as.data.frame(rownames(snps))
snpspos$chr <- as.character(lapply(strsplit(as.character(snpspos[,'rownames(snps)']), split="_"),head, n=1))
snpspos$pos <- as.character(lapply(strsplit(as.character(snpspos[,'rownames(snps)']), split="_"),tail, n=1))
snpspos$chr <- gsub('Chr','', snpspos$chr)
names(snpspos) <- c('snpid', 'chr', 'pos')
snpspos$chr <- as.numeric(snpspos$chr)
snpspos$pos <- as.numeric(snpspos$pos)

genepos = read.table(transcript_loc_name, sep='\t',header = TRUE, stringsAsFactors = FALSE);
genepos = genepos[,c('v3_gene_model','v4_chr', 'v4_start', 'v4_end')]
genepos$left <- apply(genepos[,c('v4_start','v4_end')], 1, min)
genepos$right <- apply(genepos[,c('v4_start','v4_end')], 1, max)
genepos$v4_chr <- gsub('Chr','', genepos$v4_chr)
genepos <- genepos[,c('v3_gene_model', 'v4_chr', 'left','right')]
names(genepos) <- c('geneid', 'chr', 'left', 'right')
genepos$chr <- as.numeric(genepos$chr)
genepos$left <- as.numeric(genepos$left)
genepos$right <- as.numeric(genepos$right)


## Run the analysis (took ~ 1 hour...)
me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  output_file_name      = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis  = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

print('DONE!')