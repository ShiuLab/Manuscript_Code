# calculate the TPM
library(scater)
library('stringr')
library(devtools)
library(Biobase)
library(preprocessCore)

# the transcriptomic data was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz 
Exp <- read.table('GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv',head=T,sep='\t',stringsAsFactors=F,row.names=1)
dat <- read.table('TAIR10_longest_mRNA_length_including_ncRNA.txt',head=T,sep='\t',stringsAsFactors=F,row.names=1)
rownames(dat) <- substr(rownames(dat),1,9)
Expression <- merge(Exp,dat, by="row.names")
rownames(Expression) <- Expression$Row.names
Expression <- Expression[,2:ncol(Expression)]
Len <- Expression$Length
tpm_table <- calculateTPM(as.matrix(Expression[,1:(ncol(Expression)-1)]), Len)
write.table(tpm_table,'GSE80744_TPM.txt',row.names=T,sep='\t',quote=F)
