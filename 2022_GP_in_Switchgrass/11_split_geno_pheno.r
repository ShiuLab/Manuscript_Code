library(rrBLUP)
set.seed(42)
args = commandArgs(trailingOnly=TRUE)
X_file <- args[1] # geno matrix
Y_file <- args[2] # pheno matrix
test_file <- args[3] # individuals in test set

Y <- read.csv(Y_file, row.names=1) 
X <- read.csv(X_file, row.names=1) 
Test <- read.table(test_file, head=F,sep='\t',stringsAsFactors=F) 

X_test <- X[Test[,1],]
Y_test <- Y[Test[,1],]
X <- X[!rownames(X) %in% Test[,1],, drop=False]
Y <- Y[!rownames(Y) %in% Test[,1],, drop=False]
write.csv(X,'geno_training.csv',row.names=T,quote=F)
write.csv(Y,'pheno_training.csv',row.names=T,quote=F)
write.csv(X_test,'geno_test.csv',row.names=T,quote=F)
write.csv(Y_test,'pheno_test.csv',row.names=T,quote=F)
