library(rrBLUP)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

X_file <- args[1] # geno matrix
Y_file <- args[2] # pheno matrix
test_file <- args[3] # individuals in test set


geno <- fread(X_file)
test <- scan(test_file, what='character')
training <- geno[geno$ID %in% setdiff(geno$ID,test),]
write.csv(training,'geno_training.csv',quote=F,row.names=F)

Y <- read.csv(Y_file, row.names=1) 
Y <- Y[!rownames(Y) %in% test,]
write.csv(Y,'pheno_training.csv',row.names=T,quote=F)
