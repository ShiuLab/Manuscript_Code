library(data.table)
args = commandArgs(trailingOnly=TRUE)
geno_file <- args[1]
X <- as.matrix(fread(geno_file),rownames=1)
print('Done with the reading')

X2=scale(X) # Centers (subtract the column means) and Scales (dividing the centered columns by their stdev)
G=tcrossprod(X2) # Take the cross product X transpose
G=G/mean(diag(G))
EVD=eigen(G)
rownames(EVD$vectors)=rownames(G)
save(EVD,file='EVD.RData')
write.csv(EVD$vectors,"PCA_matrix.csv",row.names=T,quote=F)
write.csv(EVD$vectors[,1:5],"PCA5_geno.csv",row.names=T,quote=F)

