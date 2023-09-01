#### get rid of the confounding effect of kinship on mCor, calculated on training instances, 383 accessions
kinship <- as.matrix(read.csv('SNP_383_accessions_kinship.csv',head=T,row.names=1))
mCor <-as.matrix(read.csv('mCor_normalized_383_accession.csv',head=T,row.names=1))
Test <- read.table('Test.txt',head=F,sep='\t')
kinship <- kinship[,!colnames(kinship) %in% paste('X',Test[,1],sep='')]
mCor <- mCor[rownames(kinship),]
mCor <- mCor[,colnames(kinship)]
Res <- mCor
for(i in 1:nrow(mCor)){
	model <- lm(mCor[i,] ~ kinship[i,])
	Res[i,] <- model$residuals
	}
write.csv(Res,'mCor_normalized_383_accession_residual_from_kinship_training.csv',row.names=T,quote=F)
