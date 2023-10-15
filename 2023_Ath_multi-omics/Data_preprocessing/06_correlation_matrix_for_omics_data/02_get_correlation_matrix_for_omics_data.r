#####################################
# remove zero variants
library(caret)
# get the eCor, using transformed TPM
setwd('D:\\Ath_GS\\Model_for_2017_Grimm')
dat <- read.csv('GSE80744_TPM_383_accessions_loge_plus1.csv',head=T,row.names=1)
datnew <- nearZeroVar(dat,uniqueCut = 5, freqCut = 95/5, saveMetrics = TRUE)
dat <- t(dat)
eCor <- cor(dat)
write.csv(eCor,'eCor_383_accession.csv',row.names=T,quote=F)
eCor_nor <- (eCor-min(eCor))/(max(eCor)-min(eCor))
write.csv(eCor_nor,'eCor_normalized_383_accession.csv',row.names=T,quote=F)

# get the eCor, using original TPM
setwd('D:\\Ath_GS\\Model_for_2017_Grimm')
dat <- read.csv('GSE80744_TPM_383_accessions.csv',head=T,row.names=1)
datnew <- nearZeroVar(dat,uniqueCut = 5, freqCut = 95/5, saveMetrics = TRUE)
dat <- t(dat)
eCor <- cor(dat)
write.csv(eCor,'eCor_non_transform_383_accession.csv',row.names=T,quote=F)
eCor_nor <- (eCor-min(eCor))/(max(eCor)-min(eCor))
write.csv(eCor_nor,'eCor_non_transform_normalized_383_accession.csv',row.names=T,quote=F)

# get the mCor
dat <- read.csv('Araport11_GB_methylation_383_accessions_knn_imputed.csv',head=T,row.names=1)
dat <- t(dat)
mCor <- cor(dat)
write.csv(mCor,'mCor_383_accession.csv',row.names=T,quote=F)
mCor_nor <- (mCor-min(mCor))/(max(mCor)-min(mCor))
write.csv(mCor_nor,'mCor_normalized_383_accession.csv',row.names=T,quote=F)

# get the pCor
dat <- read.csv('Phenotype_value_383_common_accessions_2017_Grimm.csv',head=T,row.names=1)
dat <- t(dat)
pCor <- cor(dat[1:5,])
write.csv(pCor,'pCor_383_accession.csv',row.names=T,quote=F)
pCor_nor <- (pCor-min(pCor))/(max(pCor)-min(pCor))
write.csv(pCor_nor,'pCor_normalized_383_accession.csv',row.names=T,quote=F)
