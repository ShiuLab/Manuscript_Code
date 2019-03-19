### Plots for FT benchmark genes

library('data.table')
setwd('/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/')

g <- fread('01_Data/geno.csv', header=T, sep=',')
g <- as.data.frame(g)

t <- fread('01_Data/transcripto.csv', header=T, sep=',')
t <- as.data.frame(t)
row.names(t) <- t$ID
t$ID <- NULL

p <- read.table('01_Data/pheno.csv', header=T, sep=',')
p <- p[match(rownames(g), p$ID),]

gp <- setDT(g)[p, on='ID']
row.names(g) <- g$ID
g$ID <- NULL

#### T-test for all markers by flowering time for benchmark genes
print('Running t-tests on markers & FT')
g_ttest <- data.frame(snp = character(), pval = character(), stringsAsFactors=FALSE)

for (snp in colnames(g)){
  if(snp != 'ID'){
    tmp_mod <- as.formula(paste('FT', snp, sep='~'))
    try(g_ttest <- rbind(g_ttest, list(snp = snp, pval = t.test(tmp_mod, gp)$p.value, stringsAsFactors=F)))
  }
}
write.table(g_ttest, 'G_P_lm_results.txt', sep='\t', quote=F, row.names = F, col.names = T)

#### Plot expression level by flowering time for benchmark genes

print('Running linear models on transcripts & FT')

t_r2 <- data.frame(transcript = character(), r2 = character(), pval = character(), stringsAsFactors=FALSE)
t <- cbind(FT = p$FT, t)

for(trans in colnames(t)){
  t2 <- cbind(FT = p$FT, t)
  tmp_mod <- as.formula(paste('t2$FT ~ t2$', trans, sep=''))
  model_temp <- summary(lm(as.formula(tmp_mod)))
  r2 <- model_temp$r.squared
  pval <- model_temp$coefficients[,4][2]
  t_r2 <- rbind(t_r2, list(transcript = trans, r2 = r2, pval = pval), stringsAsFactors=F)
}
t_r2 <- subset(t_r2, transcript != 'FT')
write.table(t_r2, 'T_P_lm_results.txt', sep='\t', quote=F, row.names = F, col.names = T)

print('Done!')