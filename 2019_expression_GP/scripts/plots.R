library(ggplot2)
library(reshape2)
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')


# Model accuracy using different data types

# Get summary stats for rrBLUP and BL 
d <- read.table('accuracy_withEnsemble.txt', sep='\t', header=T)
d_mean <- aggregate(d, by=list(Alg=d$model, X=d$x_file, Y=d$y), mean)
d_sd <- aggregate(d, by=list(Alg=d$model, X=d$x_file, Y=d$y), sd)
d_mean$PCC_sd <- d_sd$accuracy_PCC
d2 <- d_mean[c('Alg', 'X', 'Y','accuracy_PCC', 'PCC_sd')]
names(d2) <- c('Alg', 'X', 'Y','PCC', 'PCC_sd')
d2$X <- gsub('eCorNorm','transcripto_cor_norm', d2$X)


# Format ML results to match
rf <- read.table('05_ML/JP/results.txt', sep='\t', header=T)
#rf <- read.table('RESULTS_reg.txt', sep='\t', header=T)
rf$X <- as.character(lapply(strsplit(as.character(rf$ID), split="/"),tail, n=1))
d_rf <- rf[c('Alg', 'X', 'Y','PCC', 'PCC_sd')]
d_rf$X <- gsub('[FT/YLD/HT/_RF]',"", d_rf$X)
d_rf$X <- gsub('genoranscripto','gt', d_rf$X)
d_rf$X <- gsub('transcriptokin','kt', d_rf$X)
#d_rf$X <- gsub('','transcripto_cor_norm', d_rf$X)

# Merge results & find top performer
results <- rbind(d2, rf)
results$PCC <- as.numeric(results$PCC)
dcast(results, Alg + X ~ PCC)
max <- aggregate(results$PCC, by=list(results$Y), max)
names(max) <- c('Y','max')

x_files <- c('PC5', 'kin','geno','transcripto_cor_norm', 'transcripto','kt', 'gt')
reo_x_files <- function(x) { factor(x, levels = x_files)}
reo_Alg <- function(x) { factor(x, levels = c('ensem','RF','BL', 'rrBLUP'))}

results <- subset(results, results$X != 'PC5')
results <- subset(results, results$X != 'transcripto_shuff')
results <- subset(results, results$X != 'tshuff')
results <- subset(results, results$X != 'transcripto_cor')

library(dplyr)
normalit<-function(m){
  (m - min(m))/(max(m)-min(m))
}
results_norm <- results %>%
  group_by(Y) %>%
  mutate(PCC_norm = normalit(PCC))
unique(results_norm$X)
ggplot(results_norm, aes(x=reo_x_files(X), y=reo_Alg(Alg))) + geom_tile(aes(fill=PCC_norm), colour='white') +
  facet_grid(Y~.) + labs(x="Input Data", y="Algorithm") + theme_bw(10) +
  geom_text(aes(label = round(PCC, 2))) +
  scale_fill_gradient2(low='steelblue', mid = "white", high = "tomato",midpoint = 0.5, limits = c(0,1),na.value = "grey50")
#ggsave("00_Figures/performance_heatmap_180727.pdf", width = 5.2, height = 4)

ggplot(results_norm, aes(reo_x_files(X), y=PCC)) + geom_violin() + geom_boxplot(width=.05) + theme_bw(10) 
#ggsave("00_Figures/performance_violin_X_1800727.pdf", width = 5.2, height = 2)

results_norm$PopStru <- results_norm$Y
results_norm$PopStru <- gsub('FT',0.294, results_norm$PopStru)
results_norm$PopStru <- gsub('HT',0.318, results_norm$PopStru)
results_norm$PopStru <- gsub('YLD',0.425, results_norm$PopStru)
results_norm$PopStru <- as.numeric(results_norm$PopStru)
ggplot(results_norm, aes(x=reo_Alg(Alg), y=PCC)) + geom_violin() + 
  geom_boxplot(width=1)+ theme_bw(10) + facet_grid(Y~.) + coord_flip() +
  geom_hline(data = results_norm, aes(yintercept=PopStru, colour="red", linetype="dashed"))
#ggsave("00_Figures/performance_violin_Alg_1800727.pdf", width = 2, height = 4)



#### Compare PCCs using scatterplots

library(ggplot2)
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')

p <- read.table('01_Data/pheno.csv', header=T, sep=',', row.names = 1)
trait <- 'FT'
r <- as.data.frame(p[[trait]])
row.names(r) <- row.names(p)
names(r) <- c('true')

rr.g <- read.table('03_rrBLUP/rrBLUP_geno_FT_yhat.csv', sep=',', header=T)



##### Plot change in ro when LD changed in transcript-SNP comparison
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/06_CompareCoefImp')
ro <- read.table('ro_values.csv', sep=',', header=TRUE)
ro.m <- melt(ro, id='LD')
ro.m$window <- ro.m$LD*2
ggplot(ro.m, aes(x=log10(window), y=value, group=variable, color=variable)) +
  labs(x='log(window size)', 'ro') + facet_grid(variable~., scales='free') +
  geom_point() + geom_line() + theme_bw(8)

ggsave("~/Desktop/181002_ro_by_LD_log10.pdf", width = 3, height = 3,useDingbats=FALSE )



##### Plot importance of features in G+T model
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')
d <- fread('03_rrBLUP/rrBLUP_gt_FT_coef.csv')
dm <- colMeans(d[,5:length(d)])
dm <- as.data.frame(dm)
dm$ID <- row.names(dm)

dm$type <- row.names(dm)
dm$type <- replace(dm$type, grep("Chr", dm$type), 'G')
dm$type <- replace(dm$type, grep("B73V4", dm$type), 'G')
dm$type <- replace(dm$type, grep("GRMZ", dm$type), 'T')
dm$type <- replace(dm$type, grep("AC", dm$type), 'T')
dm$type <- replace(dm$type, grep("AF", dm$type), 'T')
dm$type <- replace(dm$type, grep("AY", dm$type), 'T')
dm$type <- replace(dm$type, grep("EF517", dm$type), 'T')

summary(lm(dm~type, data=dm))
quantile(dm$dm)
dm2 <- dm[with(dm,order(-dm)),][1:1000,]
summary(as.factor(dm2$type))
dm3 <- dm[with(dm,order(-dm)),][1001:dim(dm)[1],]
summary(as.factor(dm3$type))

# Enrichement for T in top 1000
fisher.test(matrix(c(27,973,31210,331204), nrow=2))

dm2 <- dm[with(dm,order(-dm)),][1:20,]
summary(as.factor(dm2$type))
dm3 <- dm[with(dm,order(-dm)),][21:dim(dm)[1],]
summary(as.factor(dm3$type))
# Enrichment for T in top 20?
fisher.test(matrix(c(4,16,31233,332161), nrow=2))

ggplot(dm2, aes(x=reorder(ID, -dm), y=dm, fill=type)) + 
  geom_bar(stat='identity', position= "dodge") + theme_bw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major.x = element_blank()) +
  coord_cartesian(ylim=c(0.02,0.0367))
#ggsave("00_Figures/gt_imp_top20_1800716.pdf", width = 3, height = 3)
#ggsave("00_Figures/gt_imp_top1000_1800716.pdf", width = 3, height = 3)


### Correlation between importance scores from G, T, and G+T
g <- fread('03_rrBLUP/rrBLUP_geno_FT_coef.csv')
gm <- colMeans(g[,5:length(g)])
gm <- as.data.frame(gm)
gm$ID <- row.names(gm)

gm2 <- merge(gm, dm, by='ID')
cor.test(gm2$gm, gm2$dm)

t <- fread('03_rrBLUP/rrBLUP_transcripto_FT_coef.csv')
tm <- colMeans(t[,5:length(t)])
tm <- as.data.frame(tm)
tm$ID <- row.names(tm)
tm2 <- merge(tm, dm, by='ID')
cor.test(tm2$tm, tm2$dm)
cor(tm2$tm, tm2$dm, method='spearman')

ggplot(tm2, aes(x=tm, y=dm)) + geom_point(size=0.2, alpha=0.2) + 
  theme_bw(8) + geom_smooth(method='lm',formula=y~x, se=FALSE)
#ggsave("00_Figures/cor_imp_T_GT_FTrrB_1800716.pdf", width = 3, height = 3, useDingbats=FALSE)

ggplot(gm2, aes(x=gm, y=dm)) + geom_point(size=0.2, alpha=0.2) + 
  theme_bw(8) + geom_smooth(method='lm',formula=y~x, se=FALSE)
#ggsave("00_Figures/cor_imp_G_GT_FTrrB_1800716.pdf", width = 3, height = 3, useDingbats=FALSE)





#### Plot PCC from using top 200 trans or SNPs from RF and rrBLUP for FT
# Order rrBLUP_SNP - rrBLUP_trans - RF_SNP - RF_trans

rrblup_top200 <- read.table('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/03_rrBLUP/01_TopTrans_vs_SNPs/accuracy.txt', sep=' ')
rrblup_top200_mean <- aggregate(rrblup_top200, by=list(X=rrblup_top200$V2), mean)
rrblup_top200_sd <- aggregate(rrblup_top200, by=list(X=rrblup_top200$V2), sd)

pcc = c(0.8661414, 0.7419128, 0.707616855823,0.704112568981)
sd = c(0.004561519, 0.010199201, 0.00869422280128, 0.0101994534297)
alg = c('rrBLUP','rrBLUP','RF','RF')
type = c('G','T','G','T')
top20res <- data.frame(alg,type,pcc,sd)
ggplot(top20res, aes(x=type, y=pcc, fill=type)) + coord_flip() + 
  geom_bar(stat='identity', position= "dodge") + facet_grid(.~alg) +
  geom_hline(yintercept=0.5, linetype="dashed", size=0.3) + coord_cartesian(ylim = c(0.5, 0.9)) +
  geom_errorbar(aes(ymin=pcc-sd, ymax=pcc+sd), size = 0.3, width = 0.1, position= position_dodge(0.9))




##### Feature Selection on FT #####

setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/08_FeatureSelection/')

# Get summary stats for rrBLUP and BL and RF
d <- read.table('accuracy.txt', sep='\t', header=T)
names(d) <- c('Alg','type','Y','feat_method','FeatureNum','PCC_ho','PCC_ho_sd')
dx <- d[c('Alg','Y','type','FeatureNum','PCC_ho','PCC_ho_sd')]

rf <- read.table('RESULTS_reg.txt', sep='\t', header=T)
rf$type <- substr(rf$ID, 0,1)
rfx <- rf[c('Alg','Y','type','FeatureNum','PCC_ho','PCC_ho_sd')]


# Merge results & find top performer
results <- rbind(dx, rfx)
results[results$FeatureNum > 5001,]$FeatureNum <- 7000
results$PCC <- as.numeric(results$PCC)
dcast(results, Alg + X ~ PCC)


ggplot(results, aes(x=log10(FeatureNum), y=PCC_ho, color=Alg, linetype=type)) + 
  geom_point() + geom_line() + theme_bw(10) 
#ggsave("00_Figures/performance_heatmap_180716.pdf", width = 5.2, height = 4)

ggplot(results_norm, aes(reo_x_files(X), y=PCC)) + geom_violin() + geom_boxplot(width=.05) + theme_bw(10) 
#ggsave("00_Figures/performance_violin_X_1800705.pdf", width = 5.2, height = 2)


######### Compare expression and allele patterns for FT benchmark genes
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/01_Data/')
ft_g <- fread('geno.csv', sep=',',header=T)
ft_t <- fread('transcripto.csv', sep=',')


