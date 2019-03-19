library(ggplot2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(reshape2)


# Accuracy summary
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')
k <- read.csv('01_Data/kin.csv', header=T, row.names=1)
p <- read.csv('01_Data/pheno.csv', header=T, row.names=1)
t <- read.csv('01_Data/transcripto.csv', header=T, row.names=1)

# Make sure Y is in the same order as X:
p <- p[match(rownames(t), rownames(p)),]
k <- k[match(rownames(p), rownames(k)),]

## Plot heatmaps of kinship matrix and expression correlation with phenotypes
# 5/25/18
tc <- cor(t(t))
dim(tc)

colnames(k) <- row.names(k)
colnames(tc) <- row.names(tc)
tc_med <- median(tc)
summary(as.matrix(as.vector(tc)))

#Normaliz & save normalied version for BL runs...
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
colMax <- function(data) sapply(data, max, na.rm = TRUE)
diag(tc) <- NA
diag(tc) <- colMax(as.data.frame(tc))
t_norm <- range01(tc)
diag(t_norm) <- 1
write.csv(t_norm, '01_Data/transcripto_cor_norm.csv', row.names=T, quote=F)



ht1 <- Heatmap(tc, name='Exp Correlation',  show_row_names = F, show_column_names = F) #col=colorRamp2(c(0.89,tc_med,1), c('#5ab4ac','white','#d8b365')),
ht1_order <- column_order(ht1)
htp1 <- Heatmap(p$HT, name='HT', col=colorRamp2(c(80,240), c('white', '#e6550d')), width=unit(5, 'mm')) 
htp2 <- Heatmap(p$FT, name='FT', col=colorRamp2(c(500,1300), c('white', '#31a354')), width=unit(5, 'mm')) 
htp3 <- Heatmap(p$YLD, name='YLD', col=colorRamp2(c(0,4.1), c('white', '#756bb1')), width=unit(5, 'mm'), clust) 
ht2 <- Heatmap(as.matrix(k), name='Kinship', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order) 
ht1 + htp1 + htp2 + htp3 + ht2


# Remake plot with trait distance matrix instead of raw values
pc <- cor(t(p))
dim(pc)
colnames(k) <- row.names(k)
colnames(tc) <- row.names(tc)
colnames(pc) <- row.names(pc)
tc_med <- median(tc)

ht1 <- Heatmap(tc, name='Exp Correlation',  show_row_names = F, show_column_names = F) #col=colorRamp2(c(0.89,tc_med,1), c('#5ab4ac','white','#d8b365')),
ht1_order <- column_order(ht1)
htp <- Heatmap(pc, name='Phenotype Distance', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order) 
ht2 <- Heatmap(as.matrix(k), name='Kinship', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order) 
ht1 + htp + ht2

# Remake using euclidian distance instead of correlation
#install.packages('distances')
library(distances)
p_dis <- as.matrix(distances(p))
diag(p_dis) <- min(p_dis)



ht1 <- Heatmap(tc, name='Exp Correlation',  show_row_names = F, show_column_names = F) #col=colorRamp2(c(0.89,tc_med,1), c('#5ab4ac','white','#d8b365')),
ht1_order <- column_order(ht1)
htp <- Heatmap(p_dis, name='Phenotype Distance', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order) 
ht2 <- Heatmap(as.matrix(k), name='Kinship', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order) 
ht1 + htp + ht2

ht1 <- Heatmap(t_dis, name='Exp Correlation',  show_row_names = F, show_column_names = F) #col=colorRamp2(c(0.89,tc_med,1), c('#5ab4ac','white','#d8b365')),
ht1_order <- column_order(ht1)
htp1 <- Heatmap(p$HT, name='HT', col=colorRamp2(c(80,240), c('white', '#e6550d')), width=unit(5, 'mm')) 
htp2 <- Heatmap(p$FT, name='FT', col=colorRamp2(c(500,1300), c('white', '#31a354')), width=unit(5, 'mm')) 
htp3 <- Heatmap(p$YLD, name='YLD', col=colorRamp2(c(0,4.1), c('white', '#756bb1')), width=unit(5, 'mm')) 
ht2 <- Heatmap(as.matrix(k), name='Kinship', cluster_columns=FALSE, show_row_names = F, 
               show_column_names = F, column_order = ht1_order, 
               col=colorRamp2(c(min(k),summary(as.vector(as.matrix(k)))['Median'], max(k)), c('blue', 'white', 'red'))) 
ht1 + htp1 + htp2 + htp3 + ht2

# Normalize all values between 0 and 1
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
quantile(p_norm, c(0.01, 0.5, 0.99))
p_norm <- range01(p_dis)
t_norm <- range01(tc)
k_norm <- range01(as.matrix(k))
quantile(p_dis, c(0.01, 0.5, 0.99))
quantile(tc, c(0.01, 0.5, 0.99))
quantile(as.matrix(k), c(0.01, 0.5, 0.99))

quantile(p_norm, c(0.01, 0.5, 0.99))
quantile(t_norm, c(0.01, 0.5, 0.99))
quantile(k_norm, c(0.01, 0.5, 0.99))

ht1 <- Heatmap(t_norm, name='Exp Correlation',  show_row_names = F, show_column_names = F,
               col=colorRamp2(c(0.1619918,0.3582600,0.5220454), c('blue', 'white', 'red')))
htp <- Heatmap(p_norm, name='Phenotype Distance', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order,
               col=colorRamp2(c(0.01038832,0.15763476,0.60928360), c('red', 'white', 'blue'))) 
ht2 <- Heatmap(k_norm, name='Kinship', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order,
               col=colorRamp2(c(0.0338904,0.0836813,0.3148306), c('blue', 'white', 'red'))) 
ht1 + htp + ht2
ht1 <- Heatmap(t_norm, name='Exp Correlation',  show_row_names = F, show_column_names = F,
               col=colorRamp2(c(0.08830086,0.35826004,1), c('blue', 'white', 'red')))
htp <- Heatmap(p_norm, name='Phenotype Distance', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order,
               col=colorRamp2(c(0,0.1576348,0.7955771), c('red', 'white', 'blue'))) 
ht2 <- Heatmap(k_norm, name='Kinship', cluster_columns=FALSE, show_row_names = F, show_column_names = F, column_order = ht1_order,
               col=colorRamp2(c(0.0168835,0.0836813,0.9312361), c('blue', 'white', 'red'))) 
ht1 + htp + ht2

## Plot PCC (eCor ~ Kinship) (supplemental figure)
# 6/11/18
k$line <- row.names(k)
k.m <- melt(k, id='line')
names(k.m) <- c('line1', 'line2', 'kinship')
t.cor <- as.data.frame(tc)
t.cor$line <- row.names(t.cor)
t.m <- melt(t.cor, id='line')
names(t.m) <- c('line1', 'line2', 'eCor')
kt.m <- merge(k.m, t.m, by=c('line1', 'line2') )
kt.m.s <- subset(kt.m, line1 != line2)

eCor.kin.PCC <- cor(kt.m$kinship, kt.m$eCor)
#summary(lm(kinship~eCor, data=kt.m))
#summary(lm(eCor~kinship, data=kt.m))
ggplot(kt.m.s, aes(x=kinship, y=eCor)) + geom_point(size=0.2, alpha=0.2) + 
  theme_bw(8) +
  geom_smooth(method='lm',formula=y~x, se=FALSE)
ggsave('00_Figures/180703_eCorKin_Scatter.pdf', width=3.42, height = 3.42, device='pdf', useDingbats=FALSE)


# Bin the values
binned <- kt.m.s[]
binned$k_bin <- cut(binned$kinship, breaks = c(seq(-0.22100, 1.92200, by = .15)))
binned$bin <- as.character(lapply(strsplit(as.character(binned$k_bin), split=","),head, n=1))
binned$bin <- gsub("\\(", "", binned$bin)

binned$bin2 <- as.numeric(binned$bin)
quants <- as.data.frame(do.call("rbind", tapply(binned$eCor, binned$bin2, quantile, c(0.05, 0.95))))
names(quants) <- c('x5','x95')
quants$bin2 <- as.numeric(row.names(quants))

ggplot(data = binned, aes(x=bin2, y=eCor, group=bin2))  + 
  geom_boxplot() + theme_bw(10) +
  geom_line(data=quants, aes(x=bin2, y=x95, group = 1), color='red') +
  geom_line(data=quants, aes(x=bin2, y=x5, group = 1), color='blue')
ggsave('00_Figures/180703_eCorKin_Distribution.pdf', width=3.42, height = 3.42, device='pdf', useDingbats=FALSE)

binned_count <- aggregate(cbind(count = line2) ~ bin2, data = binned, 
                          FUN = function(x){NROW(x)})
ggplot(binned_count, aes(x=bin2, y=log(count))) + theme_bw(10) +
  geom_bar(stat='identity', position= "dodge")
ggsave('00_Figures/180705_eCorVSkin_BinCount.pdf', width=3.42, height = 2.42, device='pdf', useDingbats=FALSE)


### Correlation between kinship/eCor and pheno distance
p_dis <-as.data.frame(as.matrix(distances(p)))
diag(p_dis) <- min(p_dis)
names(p_dis) <- row.names(p)
p_dis$line <- row.names(p)

p.m <- melt(p_dis, id='line')
names(p.m) <- c('line1', 'line2', 'dist')

kp.m <- merge(p.m, k.m, by=c('line1', 'line2') )
kp.m.s <- subset(kp.m, line1 != line2)
tp.m <- merge(p.m, t.m, by=c('line1', 'line2') )
tp.m.s <- subset(tp.m, line1 != line2)

cor.test(kp.m.s$kinship, kp.m.s$dist, method= 'spearman')
cor.test(tp.m.s$eCor, tp.m.s$dist, method= 'spearman')
cor.test(kt.m.s$kinship, kt.m.s$eCor, method= 'spearman')


# Bin the values for kinship vs. phenotype
summary(kp.m.s$kinship)
kp.m.s$k_bin <- cut(kp.m.s$kinship, breaks = c(seq(-0.22100, 1.92200, by = .15)))
kp.m.s$bin <- as.character(lapply(strsplit(as.character(binned$k_bin), split=","),head, n=1))
kp.m.s$bin <- gsub("\\(", "", kp.m.s$bin)

kp.m.s$bin2 <- as.numeric(kp.m.s$bin)
quants <- as.data.frame(do.call("rbind", tapply(kp.m.s$dist, kp.m.s$bin2, quantile, c(0.05, 0.95))))
names(quants) <- c('x5','x95')
quants$bin2 <- as.numeric(row.names(quants))

ggplot(data = kp.m.s, aes(x=bin2, y=dist, group=bin2))  + 
  geom_boxplot() + theme_bw(10) +
  geom_line(data=quants, aes(x=bin2, y=x95, group = 1), color='red') +
  geom_line(data=quants, aes(x=bin2, y=x5, group = 1), color='blue')
ggsave('00_Figures/180705_PhenoDistKin_Distribution.pdf', width=3.42, height = 3.42, device='pdf', useDingbats=FALSE)

binned_count <- aggregate(cbind(count = line2) ~ bin2, data = kp.m.s, 
                          FUN = function(x){NROW(x)})
ggplot(binned_count, aes(x=bin2, y=log(count))) + theme_bw(10) +
  geom_bar(stat='identity', position= "dodge")
ggsave('00_Figures/180705_PhenoDistKin_BinCount.pdf', width=3.42, height = 2.42, device='pdf', useDingbats=FALSE)

# Bin the values for kinship vs. phenotype
summary(tp.m.s$eCor)
tp.m.s$k_bin <- cut(tp.m.s$eCor, breaks = c(seq(0.8939, 0.9857, by = .006)))
tp.m.s$bin <- as.character(lapply(strsplit(as.character(tp.m.s$k_bin), split=","),head, n=1))
tp.m.s$bin <- gsub("\\(", "", tp.m.s$bin)

tp.m.s$bin2 <- as.numeric(tp.m.s$bin)
quants <- as.data.frame(do.call("rbind", tapply(tp.m.s$dist, tp.m.s$bin2, quantile, c(0.05, 0.95))))
names(quants) <- c('x5','x95')
quants$bin2 <- as.numeric(row.names(quants))

ggplot(data = tp.m.s, aes(x=bin2, y=dist, group=bin2))  + 
  geom_boxplot() + theme_bw(10) +
  geom_line(data=quants, aes(x=bin2, y=x95, group = 1), color='red') +
  geom_line(data=quants, aes(x=bin2, y=x5, group = 1), color='blue')
ggsave('00_Figures/180705_PhenoDisteCor_Distribution.pdf', width=3.42, height = 3.42, device='pdf', useDingbats=FALSE)

binned_count <- aggregate(cbind(count = line2) ~ bin2, data = tp.m.s, 
                          FUN = function(x){NROW(x)})
ggplot(binned_count, aes(x=bin2, y=log(count))) + theme_bw(10) +
  geom_bar(stat='identity', position= "dodge")
ggsave('00_Figures/180705_PhenoDisteCor_BinCount.pdf', width=3.42, height = 2.42, device='pdf', useDingbats=FALSE)


# Get expression variance across genes
colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}
# Convert from loge to log10 and get variance for each gene
medi <- apply(log2(exp(t)), 2, FUN=median)
summary(medi)
vari <- colVars(log2(exp(t)))
summary(vari)

