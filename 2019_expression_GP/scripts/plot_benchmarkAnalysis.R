### Plots for FT benchmark genes
library('ggplot2')
library('reshape')
library('data.table')
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')

ft <- read.csv('09_Candidate_Genes_CA/FT_benchmark_info.txt', sep='\t', header=T)
ftkeep <- ft[,grepl('_Z',names(ft))]
ftkeep$V3_gene <- ft$V3_gene
ftm <- melt(ftkeep, by=c('V3_gene'))

ftm$alg <- sapply(strsplit(as.character(ftm$variable),'_'), "[", 1)
ftm$type <- sapply(strsplit(as.character(ftm$variable),'_'), "[", 2)

reo_Alg <- function(x) { factor(x, levels = c('rrB','BL','RF', 'En'))}
gene_order <- c('GRMZM2G011357','GRMZM2G067921','GRMZM2G171365','GRMZM2G179264','GRMZM2G700665','GRMZM2G171650','GRMZM2G440005','AC217051.3_FG006','GRMZM2G075081','GRMZM2G381691')
gene_order <- c('GRMZM2G381691','GRMZM2G075081','AC217051.3_FG006','GRMZM2G440005','GRMZM2G171650','GRMZM2G700665','GRMZM2G179264','GRMZM2G171365','GRMZM2G067921','GRMZM2G011357')

reo_genes <- function(x) { factor(x, levels = gene_order)}
ftm$valueCO <- ftm$value
ftm$valueCO[ftm$valueCO>10] <- 10
ftm$valueCO[ftm$V3_gene == 'GRMZM2G075081' & ftm$type== 'G'] <- -10

ggplot(ftm) + geom_tile(aes(x=reo_Alg(alg), y=reo_genes(V3_gene), fill=valueCO), colour='white') +
  facet_grid(.~type) + theme_bw(12) + 
  scale_fill_gradient2(low='#0571b0', mid='white', high='#ca0020', midpoint=0)
#ggsave("180816_FT_benchmark_Z.pdf", width = 4, height = 4,useDingbats=FALSE )



#### Plot allele type by flowering time for benchmark genes
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/09_Candidate_Genes_CA/')
snps <- c('Chr1_243201519','Chr7_181107438','Chr9_156980141','Chr3_160591848','Chr8_126879558',
          'Chr7_143622186','Chr8_136009724','Chr9_59721','Chr7_85620028','Chr1_273326886',
          'Chr4_24380636')
d <- fread('geno_FTbench2.csv', header=T)
y <- as.data.frame(fread('../01_Data/pheno.csv'))
d2 <- as.data.frame(d)
d2 <- merge(d2, y[c('ID','FT')], by.x='V1', by.y='ID', all=F)
dg <- melt(d2, id=c('V1','FT'))
dg$type <- 'G'
dg$value <- as.factor(dg$value)

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
ggplot(dg, aes(x=value, y=FT)) +
  geom_violin() +
  geom_boxplot(width=.3, outlier.size = NA, coef = 0) + 
  stat_summary(fun.data = give.n, geom = "text") +
  facet_grid(.~variable, scales='free') +
  theme_bw(10)
ggsave('190311_ViolinBoxPlots_SNPs.pdf', width = 5, height = 2,useDingbats=FALSE)


stats_ttest <- data.frame(snp = character(), pval = character(), stringsAsFactors=FALSE)
for (sn in snps){
  dg_temp <- subset(dg, variable == sn)
  t_temp <- t.test(FT ~ value, dg_temp)
  stats_ttest <- rbind(stats_ttest, list(snp = sn, pval = t_temp$p.value), stringsAsFactors=F)
}


#### Plot expression level by flowering time for benchmark genes
trans <- c('GRMZM2G011357','GRMZM2G067921','GRMZM2G075081','GRMZM2G171365','GRMZM2G171650',
           'GRMZM2G179264','GRMZM2G381691','GRMZM2G440005','GRMZM2G700665','AC217051.3_FG006',
           'GRMZM2G004483', 'GRMZM2G156079','GRMZM2G026223','GRMZM2G106903')
d <- fread('../01_Data/transcripto.csv', header=T)
d <- as.data.frame(d)

d2 <- d[trans]
d2$ID <- row.names(d)
d2 <- merge(d2, y[c('ID','FT')], by.x='ID', by.y='ID', all=F)
dt <- melt(d2, id=c('ID','FT'))
dt$value2 <- as.numeric(as.character(dt$value))
dt$FT <- as.numeric(as.character(dt$FT))
dt$ft_bin <- cut(dt$FT, breaks = c(seq(500, 1300, by = 50)))
dt$bin <- as.character(lapply(strsplit(as.character(dt$ft_bin), split=","),head, n=1))
dt$bin <- gsub("\\(", "", dt$bin)
dt$bin2 <- as.numeric(dt$bin)

for(t in trans){
  dt2 <- subset(dt, variable==t)
  model_temp <- summary(lm(FT ~ value, dt2))
  r2 <- model_temp$r.squared
  pval <- model_temp$coefficients[,4][2]
  quants <- as.data.frame(do.call("rbind", tapply(dt2$value2, dt2$bin2, quantile, c(0.05, 0.95))))
  names(quants) <- c('x5','x95')
  quants$bin2 <- as.numeric(row.names(quants))
  ggplot(data = dt2, aes(x=bin2, y=value2, group=bin2))  + 
    geom_boxplot(outlier.shape = NA) + theme_bw(10) + 
    xlab('FT (GDD)') + ylab('log(FC)') + ggtitle(paste(r2, pval, sep=' : p='))+ 
    geom_line(data=quants, aes(x=bin2, y=x95, group = 1), color='red') +
    geom_line(data=quants, aes(x=bin2, y=x5, group = 1), color='blue') 
  ggsave(paste('180920_TransVSft_', t, '.pdf'), width = 2, height = 2,useDingbats=FALSE)
}



### Old code for density plots

for (tr in trans){
  dt_temp <- subset(dt, variable == tr)
  model_temp <- summary(lm(FT ~ value, dt_temp))
  r2 <- model_temp$r.squared
  pval <- model_temp$coefficients[,4][2]
  temp_title <- paste(r2, '  ', pval, sep='')
  save_name <- paste('denPlot_', tr, '.pdf', sep='')
  pdf(save_name, width=4.5, height=4)
  
  smoothScatter(dt_temp$value, dt_temp$FT, xlab='Expression Level', ylab='FT (GDD)', useRaster=TRUE)
  title(temp_title)
  dev.off()
}

for (tr in trans){
  dt_temp <- subset(dt, variable == tr)
  model_temp <- summary(lm(FT ~ value, dt_temp))
  r2 <- model_temp$r.squared
  pval <- model_temp$coefficients[,4][2]
  temp_title <- paste(r2, '  ', pval, sep='')
  save_name <- paste('denPlot_', tr, '.pdf', sep='')
  pdf(save_name, width=4.5, height=4)
  smoothScatter(dt_temp$value, dt_temp$FT, xlab='Expression Level', ylab='FT (GDD)', useRaster=TRUE)
  title(temp_title)
  dev.off()
}


###################
# Plot percentile #
###################
d <- read.csv('09_Candidate_Genes_CA/use_scores.csv')
d <- read.csv('12_eQTL//use_scores.csv')
key <- read.csv('09_Candidate_Genes_CA/use_key.txt', sep='\t')

name_order <- c('mads69','zag6','mads1','pebp24','pebp8','id1','cct2','rap2','cct1','pebp4','pebp5',
                'dlf1','pebp2','unk')

name_order <- c('unk','pebp2','dlf1', 'pebp5','pebp4','cct1','rap2','cct2','id1',
                'pebp8','pebp24','mads1','zag6','mads69')
reo_na <- function(x) { factor(x, levels = name_order)}
reo_Alg <- function(x) { factor(x, levels = c('rrBLUP','BL','RF', 'Ensemble'))}

d2 <- merge(d, key, by='Name', all.x=TRUE)
d2$perc2 <- d2$perc
d2$perc2[d2$perc2 <= 0.5] <- 0.5
ggplot(d2, aes(x=reo_Alg(alg), y=reo_na(Name))) + 
  geom_tile(aes(fill = perc), color='white') +
  facet_grid(.~Type.x) + theme_bw(12) +
  scale_fill_gradientn(colors = c('white','#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026'),
                         values = c(0,0,0.5,0.8,0.9,0.95,0.99,1))
ggsave('00_Figures/190308_FTBM_heatmap.pdf', width = 4, height = 3,useDingbats=FALSE)



######### Background t-test and lm expectations ##########
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/09_Candidate_Genes_CA/')
bg_t <- read.csv('T_P_lm_results.txt', sep='\t')

p99 <- quantile(-log10(bg_t$pval), 0.99, na.rm=T)
p999 <- quantile(-log10(bg_t$pval), 0.999, na.rm=T)
bg_t$log10pval <- -log10(bg_t$pval)
ggplot(bg_t, aes(x=-log10(pval))) + geom_histogram(fill='skyblue', color='gray50', binwidth=0.5) + scale_y_log10() +
  stat_bin(geom="text", size=3.5, aes(label=..count.., y=1.2*(..count..))) + theme_bw() + 
  geom_vline(xintercept = p99) + geom_vline(xintercept = p999)
ggsave('../00_Figures/190315_FT_background_T.pdf', width = 4, height = 3,useDingbats=FALSE)

bg_g <- read.csv('G_P_lm_results.txt', sep='\t')

p99 <- quantile(-log10(bg_g$pval), 0.99, na.rm=T)
p999 <- quantile(-log10(bg_g$pval), 0.999, na.rm=T)
bg_g$log10pval <- -log10(bg_g$pval)
ggplot(bg_g, aes(x=-log10(pval))) + geom_histogram(fill='skyblue', color='gray50', binwidth=0.25) + scale_y_log10() +
  stat_bin(geom="text", size=3.5, aes(label=..count.., y=1.2*(..count..))) + theme_bw() + 
  geom_vline(xintercept = p99) + geom_vline(xintercept = p999)
ggsave('../00_Figures/190315_FT_background_G.pdf', width = 4, height = 3,useDingbats=FALSE)


