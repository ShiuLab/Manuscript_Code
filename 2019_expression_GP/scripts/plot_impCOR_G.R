#####################################
# Importance Analysis from Genomic Prediction analysis
# Comparing top SNP within the desired LD (see: /mnt/research/ShiuLab/17_GP_SNP_Exp/scripts/filter_SNPs_LD.py)
# to transcript importance scores
#
#     arg1 = output from filter_SNPs_LD.py
#     arg2 = save_name (script adds .pdf)
#
# Written by: Christina Azodi
# Original: 6.14.18
#####################################
library('ggplot2')
library('reshape')
library('MASS')
library('viridis')

# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
start.time <- Sys.time()


# Read in arguments with 5 as the default PCs to include
if (length(args)==0) {
  stop("Need to provide 2 arguments (output from filter_SNPs_LD.py and save name)", call.=FALSE)
}

# Read in importance scores for genotype and transcriptome data
coef_file <- args[1]
save_name <- args[2]

print(save_name)
print('Reading in coefficients...')
imp <- read.table(coef_file, sep=',', head=T, row.names = 1, na.strings = "na")
print(head(imp))
print('Number of transcripts with SNP available within window size:')
print(sum(!is.na(imp$SNP)))

imp <- imp[!is.na(imp$SNP_coef),]
imp$Chromo <- gsub('Chr','', imp$Chromo)

# Get 99th percentile 
cutoff <- 0.99
g.per95 <- quantile(imp$SNP_coef, cutoff, na.rm=TRUE)
t.per95 <- quantile(imp$gene_coef, cutoff, na.rm=TRUE)
print('Cutoff for transcript (coral)')
print(t.per95)
print('Cutoff for SNPs (green)')
print(g.per95)

data <- melt(imp, id.vars = c("SNP",'gene_mid', 'v3_gene_model', 'Chromo'), measure.vars = c("SNP_coef", "gene_coef"))
data$thresh <- g.per95
data$thresh[data$variable == 'gene_coef'] <- t.per95

print('Calculating cumulative position...')
data$cumBP <- data$gene_mid
cumu <- max(subset(data, Chromo == '1')$gene_mid)
data$cumBP[data$Chromo == '2'] <- data$cumBP[data$Chromo == '2'] + cumu
cumu <- cumu + max(subset(data, Chromo == '2')$gene_mid)
data$cumBP[data$Chromo == '3'] <- data$cumBP[data$Chromo == '3'] + cumu
cumu <- cumu + max(subset(data, Chromo == '3')$gene_mid)
data$cumBP[data$Chromo == '4'] <- data$cumBP[data$Chromo == '4'] + cumu
cumu <- cumu + max(subset(data, Chromo == '4')$gene_mid)
data$cumBP[data$Chromo == '5'] <- data$cumBP[data$Chromo == '5'] + cumu
cumu <- cumu + max(subset(data, Chromo == '5')$gene_mid)
data$cumBP[data$Chromo == '6'] <- data$cumBP[data$Chromo == '6'] + cumu
cumu <- cumu + max(subset(data, Chromo == '6')$gene_mid)
data$cumBP[data$Chromo == '7'] <- data$cumBP[data$Chromo == '7'] + cumu
cumu <- cumu + max(subset(data, Chromo == '7')$gene_mid)
data$cumBP[data$Chromo == '8'] <- data$cumBP[data$Chromo == '8'] + cumu
cumu <- cumu + max(subset(data, Chromo == '8')$gene_mid)
data$cumBP[data$Chromo == '9'] <- data$cumBP[data$Chromo == '9'] + cumu
cumu <- cumu + max(subset(data, Chromo == '9')$gene_mid)
data$cumBP[data$Chromo == '10'] <- data$cumBP[data$Chromo == '10'] + cumu


print('Plotting results...')
# Generate manhattan plots
data$Chromof = factor(data$Chromo, levels=c('1','2','3','4','5','6','7','8','9','10'))
ggplot(data, aes(x=cumBP, y=value)) + theme_bw(10) +
  geom_point( aes(color=Chromof), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("gray30", "gray60"), 22 )) +
  geom_hline(data=data, aes(yintercept = thresh), color='red', linetype='dashed') +
  facet_grid(variable~., scales = "free") +
  geom_point(data=subset(data, value > thresh), color="blue", size=0.8) +
  #scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  
  theme_bw() + 
  theme(legend.position="none", text=element_text(size=8, color='black'),
    panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank())


save_name1 <- paste(save_name, '.pdf', sep='')
ggsave(save_name1, width = 7, height = 3, useDingbats=FALSE)


# Pull tran:SNP pairs that were both in the top 1%
both_top <- subset(imp, (SNP_coef > g.per95 & gene_coef > t.per95))

print('')
print('Transcripts where both Transcript & associated SNP were in the top 1%')
print(dim(both_top))
print(both_top$v3_gene_model)
print('')

print(head(imp[order(-imp$SNP_coef),]))
print(head(imp[order(-imp$gene_coef),]))

# Generate correlation plots
slope <- max(imp$SNP_coef, na.rm=TRUE)/max(imp$gene_coef, na.rm=TRUE)
correlation = cor(imp$SNP_coef, imp$gene_coef, use='complete.obs', method ='spearman')

# Get density of points in 2 dimensions.
theme_set(theme_bw(base_size = 16))
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
set.seed(1)
dat <- data.frame(
  x = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0, sd = 0.1)
  ),
  y = c(
    rnorm(1e4, mean = 0, sd = 0.1),
    rnorm(1e3, mean = 0.1, sd = 0.2)
  )
)

imp$density <- get_density(imp$gene_coef, imp$SNP_coef, n=100)
ggplot(imp) + geom_point(aes(x=gene_coef, y=SNP_coef, color=density),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  geom_hline(yintercept=g.per95, color='red', linetype='dashed') + 
  geom_vline(xintercept=t.per95, color='red', linetype='dashed') +
  ggtitle(paste('ro = ', round(correlation, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save_name, '_cor.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)


print(coef_file)
print('ro:')
print(correlation)
print('Number of transcripts with SNP available within window size:')
print(sum(!is.na(imp$SNP)))
print('finished!')


