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
#library('data.table')
# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
start.time <- Sys.time()
#install.packages('ggtern')
#library(ggtern)

# Read in arguments with 5 as the default PCs to include
if (length(args)==0) {
  stop("Need to provide 2 arguments (output from filter_SNPs_LD.py and save name)", call.=FALSE)
}
#setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/06_CompareCoefImp/')


# Read in importance scores for genotype and transcriptome data

coef_file <- args[1]
eqtl_file <- args[2]
save_name <- args[3]

print(save_name)
print('Reading in coefficients...')
imp <- read.table(coef_file, sep=',', head=T, row.names = 1, na.strings = "na")
print(head(imp))

print('Number of transcripts with SNP available within transcript region:')
print(sum(!is.na(imp$SNP)))
imp <- imp[!is.na(imp$SNP_coef),]

# Read in eQTL pairs and subset distance pairs to have transcripts with eQTL only
eqtl <- read.table(eqtl_file, sep=',', head=T, row.names = 1, na.strings = "na")
imp_dist_eq <- subset(imp, v3_gene_model %in% eqtl$v3_gene_model)

cutoff <- 0.99
g.per95.d <- quantile(imp_dist_eq$SNP_coef, cutoff, na.rm=TRUE)
t.per95.d <- quantile(imp_dist_eq$gene_coef, cutoff, na.rm=TRUE)
g.per95.e <- quantile(eqtl$eqtl_imp, cutoff, na.rm=TRUE)
t.per95.e <- quantile(eqtl$gene_imp, cutoff, na.rm=TRUE)

# Calculate correlation scores
print('eQTL results')

corr_eqtl = cor(eqtl$eqtl_imp, eqtl$gene_imp, use='complete.obs', method ='spearman')
corr_dist = cor(imp_dist_eq$SNP_coef, imp_dist_eq$gene_coef, use='complete.obs', method ='spearman')

imp_dist_eq <- imp_dist_eq[,c('SNP_coef', 'gene_coef')]
imp_dist_eq$type <- 'distance'
eqtl2 <- eqtl[,c('eqtl_imp', 'gene_imp')]
eqtl2$type <- 'eQTL'
names(eqtl2) <- c('SNP_coef', 'gene_coef', 'type')

eqtl_snp <- rbind(imp_dist_eq, eqtl2)

# Generate Density plots

library(ggplot2)
library(MASS)
library(viridis)

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

eqtl_snp$density <- get_density(eqtl_snp$gene_coef, eqtl_snp$SNP_coef, n=100)
eqtl_snp_thresh <- data.frame(type=c('distance', 'eQTL'), SNP_thresh=c(g.per95.d, g.per95.e), gene_thresh=c(t.per95.d,t.per95.e))
print(max(eqtl_snp$gene_coef))

ggplot(eqtl_snp) + geom_point(aes(x=gene_coef, y=SNP_coef, color=density),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + facet_grid(type~.) + theme_bw() +
  geom_hline(data=eqtl_snp_thresh, aes(yintercept=SNP_thresh), color='red', linetype='dashed') + 
  geom_vline(data=eqtl_snp_thresh, aes(xintercept=gene_thresh), color='red', linetype='dashed') +
  ggtitle(paste('eQTL = ', round(corr_eqtl, 3),'dist = ', round(corr_dist, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save_name, '_corplot.pdf', sep=''), width=3, height = 4.2, device='pdf', useDingbats=FALSE)
