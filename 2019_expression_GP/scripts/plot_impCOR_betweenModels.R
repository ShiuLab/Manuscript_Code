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
#install.packages('ggtern')
#library(ggtern)

# Read in arguments with 5 as the default PCs to include
if (length(args)==0) {
  stop("Need to provide 2 arguments (output from filter_SNPs_LD.py and save name)", call.=FALSE)
}

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


# Read in importance scores for genotype and transcriptome data

save <- args[1]
rrB_file = args[2]
BL_file = args[3]
RF_file = args[4]

#coef_file <- 'RF_FT_SNP2kb.csv'
print(save)
print('Reading in coefficients...')
rrB <- read.table(rrB_file, sep=',', head=T, row.names = 1, na.strings = "na")
BL <- read.table(BL_file, sep=',', head=T, row.names = 1, na.strings = "na")
RF <- read.table(RF_file, sep=',', head=T, row.names = 1, na.strings = "na")

print('Number of transcripts with SNP available within window size:')
print(sum(!is.na(rrB$SNP)))


rrB <- rrB[!is.na(rrB$SNP_coef),c('v3_gene_model','SNP_coef', 'gene_coef')]
names(rrB) <- c('gene','rrBLUP_SNP', 'rrBLUP_trans')
BL <- BL[!is.na(BL$SNP_coef),c('v3_gene_model','SNP_coef', 'gene_coef')]
names(BL) <- c('gene','BL_SNP', 'BL_trans')
RF <- RF[!is.na(RF$SNP_coef),c('v3_gene_model','SNP_coef', 'gene_coef')]
names(RF) <- c('gene','RF_SNP', 'RF_trans')


data <- merge(rrB, BL, by='gene')
data <- merge(data, RF, by='gene')
data$d_G_rrB_BL <- get_density(data$rrBLUP_SNP, data$BL_SNP, n=100)
data$d_T_rrB_BL <- get_density(data$rrBLUP_trans, data$BL_trans, n=100)
data$d_G_rrB_RF <- get_density(data$rrBLUP_SNP, data$RF_SNP, n=100)
data$d_T_rrB_RF <- get_density(data$rrBLUP_trans, data$RF_trans, n=100)
data$d_G_BL_RF <- get_density(data$RF_SNP, data$BL_SNP, n=100)
data$d_T_BL_RF <- get_density(data$RF_trans, data$BL_trans, n=100)

print(head(data))


# rrBLUP vs. BL
cor <- cor(data$rrBLUP_SNP, data$BL_SNP, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=rrBLUP_SNP, y=BL_SNP, color=d_G_rrB_BL),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_rrB_BL_SNPs.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)

cor <- cor(data$rrBLUP_trans, data$BL_trans, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=rrBLUP_trans, y=BL_trans, color=d_T_rrB_BL),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_rrB_BL_trans.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)



# rrBLUP vs. RF
cor <- cor(data$rrBLUP_SNP, data$RF_SNP, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=rrBLUP_SNP, y=RF_SNP, color=d_G_rrB_RF),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_rrB_RF_SNPs.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)

cor <- cor(data$rrBLUP_trans, data$RF_trans, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=rrBLUP_trans, y=RF_trans, color=d_T_rrB_RF),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_rrB_RF_trans.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)


# BL vs RF
cor <- cor(data$BL_SNP, data$RF_SNP, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=BL_SNP, y=RF_SNP, color=d_G_BL_RF),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_BL_RF_SNPs.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)

cor <- cor(data$BL_trans, data$RF_trans, use='complete.obs', method ='spearman')
ggplot(data) + geom_point(aes(x=BL_trans, y=RF_trans, color=d_T_BL_RF),shape=16, alpha=0.75, size=1) + 
  scale_color_viridis() + theme_bw() +
  ggtitle(paste('ro = ', round(cor, 3))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste(save, '_BL_RF_trans.pdf', sep=''), width=3, height = 2.5, device='pdf', useDingbats=FALSE)


print('finished!')


