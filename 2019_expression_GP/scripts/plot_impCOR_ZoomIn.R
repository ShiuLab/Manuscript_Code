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
library('data.table')

setwd('/Users/cazodi/Desktop/')
#library('data.table')
# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args)==0) {
  stop("Need to provide 2 arguments (output from filter_SNPs_LD.py and save name)", call.=FALSE)
}
setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')


# Read in importance scores for genotype and transcriptome data
g_file <- args[1]
t_file <- args[2]
zoom_to <- args[3]
save_name <- args[4]
t_file <- '03_rrBLUP/rrBLUP_transcripto_FT_coef.csv'

print(save_name)
print('Reading in coefficients...')
g_imp <- fread('ensem_geno_FT_imp_stats.csv', sep=',', head=T, na.strings = "na")
t_imp <- fread('ensem_transcripto_FT_imp_stats.csv', sep=',', head=T, na.strings = "na")
key <- fread('v3_v4_xref.txt', header = T)
key <- key[,c('v4_chr', 'v4_start','v4_end','v3_gene_model')]
key <- unique(key)

# Process genotype data
g_imp$v4_chr <- as.character(lapply(strsplit(as.character(g_imp$ID), split='_'), head, n=1))
g_imp$loc <- as.numeric(lapply(strsplit(as.character(g_imp$ID), split='_'), tail, n=1))
t_imp$type <- 'G'

# Process transcript data
t_imp$type <- 'T'
t_imp <- merge(t_imp, key, by.x='ID', by.y='v3_gene_model', all.x=T, all.y=F)
t_imp$loc <-  t_imp$v4_start + (abs(t_imp$v4_start - t_imp$v4_end)/2)

# Get top X percent
cutoff <- 0.99
t.per95 <- quantile(t_imp$zscore, cutoff)
g.per95 <- quantile(g_imp$zscore, cutoff)
t_imp$thresh <- t.per95
g_imp$thresh <- g.per95

# Top ensemble SNP: Chr3_160559109
target_C <- 'Chr3'
target_loc <- 160559109

# Top ensemble transcript- mads69: GRMZM2G171650
target_C <- 'Chr3'
target_loc <- 160591595

# Second transcript- GRMZM5G865543
target_C <- 'Chr2'
target_loc <- 206591808

# Second SNP- Chr8_136012624
target_C <- 'Chr8'
target_loc <- 136012624

window <- 1000000

t_imp_s <- subset(t_imp, v4_chr == target_C & loc > target_loc-window & loc < target_loc+window)
g_imp_s <- subset(g_imp, v4_chr == target_C & loc > target_loc-window & loc < target_loc+window)

test_c <- 'Chr7'
test_loc <- 85621313

test_w <- 1500
g_imp_s <- subset(g_imp, v4_chr == test_c & loc > test_loc-test_w & loc < test_loc+test_w)


print('Number of transcript and SNP within window size:')
print(dim(t_imp_s)[1])
print(dim(g_imp_s)[1])

print('Plotting results...')
# Generate manhattan plots

ggplot(g_imp_s, aes(x=loc, y=zscore)) + theme_bw() +
  geom_point(color='red', size=1) +
  geom_hline(data=g_imp_s, aes(yintercept = thresh), color='red', linetype='dashed') +
  geom_point(data=t_imp_s, aes(x=v4_start, y=zscore), color='green', size=1, alpha=0.7) +
  geom_point(data=t_imp_s, aes(x=v4_end, y=zscore), color='blue', size=1, alpha=0.7) +
  geom_hline(data=t_imp_s, aes(yintercept = thresh), color='green', linetype='dashed') +
  geom_hline(yintercept = 0, color='black') +
  theme(legend.position="none", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())
#ggsave('180920_MPzsco_Chr3160559109.pdf', width = 4, height = 2.5, useDingbats=FALSE)
ggsave('180920_MPzsco_GRMZM5G865543.pdf', width = 4, height = 2.5, useDingbats=FALSE)
ggsave('180920_MPzsco_Chr8136012624.pdf', width = 4, height = 2.5, useDingbats=FALSE)


ggplot(g_imp_s, aes(x=loc, y=-log(1-perc))) + theme_bw() +
  geom_point(color='red', size=1) +
  geom_point(data=t_imp_s, aes(x=v4_start, y=-log(1-perc)), color='green', size=1, alpha=0.7) +
  geom_point(data=t_imp_s, aes(x=v4_end, y=-log(1-perc)), color='blue', size=1, alpha=0.7) +
  geom_hline(yintercept = -log(1-0.99), color='red', linetype='dashed') +
  geom_hline(yintercept = -log(1), color='black') +
  theme(legend.position="none", panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank())

#ggsave('180920_MPperc_Chr3160559109.pdf', width = 4, height = 2.5, useDingbats=FALSE)
#ggsave('180920_MPper_GRMZM5G865543.pdf', width = 4, height = 2.5, useDingbats=FALSE)
#ggsave('180920_MPper_Chr8136012624.pdf', width = 4, height = 2.5, useDingbats=FALSE)


print('finished!')


