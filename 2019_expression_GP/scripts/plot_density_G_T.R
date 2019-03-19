### Plot SNP vs. transcript density across the genome
library(data.table)
library(ggplot2)
setwd('/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/00_RawData')

# Transcript density
genes <- fread('v3_v4_xref.txt', sep='\t')
genes$v4_chr <- gsub('Chr', '', genes$v4_chr)
genes$v4_chr <- with(genes, factor(v4_chr, levels=c(1:10), ordered=TRUE))

genes <- genes[!is.na(genes$v4_chr)]
genes <- genes[!is.na(genes$v4_start)]
genes <- genes[!duplicated(genes$v4_gene_model),]

print('Calculating cumulative position...')
genes$cumBP <- genes$v4_start
cumu <- max(subset(genes, v4_chr == 1)$v4_start)
genes$cumBP[genes$v4_chr == 2] <- genes$cumBP[genes$v4_chr == 2] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 2)$v4_start)
genes$cumBP[genes$v4_chr == 3] <- genes$cumBP[genes$v4_chr == 3] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 3)$v4_start)
genes$cumBP[genes$v4_chr == 4] <- genes$cumBP[genes$v4_chr == 4] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 4)$v4_start)
genes$cumBP[genes$v4_chr == 5] <- genes$cumBP[genes$v4_chr == 5] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 5)$v4_start)
genes$cumBP[genes$v4_chr == 6] <- genes$cumBP[genes$v4_chr == 6] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 6)$v4_start)
genes$cumBP[genes$v4_chr == 7] <- genes$cumBP[genes$v4_chr == 7] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 7)$v4_start)
genes$cumBP[genes$v4_chr == 8] <- genes$cumBP[genes$v4_chr == 8] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 8)$v4_start)
genes$cumBP[genes$v4_chr == 9] <- genes$cumBP[genes$v4_chr == 9] + cumu
cumu <- cumu + max(subset(genes, v4_chr == 9)$v4_start)
genes$cumBP[genes$v4_chr == 10] <- genes$cumBP[genes$v4_chr == 10] + cumu


# Density plot
print('Plotting...')
ggplot(genes) + geom_histogram(aes(x=cumBP, color=v4_chr), binwidth=1000000) + 
  xlab('Genomic Position (1Mb bins)') +
  scale_y_continuous(expand = c(0, 0) ) +  
  scale_color_manual(values = rep(c("gray30", "gray60"), 10)) +
  theme_bw() + 
  theme(legend.position="none", text=element_text(size=8, color='black'),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

#save_name1 <- paste(save_name, '.pdf', sep='')
ggsave('Density_Transcripts.pdf', width = 14, height = 2, useDingbats=FALSE)



# For SNPs
snps <- fread('zm_v4_503_snp_w_v2_pos.sort.header.txt', sep='\t')
snps$CHROM_V4 <- gsub('Chr', '', snps$CHROM_V4)
snps <- subset(snps, CHROM_V4 %in% c(1:10))
print(head(snps))
snps$CHROM_V4 <- with(snps, factor(CHROM_V4, levels=c(1:10), ordered=TRUE))

snps <- snps[!is.na(snps$CHROM_V4)]
snps <- snps[!is.na(snps$POS_V4)]

print('Calculating cumulative position...')
snps$cumBP <- snps$POS_V4
cumu <- max(subset(snps, CHROM_V4 == 1)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 2] <- snps$cumBP[snps$CHROM_V4 == 2] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 2)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 3] <- snps$cumBP[snps$CHROM_V4 == 3] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 3)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 4] <- snps$cumBP[snps$CHROM_V4 == 4] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 4)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 5] <- snps$cumBP[snps$CHROM_V4 == 5] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 5)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 6] <- snps$cumBP[snps$CHROM_V4 == 6] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 6)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 7] <- snps$cumBP[snps$CHROM_V4 == 7] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 7)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 8] <- snps$cumBP[snps$CHROM_V4 == 8] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 8)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 9] <- snps$cumBP[snps$CHROM_V4 == 9] + cumu
cumu <- cumu + max(subset(snps, CHROM_V4 == 9)$POS_V4)
snps$cumBP[snps$CHROM_V4 == 10] <- snps$cumBP[snps$CHROM_V4 == 10] + cumu


# Density plot
print('Plotting...')
ggplot(snps) + geom_histogram(aes(x=cumBP, color=CHROM_V4), binwidth=1000000) + 
  xlab('Genomic Position (1Mb bins)') +
  scale_y_continuous(expand = c(0, 0) ) +  
  scale_color_manual(values = rep(c("gray30", "gray60"), 10)) +
  theme_bw() + 
  theme(legend.position="none", text=element_text(size=8, color='black'),
        panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())

#save_name1 <- paste(save_name, '.pdf', sep='')
ggsave('Density_SNPs.pdf', width = 14, height = 2, useDingbats=FALSE)

