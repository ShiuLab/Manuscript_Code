setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')
library(data.table)

### Format each dataset individually

# Phenotype data
# Get mean phenotype values across years
library(dplyr)
pheno <- fread('WiDiv_Phenotypic_Data_08_09.txt', sep='\t', na.strings=".")
pheno <- subset(pheno, select=c(Entry, plant_height, GDD, stover_yield))
p_agg <- pheno %>%
  group_by(Entry) %>% 
  summarise_all(funs(mean(., na.rm = TRUE)))
names(p_agg) <- c('ID','HT','FT','YLD')
write.table(p_agg, file='WiDiv_Phenotypic_Data_08_09.txt_meansCA.csv', sep=",", quote = FALSE, row.names=FALSE)

# Transcriptome data
# Processing done by Jeremy
# Data to use: maize503_logPlus1FPKM_B73setFiltered.txt
trans <- fread('maize503_logPlus1FPKM_B73setFiltered.txt', sep='\t')

# Genotype data
# 1. save GAPIT.hmp with only the lines wer're going to use in the study
# 2. Convert to [-1,0,1] format
# 3. Convert AGPv2 SNP names to AGPv4 & remove SNPs that don't map to v4

geno <- fread('GAPIT.RNAseq.hmp_438K_imputed2.csv', sep=',')

### 1. Save GAPIT.hmp with only the lines we're going to use in this study
g_names <- c('rs','alleles','chrom','pos','strand','assembly','center','protLSID','assayLSID','panel','Qccode')
g_names <- c(g_names, as.vector(key$geno))
geno_hmp <- subset(geno, select=g_names)
write.table(geno_hmp, file='GAPIT.RNAseq.hmp_438K_imputed2_keepCA.csv', sep=",", quote = FALSE, row.names=FALSE)


### 2. Convert genotypes from HapMap to [-1,0,1] format
genot <- geno[]
genot$alleles <- substr(genot$alleles, 1,1)
genot$allelesX <- paste(genot$alleles, genot$alleles, by="")
genot$allelesX <- gsub(" ","", genot$allelesX)

# Replace allele 1 with 1 and allele 2 with -1 
genot[,2:ncol(genot)] <- as.data.frame(lapply(genot[,2:ncol(genot)], function(x) ifelse(x == genot$allelesX, 1, -1)))
genot <- subset(genot, select = -c(alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panel, Qccode))
genot <- na.omit(genot)

# Transpose
genott <- as.data.frame(t(genot[,-1]))
colnames(genott) <- as.character(unlist(genot$rs))
genott <- genott[!rownames(genott) %in% c('alleles', 'allelesX'), ]

# Make sure dominant allele is 1 and minor allele is -1
genott['Sum', ] <- colSums(genott)
count <- 0
for(i in names(genott)){
  if(sum(as.numeric(unlist(genott[i]))) < 0){
    genott[i] <- genott[i] * -1
    count <- count + 1
    if(count %% 10000==0) {
      # Print on the screen some message
      print(paste0("progress: ", count, "\n"))
    }
  }
}


# Remove SNPs that have less than 5% minor allele frequency
genott2 <- genott[!rownames(genott) %in% c('Sum'), ]
geno95 <- genott2[,colMeans(genott2) < 0.95 & colMeans(genott2) > -0.95]
summary(genott2[,'rna1_509'])

# Geno file in [-1,0,1] format with only SNPs with >5% allele frequency
#fwrite(geno95, file = "GAPIT.RNAseq.hmp_438K_imputed2_AlleleConvertedCA.csv", append = FALSE, quote = "auto", sep = ",", row.names=TRUE)
write.table(geno95, 'GAPIT.RNAseq.hmp_438K_imputed2_AlleleConvertedCA.csv', quote = FALSE, row.names=TRUE, col.names=TRUE, append=FALSE)



### 3. Convert AGPv2 SNP names to AGPv4
library(tidyr)
geno95 <- fread('GAPIT.RNAseq.hmp_438K_imputed2_AlleleConvertedCA.csv', sep=',', header=T) 
snp_key <- fread('zm_v4_503_snp_w_v2_pos.sort.header.txt', sep='\t', header=T)
snp_key$id <- paste(snp_key$CHROM_V2, snp_key$POS_V2, sep='_')

geno_names <- as.data.frame(names(geno95))
geno_names <- as.data.frame(geno_names[2:nrow(geno_names),])
names(geno_names) <- 'id'

geno_names$id2 <- gsub("rna", "chr", geno_names$id)
g_names <- separate(data = geno_names, col = id2, into = c("CHROM_V2", "POS_V2"), sep = "_")
g_names$id <- gsub("rna", "chr", g_names$id)

names_m <- merge(snp_key, g_names, by='id')
names_m$V4_name <- paste(names_m$CHROM_V4, names_m$POS_V4, sep='_')
names_m$V2_name <- gsub('chr', 'rna', names_m$id)

genot <- as.data.frame(t(geno95))
geno_newnames <- merge(genot, names_m[,c('V2_name','V4_name')], by.x='row.names', by.y='V2_name')

geno_newnamesT <- as.data.frame(t(geno_newnames))
names(geno_newnamesT) <- as.character(unlist(geno_newnamesT['V4_name',]))
geno_newnamesT <- geno_newnamesT[!rownames(geno_newnamesT) %in% c('V4_name', 'Row.names'), ]
rownames(geno_newnamesT) <- rownames(geno95)

write.table(geno_newnamesT, 'GAPIT.RNAseq.hmp_438K_imputed2_AlleleConverted_V4_CA.csv', sep=",", quote = FALSE, row.names=TRUE)




# Merge all three data types
# Need to standardize the line names and get the 388 lines w/ complete data
g <- fread('GAPIT.RNAseq.hmp_438K_imputed2_AlleleConverted_V4_CA.csv', sep=',')
t <- fread('maize503_logPlus1FPKM_B73setFiltered.txt', sep='\t')
p <- fread('WiDiv_Phenotypic_Data_08_09.txt_meansCA.csv', sep=',')
k <- read.table('../01_Data/kinship.txt', sep='\t', skip=3, row.names=1)

key <- read.csv('maize_lineID_key.txt', header=T, sep='\t', na='')

# Phenotypes
p_named <- merge(key, p, by.x = 'pheno', by.y='ID')
p_named[,c('transcripto', 'pheno', 'geno')] <- NULL
write.table(p_named, file='../01_Data/pheno.csv', sep=",", quote = FALSE, row.names=FALSE)

# Transcripto
t_named <- merge(key, t, by.x = 'transcripto', by.y='ID')
t_named[,c('transcripto', 'pheno', 'geno')] <- NULL
write.table(t_named, file='../01_Data/transcripto.csv', sep=",", quote = FALSE, row.names=FALSE)

# Geno
g_named <- merge(key, g, by.x = 'geno', by.y='row.names')
g_named[,c('transcripto', 'pheno', 'geno')] <- NULL
write.table(g_named, file='../01_Data/geno.csv', sep=",", quote = FALSE, row.names=FALSE)

# Kinship matrix (made by Jeremy)
k_named <- merge(key, k, by.x = 'geno', by.y='row.names')
k_named[,c('transcripto', 'pheno', 'geno')] <- NULL
write.table(k_named, file='../01_Data/kin.csv', sep=",", quote = FALSE, row.names=FALSE, col.names=TRUE)



# Make randomized transcriptome dataframe
t <- fread('01_Data/transcripto.csv', sep=',', header=T)
row.names(t) <- t$ID
t$ID <- NULL
t.shuf <- t[sample(nrow(t)),] # shuffle row-wise (i.e. genes have same exp values, just shuffled)
row.names(t.shuf) <- row.names(t)
write.table(t.shuf, file='01_Data/transcripto_shuff.csv', sep=',', quote=F, row.names=T, col.names=T)

# Make eC matrix-- expression correlation matrix
t.cor <- as.data.frame(cor(t(t)))
row.names(t.cor) <- row.names(t)
write.table(t.cor, file='01_Data/transcripto_cor.csv', sep=',', quote=F, row.names=T, col.names=T)



######################################
## Format HMP321 SNPs (genome wide) ##
######################################

library(data.table)
g <- fread('hmp321_withDPGL_imputed_numeric_NoDups_merged.csv')
g <- as.data.frame(g)
g_names <- fread('c10_hmp321_withDPGL_imputed.vcf.gz_Numerical.txt')
g$V1 <- NULL

g$ID <- unlist(g_names[,"<Marker>"])

key <- read.csv('../../00_RawData/maize_lineID_key_withHMP321.txt', header=T, sep='\t', na='')
g_named <- merge(key, g, by.x = 'HMP321', by.y='ID')
g_named[,c('transcripto', 'pheno', 'geno', 'HMP321')] <- NULL
row.names(g_named) <- g_named$ID
g_named$ID <- NULL

p <- read.csv('hm_pheno.csv', sep=',', row.names=1)
p2 <- p[match(rownames(g_named), rownames(p)),]

write.table(g_named, file='../hmp_geno.csv', sep=',', quote=F, row.names=T, col.names=T)
write.table(p2, file='../hmp_pheno.csv', sep=',', quote=F, row.names=T, col.names=T)

# Format transcripts from lines with HMP data only
hmp_p <- read.csv('01_Data/hmp_pheno.csv', sep=',',header=T, row.names=1)
t <- fread('01_Data/transcripto.csv')
t <- as.data.frame(t)
row.names(t) <- t$ID
t$ID <- NULL
hmp_t <- t[match(rownames(hmp_p), rownames(t)),]
write.table(hmp_t, file='01_Data/hmp_transcripto.csv', sep=',', quote=F, row.names=T, col.names=T)


