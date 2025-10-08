# Description: Estimate the narrow-sense heritability of diploid S. cerevisiae
# isolates using biallelic SNPs data

### 0. Install necessary packages
#install.packages('sommer') # v.4.3.0 (in R v4.0.3)


### 1. Load packages and data
library(sommer)

geno <- read.csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/geno.csv", header=T, row.names=1) # Peter 2018 diploid S. cerevisiae isolate biallelic SNPs data
pheno <- read.csv("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", header=T, row.names=1) # Peter 2018 diploid S. cerevisiae isolate fitness data

### 2. Compute relationship matrices
A <- A.mat(as.matrix(geno)) # additive relationship matrix
D <- D.mat(as.matrix(geno)) # dominance relationship matrix
# E <- E.mat(as.matrix(geno)) # epistatic relationship matrix

### 3. Fit the mixed model 
traits <- colnames(pheno)
results <- data.frame(Conditions=traits,h2=0,h2_SE=0) #,H2_ADE=0,H2_ADE_SE=0, H2_AD=0, H2_AD_SE=0)
for (i in 1:length(traits)){
  p <- pheno[i] # subset pheno matrix
  p$ID <- rownames(geno)
  p$IDD <- rownames(geno)
  p$IDE <- rownames(geno)
  colnames(p) <- c("trt", "ID", "IDD", "IDE")
  head(p)
  
  ### 4. Estimate the narrow-sense heritability
  ans.ADE <- mmer(trt~1, random=~vs(ID,Gu=A) + vs(IDD,Gu=D), rcov=~units, data=p, verbose = FALSE)
  summary(ans.ADE)$varcomp
  h2 <- vpredict(ans.ADE, h2 ~ (V1) / (V1+V3))
  results$h2[i] <- h2$Estimate[1] # narrow sense heritability
  results$h2_SE[i] <- h2$SE[1]
}

### 5. Write to file
write.csv(results, "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Data/Peter_2018/Heritability_h2_H2_sommer_CORRECTED.csv")
