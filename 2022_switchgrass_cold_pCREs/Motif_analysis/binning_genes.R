#This is a code to combine all the genes in upregulation sets and binning them to subsets for mapping
#####################################################################################################################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#####################################################################################################################################################################
min_30 <- read.csv("upregulated_genes_30min.txt", header = F, sep ="\t")
hr_1 <- read.csv("upregulated_genes_1hr.txt", header = F, sep ="\t")
hr_3 <- read.csv("upregulated_genes_3hr.txt", header = F, sep ="\t")
hr_6 <- read.csv("upregulated_genes_6hr.txt", header = F, sep ="\t")
hr_16 <- read.csv("upregulated_genes_16hr.txt", header = F, sep ="\t")
hr_24 <- read.csv("upregulated_genes_24hr.txt", header = F, sep ="\t")
non_res <- read.csv("nonresponsive_genes.txt", header = F, sep ="\t")

all_genes <- c(min_30$V1,hr_1$V1,hr_3$V1,hr_6$V1,hr_16$V1,hr_24$V1, non_res$V1)
all_genes <- unique(all_genes)

all_gene_25_bins <- split(all_genes,             # Applying split() function
      ceiling(seq_along(all_genes) / 25))
x =1
for (list_bin in all_gene_25_bins) {
  #print(list_bin)
  df <- as.data.frame(list_bin)
  write.table(as.data.frame(list_bin), paste(paste("gene_set",x,sep = "_"),"txt",sep = "."),col.names = F, quote = F, row.names = F)
  x = x +1
}

