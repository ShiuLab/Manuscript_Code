rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(7)
##########################Libraries##############################
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
###########################Data##########################################################################################################################################
TPM <- read.csv("TPM_table_non_T_tissues.txt",row.names=1, header = TRUE, sep = "\t")
colnames(TPM) 
#remove Mature.anthers
TPM <- TPM[,-which(colnames(TPM) == "Mature.anthers")]
dff_TPM <- read.table('Tissue_cluster_groups_TPM.txt',head=T,sep='\t',stringsAsFactors=F)
dff_nont_T <- dff_TPM[match(colnames(TPM) , dff_TPM[,1]),]
#TPM["Solyc00g014290.1.1",]
paralogs_with_out <- read.table("../alignmnts_guide_with_groups.txt",sep = "\t",header = TRUE)
paralogs <- paralogs_with_out[,c(1,3)]
head(paralogs)
#paralogs <- paralogs[order(paralogs$in_group),]
#####################################################################################################################################################################
#categorize the expression level
df_TPM <- TPM
df_TPM[df_TPM >= 10] = 10
df_TPM[df_TPM >= 2 & df_TPM < 10 ] = 2
df_TPM[df_TPM < 2] = 0

#change paralog dataframe into single gene per row

paralogs_df <- paralogs %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
paralogs_df[,1] <- gsub("\\.[^.]*$","",paralogs_df[,1])
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
#paralogs_df["Solyc01g007810.1",]
#rownames(paralogs_df) <- paralogs_df$in_group
head(paralogs_df)
#####################################################################################################################################################################
#prase TPM values
paralogs_TPM <- df_TPM[paralogs_df$in_group,]
head(paralogs_TPM)

#####################################################################################################################################################################
gene_not_exp_both <-list()
gene_exp_both <-list()
gene_exp_oly_one <-list()
single_exp_gene <-list()

for (group in unique(paralogs_df$Group)) {
    #group = 26
  TPM_sub_df <- paralogs_TPM[rownames(paralogs_df[paralogs_df[,2] == group,]),]
  TPM_sub_df[TPM_sub_df >= 2] = 1
  if (nrow(TPM_sub_df) == 2) {
    #mat_pair_exp[as.character(rowSums(TPM_sub_df[1,])),as.character(rowSums(TPM_sub_df[2,]))] <- mat_pair_exp[as.character(rowSums(TPM_sub_df[1,])),as.character(rowSums(TPM_sub_df[2,]))] +1
    #check if both genes are expressed in all the tissues
    if (sum(rowSums(TPM_sub_df) == ncol(TPM_sub_df)) == 2) {
      #append gene names to the list
        gene_exp_both[[group]] <- rownames(TPM_sub_df)
    } else if (sum(rowSums(TPM_sub_df) == 0) == 2) {
      gene_not_exp_both[[group]] <- rownames(TPM_sub_df)
    }
    #check if one paralog is expressed in all the tissues while the other one is  not expressed in any of the tissues
    else if (sum(rowSums(TPM_sub_df) == ncol(TPM_sub_df)) == 1 & sum(rowSums(TPM_sub_df) == 0) == 1) {
      gene_exp_oly_one[[group]] <- rownames(TPM_sub_df)
    }
    #heck if one paralog is expressed only in one tissues while the other one is  not expressed in any of the tissues
    else if (sum(rowSums(TPM_sub_df) == 1) == 1 & sum(rowSums(TPM_sub_df) == 0) == 1) {
      single_exp_gene[[group]] <- rownames(TPM_sub_df)
    }
  }
}

#get all the gene names into a single list
gene_exp_both <- unlist(gene_exp_both)
length(gene_exp_both)
gene_not_exp_both <- unlist(gene_not_exp_both)
length(gene_not_exp_both)
gene_exp_oly_one <- unlist(gene_exp_oly_one)
length(gene_exp_oly_one)
single_exp_gene <- unlist(single_exp_gene)
length(single_exp_gene)

#spot checking
df_TPM[gene_exp_oly_one[4],]


#write the gene names to a file without quotes
write.table(gene_exp_both, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_expressed_in_all_tissues.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gene_not_exp_both, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_all_tissues.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gene_exp_oly_one, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/only_one_paralog_expressed_in_all_tissues.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(single_exp_gene, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/one_paralog_expressed_in_one_tissue_while_other_non.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#######################
#Checking if the genes that does not expresse in any tissues express in the tissues we did not check

TPM_OG <- read.csv("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/TMP_for_all_tissues.txt",row.names=1, header = TRUE, sep = "\t")
head(TPM_OG)
#remove the columns that are in TPM data frame
TPM_non_looked <- TPM_OG[,!colnames(TPM_OG) %in% colnames(TPM)]
head(TPM_non_looked)
#check if the genes that are not expressed in any tissues are expressed in the tissues we did not check

TPM_non_exp_genes <- TPM_non_looked[gene_not_exp_both,]
head(TPM_non_exp_genes)
#save the data frame
write.table(TPM_non_exp_genes, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_all_tissues_not_exp_in_any_tissues.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#get the genes whhoes expression is less than 2 in all the tissues
subdat_all <- TPM_non_exp_genes
subdat_all[subdat_all < 2] = 0
subdat_all[subdat_all >= 2] = 1

genes_not_exp_any_tissues <- rownames(subdat_all[rowSums(subdat_all) == 0,])
length(genes_not_exp_any_tissues)
length(gene_not_exp_both)
#get the genes that are not common in both lists
genes_exp_in_OG_not_in_23<- genes_not_exp_any_tissues[!genes_not_exp_any_tissues %in% gene_not_exp_both]
length(genes_exp_in_OG_not_in_23)
#save the data frame
write.table(genes_not_exp_any_tissues, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_al_46_different_tissues.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(genes_exp_in_OG_not_in_23, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_tissues_contain_trichome_but_not_in_23.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
