#####################################################################################################################################################################
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
#paralogs <- paralogs[order(paralogs$in_group),]
#####################################################################################################################################################################
#categorize the expression level
df_TPM <- TPM
df_TPM[df_TPM >= 10] = 10
df_TPM[df_TPM >= 2 & df_TPM < 10 ] = 2
df_TPM[df_TPM < 2] = 0
#####################################################################################################################################################################
#change paralog dataframe into single gene per row

paralogs_df <- paralogs %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
paralogs_df[,1] <- gsub("\\.[^.]*$","",paralogs_df[,1])
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
#paralogs_df["Solyc01g007810.1",]
#rownames(paralogs_df) <- paralogs_df$in_group
#####################################################################################################################################################################
#prase TPM values
paralogs_TPM <- df_TPM[paralogs_df$in_group,]

#add group column to the paralogs_TPM
#paralogs_TPM$Group <- paralogs_df$Group
#save paraologs_TPM as a file
#write.table(paralogs_TPM, "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/modling_in_python/paralogs_TPM_with_groups.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#Here we are tring to look at the expression of the homologous gene pairs. 
#1. First compare the number of tissues a gene expressed between gene pair and calculate frequncy of combinaions of tissue expression 
#   among all the gene pairs frequency of combinations of tissue expression 

#     1a.getting basic distribution for frequncy of combinaions of tissue expression 

mat_pair_exp <- data.frame(matrix(ncol = ncol(df_TPM)+1, nrow = ncol(df_TPM)+1))
colnames(mat_pair_exp) <- 0:ncol(df_TPM)
rownames(mat_pair_exp) <- 0:ncol(df_TPM)
mat_pair_exp[is.na(mat_pair_exp)] <- 0
for (group in unique(paralogs_df$Group)) {
  TPM_sub_df <- paralogs_TPM[rownames(paralogs_df[paralogs_df[,2] == group,]),]
  TPM_sub_df[TPM_sub_df >= 2] = 1
  if (nrow(TPM_sub_df) == 2) {
    mat_pair_exp[as.character(rowSums(TPM_sub_df[1,])),as.character(rowSums(TPM_sub_df[2,]))] <- mat_pair_exp[as.character(rowSums(TPM_sub_df[1,])),as.character(rowSums(TPM_sub_df[2,]))] +1
  }
}

#calculate sum of all the values in the matrix

(sum(mat_pair_exp[1:5,19:24]) + sum(mat_pair_exp[19:24,1:5]) )/sum(mat_pair_exp)
394/sum(mat_pair_exp)
sum(mat_pair_exp[1,1])/sum(mat_pair_exp)
#####################################################################################################################################################################
#Drawing heatmap for the basic distribution
#col_fun = colorRamp2(c( 0, max(mat_pair_exp)), c("blue", "deeppink2"))
#col_fun(seq(-2, 2))
#col_fun = colorRamp2(c(0,50), hcl_palette = 'Mako', reverse = TRUE)
col_fun = colorRamp2(c(0,25,50,75,100), c("#ffffff", "#3da41a", "#0e98e3", "#77059a", "#060c50"))
HM <- Heatmap(as.matrix(mat_pair_exp), name = "Count", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(23:0),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "number of tissues that gene 1 is expresssed",column_title = "number of tissues that gene 2 is expresssed")
pdf("expression_frequncy_of_gene_paires_in_different_tissues_anters_removed.pdf",height = 10,width = 10)
HM
dev.off()
#####################################################################################################################################################################
#     2a.simulating background distribution for the frequncy


sim_df_list <- list()
itrs <- 1000
for (i in 1:itrs) {
  mat_pair_sim<- data.frame(matrix(ncol = ncol(df_TPM)+1, nrow = ncol(df_TPM)+1))
  #mat_pair_sim<- matrix(ncol = ncol(df_TPM)+1, nrow = ncol(df_TPM)+1)
  colnames(mat_pair_sim) <- 0:ncol(df_TPM)
  rownames(mat_pair_sim) <- 0:ncol(df_TPM)
  mat_pair_sim[is.na(mat_pair_sim)] <- 0
  for (group in unique(paralogs_df$Group)) {
    TPM_sub_df <- paralogs_TPM[rownames(paralogs_df[paralogs_df[,2] == group,]),]
    TPM_sub_df[TPM_sub_df >= 2] = 1
    if (nrow(TPM_sub_df) == 2) {
      sim_df <- sapply(1:ncol(TPM_sub_df), function (row) TPM_sub_df[,row] <- sample(TPM_sub_df[,row]))
      mat_pair_sim[as.character(sum(sim_df[1,])),as.character(sum(sim_df[2,]))] <- mat_pair_sim[as.character(sum(sim_df[1,])),as.character(sum(sim_df[2,]))] +1
      #z = ((as.numeric(rowSums(TPM_sub_df[1,])) -  as.numeric(rowSums(TPM_sub_df[2,]))) - median(sim_exp_pairs))/sd(sim_exp_pairs)
      #mat_pair_Z_score[as.character(rowSums(TPM_sub_df[1,])),as.character(rowSums(TPM_sub_df[2,]))] <- z
    }
  }
  sim_df_list[[i]] <- mat_pair_sim
}
length(sim_df_list)
z_val_mat <- data.frame(matrix(ncol = ncol(df_TPM)+1, nrow = ncol(df_TPM)+1))
#mat_pair_sim<- matrix(ncol = ncol(df_TPM)+1, nrow = ncol(df_TPM)+1)
colnames(z_val_mat) <- 0:ncol(df_TPM)
rownames(z_val_mat) <- 0:ncol(df_TPM)
z_val_mat[is.na(z_val_mat)] <- 0

for (row_pos in 1:ncol(df_TPM)+1) {
  for (col_pos in 1:ncol(df_TPM)+1) {
    res <- sapply(sim_df_list, function(x) x[[row_pos,col_pos]])
    act_val <- mat_pair_exp[row_pos,col_pos]
    z <- (act_val - mean(res))/sd(res)
    z_val_mat[row_pos,col_pos] <- z
  }
  
}

#test
row_pos <- 1
col_pos <- 23
res <- sapply(sim_df_list, function(x) x[[row_pos,col_pos]])
act_val <- mat_pair_exp[row_pos,col_pos]
z <- (act_val - mean(res))/sd(res)


z_val_mat[0,1]
z_val_mat[sapply(z_val_mat, is.infinite)] <- NA
mean(res)
sd(res)
  
  
  
  
sim_df_list[[1]][1,1]
as.data.frame(mat_pair_sim)
#unique(sim_exp_pairs)
#plot( as.factor(sim_exp_pairs))
#axis(1, at=c(-24:24),labels=c(-24:24))
#####################################################################################################################################################################
#Drawing heatmap for the basic distribution
#col_fun = colorRamp2(c( 0, max(mat_pair_exp)), c("blue", "deeppink2"))
#col_fun(seq(-2, 2))less
col_fun = colorRamp2(c(min(z_val_mat,na.rm = TRUE),5), hcl_palette = 'Tropic', reverse = TRUE)
HM_Z <- Heatmap(as.matrix(z_val_mat), name = "Z_score", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "gray",
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(23:0),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "number of tissues that gene 1 is expresssed",column_title = "number of tissues that gene 2 is expresssed")
pdf("expression_frequncy_of_gene_paires_in_different_tissues_after_simulation.pdf",height = 10,width = 10)
HM_Z
dev.off()  
  
#####################################################################################################################################################################
#     3a.calculate probebilities for each cell in the frequncy matrix
  