#THIS IS THE UPDATED VERSION OF THE SCRIPT deviation_between_ancestral_vs_extant_vals.R 
#See the notes in Ancestral state interpritation to see the changes made to the script
#_____________________________________________________________________________________________________________________#
#This is a ssript to calculate the deviation between the ancestral and extant states for 
#1. each gene in the datasets
#2. for the syntanic homelogs that constain two genes
#This is basically replicating the analysis on Fig  of Panchy et al 2019 paper
########################################################################################
#####################################################################################################################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##########################Libraries##############################
##########################Libraries##############################
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
##########################Functions##############################
################################################################
#get the upper triangle of the partition_table# Get lower triangle of the correlation matrix
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}
#####################################################################################################################################################################
TPM <- read.csv("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/TPM_table_non_T_tissues.txt",row.names=1, header = TRUE, sep = "\t")
dff_TPM <- read.table('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Tissue_cluster_groups_TPM.txt',head=T,sep='\t',stringsAsFactors=F)
dff_nont_T <- dff_TPM[match(colnames(TPM) , dff_TPM[,1]),]
#TPM["Solyc00g014290.1.1",]
paralogs_with_out <- read.table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Run_7_modified_alignmnts_guide_with_key_info_with_groups.txt",sep = "\t",header = TRUE)
head(paralogs_with_out)
paralogs_with_out <- paralogs_with_out[paralogs_with_out[,"toplogy"]=="abide",c(1,3)]

ancestral_vals <- read.table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Trichome_ancestral_state_of_genome_T_paralogs.txt",sep = "\t",header = TRUE)
head(ancestral_vals)
#remove 'Group' from  from group column
ancestral_vals$group <- gsub("Group","",ancestral_vals$group)
#####################################################################################################################################################################
#categorize the expression level
high_exp_thresh <- 10

df_TPM <- TPM
df_TPM[df_TPM >= high_exp_thresh] = high_exp_thresh
df_TPM[df_TPM >= 2 & df_TPM < high_exp_thresh ] = 2
df_TPM[df_TPM < 2] = 0
#####################################################################################################################################################################
#change paralog dataframe into single gene per row
#Solyc01g110440.5.1
paralogs_df <- paralogs_with_out %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
head(paralogs_df)
paralogs_df[paralogs_df[,1]== 'Solyc02g038817.1.1',]
paralogs_df[paralogs_df[,2]== 2117 | paralogs_df[,2]== 1261,]
paralogs_df[,1] <- gsub("\\.[^.]*$","",paralogs_df[,1])
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)

#chenge ancestral_vals dataframe into single gene per row
ancestral_vals_df <- ancestral_vals %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
ancestral_vals_df <- as.data.frame(ancestral_vals_df, row.names = ancestral_vals_df$in_group)
ancestral_vals_df <- ancestral_vals_df[,c(2,1,3)]
#ancestral_vals_df[,1] <- gsub("\\.[^.]*$","",ancestral_vals_df[,1])
head(ancestral_vals_df)
#cataegorize the expression level in the ancestral state
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= high_exp_thresh] = high_exp_thresh
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= 2 & ancestral_vals_df$ancestral_state < high_exp_thresh ] = 2
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state < 2] = 0
#ancestral_vals_df[,2] <- gsub("\\.[^.]*$","",ancestral_vals_df[,2])
head(ancestral_vals_df)

#####################################################################################################################################################################
#get extarnt state TPM values for in_group genes in ancestral_vals_df in a new column
colnames(df_TPM)   
tissue <- "ZZTrichome"
#rownames(ancestral_vals_df) <- ancestral_vals_df$in_group
#ancestral_vals_df[ancestral_vals_df[,1] ==  'Solyc02g038817.1.1',]
#ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",rownames(ancestral_vals_df)),tissue]
ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",ancestral_vals_df$in_group),tissue]
head(ancestral_vals_df)
df_TPM['Solyc01g109320.4',tissue]

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#calculate the deviation between the ancestral and extant states for each gene paire
partition_table_0 <- data.frame(matrix(NA,nrow = 5,ncol = 5))

#make colnames and rownames for the partition_table to be between -2 and 2
colnames(partition_table_0) <- c(-2,-1,0,1,2)
rownames(partition_table_0) <- c(-2,-1,0,1,2)
partition_table_0[is.na(partition_table_0)] <- 0
#deep copy the partition_table_0
partition_table_1 <- partition_table_0
partition_table_2 <- partition_table_0


for (group in unique(ancestral_vals_df$group)) {
    #group <- 200
    ancestral_vals_sub_df <- ancestral_vals_df[ancestral_vals_df$group == group,]
    
    if (nrow(ancestral_vals_sub_df) == 2) {
        #order rows by extant state
        ancestral_vals_sub_df <- ancestral_vals_sub_df[order(ancestral_vals_sub_df$extant_state,decreasing = TRUE),]
        #chaning the ancestral state values to 0,1,2
        ancestral_vals_sub_df$ancestral_state[ancestral_vals_sub_df$ancestral_state == 2] = 1
        ancestral_vals_sub_df$extant_state[ancestral_vals_sub_df$extant_state == 2] = 1
        ancestral_vals_sub_df$ancestral_state[ancestral_vals_sub_df$ancestral_state == high_exp_thresh] = 2
        ancestral_vals_sub_df$extant_state[ancestral_vals_sub_df$extant_state == high_exp_thresh] = 2
        #new column to store the deviation between the ancestral and extant states
        ancestral_vals_sub_df$deviation <- ancestral_vals_sub_df$extant_state - ancestral_vals_sub_df$ancestral_state
        #count the frequency of each deviation
        if (unique(ancestral_vals_sub_df$ancestral_state) == 0) {
            if(length(unique(ancestral_vals_sub_df$deviation)) == 1){
               partition_table_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            }else {           
            partition_table_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            partition_table_0[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- partition_table_0[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
            }
        }else if (unique(ancestral_vals_sub_df$ancestral_state) == 1) {
            if(length(unique(ancestral_vals_sub_df$deviation)) == 1){
               partition_table_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            }else {           
            partition_table_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            partition_table_1[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- partition_table_1[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
            }
        }else if (unique(ancestral_vals_sub_df$ancestral_state) == 2) {
            if(length(unique(ancestral_vals_sub_df$deviation)) == 1){
               partition_table_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            }else {           
            partition_table_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- partition_table_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            partition_table_2[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- partition_table_2[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
            }
        }
        
    
    }
}

up_partition_table_0 <- get_upper_tri(partition_table_0)
#get sum of all numbers in the upper triangle
sum(up_partition_table_0,na.rm = T)
#presentage retention
sum(up_partition_table_0[3,3],na.rm = T)/sum(partition_table_0,na.rm = T)


up_partition_table_2 <- get_upper_tri(partition_table_2)
#get sum of all numbers in the upper triangle
sum(up_partition_table_2,na.rm = T)
#presentage retention
sum(up_partition_table_2[3,3],na.rm = T)/sum(up_partition_table_2,na.rm = T)

up_partition_table_1 <- get_upper_tri(partition_table_1)
#get sum of all numbers in the upper triangle
96/sum(up_partition_table_1,na.rm = T)

#log10 partition_table
l10_partition_table <- log10(partition_table)
upper_tri <- get_upper_tri(l10_partition_table)
#make infinities to NA
upper_tri[upper_tri == -Inf] <- NA
col_fun = colorRamp2(c(0,max(upper_tri,na.rm = T)), hcl_palette = 'Mint', reverse = TRUE)
HM <- Heatmap(as.matrix(upper_tri), name = "log_10_frequency", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "white",
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(-2:2),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "Duplicate_2_deviation",column_title = "Duplicate_1_deviation"
              #add numbers to the heatmap
              )
pdf("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Figures/Ancestral_deviation_true_Young.seeds.pdf",height = 10,width = 10)
HM
dev.off()                         




#simulation of the deviation between the ancestral and extant states for the syntanic homelogs that constain two genes

sim_df_list_0 <- list()
sim_df_list_1 <- list()
sim_df_list_2 <- list()
iter <- 10000
for(i in 1:iter){
    sim_df_1 <- data.frame(matrix(NA,nrow = 5,ncol = 5))
    colnames(sim_df_1) <- c(-2,-1,0,1,2)
    rownames(sim_df_1) <- c(-2,-1,0,1,2)
    sim_df_1[is.na(sim_df_1)] <- 0
    sim_df_0 <- sim_df_1
    sim_df_2 <- sim_df_1
    for (group in unique(ancestral_vals_df$group)) {
        #group <- 210
        ancestral_vals_sub_df <- ancestral_vals_df[ancestral_vals_df$group == group,]
        
        if (nrow(ancestral_vals_sub_df) == 2) {

            #replace the extant state values with random values
            ancestral_vals_sub_df$extant_state <- sample(c(0,2,high_exp_thresh),2,replace = TRUE)
            #replace the ancestral state values with random values
            #ancestral_vals_sub_df$ancestral_state <- sample(c(0,2,10),2,replace = TRUE)
            #chaning the ancestral state values to 0,1,2
            ancestral_vals_sub_df$ancestral_state[ancestral_vals_sub_df$ancestral_state == 2] = 1
            ancestral_vals_sub_df$extant_state[ancestral_vals_sub_df$extant_state == 2] = 1
            ancestral_vals_sub_df$ancestral_state[ancestral_vals_sub_df$ancestral_state == high_exp_thresh] = 2
            ancestral_vals_sub_df$extant_state[ancestral_vals_sub_df$extant_state == high_exp_thresh] = 2
            #new column to store the deviation between the ancestral and extant states
            ancestral_vals_sub_df$deviation <- ancestral_vals_sub_df$extant_state - ancestral_vals_sub_df$ancestral_state
            #count the frequency of each deviation
            #sim_df[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
            #sim_df[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- sim_df[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
            if (unique(ancestral_vals_sub_df$ancestral_state) == 0) {
                if(length(unique(ancestral_vals_sub_df$deviation)) == 1){
                   sim_df_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                }else {           
                sim_df_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_0[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                sim_df_0[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- sim_df_0[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
                }
            }else if (unique(ancestral_vals_sub_df$ancestral_state) == 1) {
                if(length(unique(ancestral_vals_sub_df$deviation)) == 1){
                   sim_df_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                }else {           
                sim_df_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_1[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                sim_df_1[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- sim_df_1[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
                }
            }else if (unique(ancestral_vals_sub_df$ancestral_state) == 2) {
                if (length(unique(ancestral_vals_sub_df$deviation)) == 1){
                   sim_df_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                }else {
                sim_df_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] <- sim_df_2[as.character(ancestral_vals_sub_df[2,"deviation"]),as.character(ancestral_vals_sub_df[1,"deviation"])] + 1
                sim_df_2[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] <- sim_df_2[as.character(ancestral_vals_sub_df[1,"deviation"]),as.character(ancestral_vals_sub_df[2,"deviation"])] + 1
                }
            }
        }
    } 
    sim_df_list_0[[i]] <- sim_df_0
    sim_df_list_1[[i]] <- sim_df_1
    sim_df_list_2[[i]] <- sim_df_2  
}

length(sim_df_list_0)
#calculate the z score for each deviation
z_val_mat_0 <- data.frame(matrix(NA,nrow = 5,ncol = 5))
colnames(z_val_mat_0) <- c(-2,-1,0,1,2)
rownames(z_val_mat_0) <- c(-2,-1,0,1,2)
z_val_mat_0[is.na(z_val_mat_0)] <- 0
z_val_mat_1 <- z_val_mat_0
z_val_mat_2 <- z_val_mat_0

for(col_pos in 1:5){
    for(row_pos in 1:5){
        res_0 <- sapply(sim_df_list_0, function(x) x[row_pos,col_pos])
        act_val_0 <- partition_table_0[row_pos,col_pos]
        z_0 <- (act_val_0 - mean(res_0))/sd(res_0)
        z_val_mat_0[row_pos,col_pos] <- z_0

        res_1 <- sapply(sim_df_list_1, function(x) x[row_pos,col_pos])
        act_val_1 <- partition_table_1[row_pos,col_pos]
        z_1 <- (act_val_1 - mean(res_1))/sd(res_1)
        z_val_mat_1[row_pos,col_pos] <- z_1

        res_2 <- sapply(sim_df_list_2, function(x) x[row_pos,col_pos])
        act_val_2 <- partition_table_2[row_pos,col_pos]
        z_2 <- (act_val_2 - mean(res_2))/sd(res_2)
        z_val_mat_2[row_pos,col_pos] <- z_2
    }
}

#order rows of z_val_mat
#z_val_mat_0 <- z_val_mat_0[c(5,4,3,2,1),]
#z_val_mat_1 <- z_val_mat_1[c(5,4,3,2,1),]
#z_val_mat_2 <- z_val_mat_2[c(5,4,3,2,1),]

#save Z val matrix
#write.csv(z_val_mat,"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Deviation_from_ancestors_Z_val_matrix_10k_sims.csv")
write.csv(z_val_mat_0,"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Shaved_stem_Deviation_from_ancestors_where_0_ancestral_expression_keeping_ancestral_exp_constant_in_sim_Z_val_matrix_10k_sims.csv")
write.csv(z_val_mat_1,"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Shaved_stem_Deviation_from_ancestors_where_2_ancestral_expression_keeping_ancestral_exp_constant_in_sim_Z_val_matrix_10k_sims.csv")
write.csv(z_val_mat_2,"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Shaved_stem_Deviation_from_ancestors_where_10_ancestral_expression_keeping_ancestral_exp_constant_in_sim_Z_val_matrix_10k_sims.csv")

#read Z val matrix
#z_val_mat <- read.csv("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Deviation_from_ancestors_keeping_ancestral_exp_constant_in_sim_Z_val_matrix_10k_sims.csv",row.names = 1)
#colnames(z_val_mat) <- c(-2,-1,0,1,2)
#head(z_val_mat)

Z_upper_tri_0 <- get_upper_tri(z_val_mat_0)
Z_upper_tri_0[Z_upper_tri_0 == -Inf] <- NA
Z_upper_tri_1 <- get_upper_tri(z_val_mat_1)
Z_upper_tri_1[Z_upper_tri_1 == -Inf] <- NA
Z_upper_tri_2 <- get_upper_tri(z_val_mat_2)
Z_upper_tri_2[Z_upper_tri_2 == -Inf] <- NA



col_fun = colorRamp2(c(-20,0,20),colors =  c("#0008ff", "white", "#e40a0a"))
HM_Z_0 <- Heatmap(as.matrix(Z_upper_tri_0), name = "Z_score", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "gray",
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(c(-2,-1,0,1,2)),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "Duplicate_2_deviation",column_title = "Duplicate_1_deviation",
              cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", Z_upper_tri_0[i, j]), x, y, gp = gpar(fontsize = 10))})
HM_Z_1 <- Heatmap(as.matrix(Z_upper_tri_1), name = "Z_score", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "gray",
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(c(-2,-1,0,1,2)),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "Duplicate_2_deviation",column_title = "Duplicate_1_deviation",
              cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", Z_upper_tri_1[i, j]), x, y, gp = gpar(fontsize = 10))})
HM_Z_2 <- Heatmap(as.matrix(Z_upper_tri_2), name = "Z_score", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "gray",
                row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(c(-2,-1,0,1,2)),
                row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
                row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
                row_title = "Duplicate_2_deviation",column_title = "Duplicate_1_deviation",
                cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", Z_upper_tri_2[i, j]), x, y, gp = gpar(fontsize = 10))})

HM_list <- HM_Z_0 + HM_Z_1 + HM_Z_2

pdf("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Figures/Shaved_stem_Deviation_from_ancestors_where_0_ancestral_expression_grouped_keeping_ancestral_exp_constant_in_sim_Z_val_matrix_10k_sims.pdf",height = 10,width = 30)
HM_list
dev.off()  
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#2A. Looking at the deviation between the ancestral and extant states for all the syntanic homelogs
single branches
#################################
#####Redefine TPM thresholds#####
#################################
high_exp_thresh <- 10

df_TPM <- TPM
df_TPM[df_TPM >= high_exp_thresh] = high_exp_thresh
df_TPM[df_TPM >= 2 & df_TPM < high_exp_thresh ] = 2
df_TPM[df_TPM < 2] = 0

#for ZZtrichome count the presentage of gene that are expressed at high level
#high_exp_genes <- nrow(df_TPM[df_TPM[,c("ZZTrichome")] >= high_exp_thresh,"ZZTrichome",drop =F])/nrow(df_TPM)
#####################################################################################################################################################################
#change paralog dataframe into single gene per row
#Solyc01g110440.5.1
paralogs_df <- paralogs_with_out %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
head(paralogs_df)
paralogs_df[paralogs_df[,1]== 'Solyc02g038817.1.1',]
paralogs_df[paralogs_df[,2]== 2117 | paralogs_df[,2]== 1261,]
paralogs_df[,1] <- gsub("\\.[^.]*$","",paralogs_df[,1])
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)

#chenge ancestral_vals dataframe into single gene per row
ancestral_vals_df <- ancestral_vals %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
ancestral_vals_df <- as.data.frame(ancestral_vals_df, row.names = ancestral_vals_df$in_group)
ancestral_vals_df <- ancestral_vals_df[,c(2,1,3)]
#ancestral_vals_df[,1] <- gsub("\\.[^.]*$","",ancestral_vals_df[,1])
head(ancestral_vals_df)
#cataegorize the expression level in the ancestral state
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= high_exp_thresh] = high_exp_thresh
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= 2 & ancestral_vals_df$ancestral_state < high_exp_thresh ] = 2
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state < 2] = 0
#ancestral_vals_df[,2] <- gsub("\\.[^.]*$","",ancestral_vals_df[,2])
head(ancestral_vals_df)
#ancestral_vals_df[ancestral_vals_df[,3]== high_exp_thresh,]
#####################################################################################################################################################################
#get extarnt state TPM values for in_group genes in ancestral_vals_df in a new column
colnames(df_TPM)   
tissue <- "ZZTrichome"
#rownames(ancestral_vals_df) <- ancestral_vals_df$in_group
#ancestral_vals_df[ancestral_vals_df[,1] ==  'Solyc02g038817.1.1',]
#ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",rownames(ancestral_vals_df)),tissue]
ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",ancestral_vals_df$in_group),tissue]
head(ancestral_vals_df)
df_TPM['Solyc01g109320.4',tissue]

ancestral_vals_diff_df <- ancestral_vals_df 
ancestral_vals_diff_df$ancestral_state[ancestral_vals_diff_df$ancestral_state == 2] = 1
ancestral_vals_diff_df$extant_state[ancestral_vals_diff_df$extant_state == 2] = 1
ancestral_vals_diff_df$ancestral_state[ancestral_vals_diff_df$ancestral_state == high_exp_thresh] = 2
ancestral_vals_diff_df$extant_state[ancestral_vals_diff_df$extant_state == high_exp_thresh] = 2
#get the deviation between the ancestral and extant states
#ancestral_vals_diff_df$deviation <- ancestral_vals_diff_df$extant_state - ancestral_vals_diff_df$ancestral_state
head(ancestral_vals_diff_df)
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#calcullate ancestral vs extant expression quartiles
partition_table <- data.frame(matrix(NA,nrow = 3,ncol = 3))
#make colnames and rownames for the partition_table to be between -2 and 2
colnames(partition_table) <- c(0,1,2)
rownames(partition_table) <- c(0,1,2)
partition_table[is.na(partition_table)] <- 0

#go through each row of the ancestral_vals_diff_df and fill partition_table, where rows are the ancestral state and columns are the extant state
for (i in 1:nrow(ancestral_vals_diff_df)) {
    #i <- 1
    row <- ancestral_vals_diff_df[i,]
    partition_table[as.character(row$ancestral_state),as.character(row$extant_state)] <- partition_table[as.character(row$ancestral_state),as.character(row$extant_state)] + 1
}

#get presentage retention
1753/sum(partition_table[1,1],partition_table[1,2],partition_table[1,3])

2210/sum(partition_table[3,1],partition_table[3,2],partition_table[3,3])

428/sum(partition_table[2,1],partition_table[2,2],partition_table[2,3])
231/sum(partition_table[2,1],partition_table[2,2],partition_table[2,3])
#####################################################################################################################################################################
#simulation of the deviation between the ancestral and extant states for the syntanic homelogs that constain two genes
#REPLACING BELOW FUNCTION WITH LAPPY
#sim_df_list <- list()
#iter <- 10000
#for(i in 1:iter){
#    sim_df <- data.frame(matrix(NA,nrow = 3,ncol = 3))
#    colnames(sim_df) <- c(0,1,2)
#    rownames(sim_df) <- c(0,1,2)
#    sim_df[is.na(sim_df)] <- 0
#    for (j in 1:nrow(ancestral_vals_diff_df)) {
#        #i <- 1
#        row <- ancestral_vals_diff_df[j,]
#        #replace the extant state values with random values
#        row$extant_state <- sample(c(0,1,2),1,replace = TRUE)
        #replace the ancestral state values with random values
        #row$ancestral_state <- sample(c(0,1,2),1,replace = TRUE)
        #new column to store the deviation between the ancestral and extant states
        #row$deviation <- row$extant_state - row$ancestral_state
        #count the frequency of each deviation
 #       sim_df[as.character(row$ancestral_state),as.character(row$extant_state)] <- sim_df[as.character(row$ancestral_state),as.character(row$extant_state)] + 1
 #   }
 #   sim_df_list[[i]] <- sim_df
#}


sim_df_list <- list()
iter <- 10000
sim_df_list <- lapply(1:iter, function(x) {
  sim_df <- matrix(0, nrow = 3, ncol = 3, dimnames = list(c(0, 1, 2), c(0, 1, 2)))
  ancestral_vals_diff_df$extant_state <- sample(c(0, 1, 2), nrow(ancestral_vals_diff_df), replace = TRUE)
  
  for (j in 1:nrow(ancestral_vals_diff_df)) {
    row <- ancestral_vals_diff_df[j, ]
    sim_df[as.character(row$ancestral_state), as.character(row$extant_state)] <- sim_df[as.character(row$ancestral_state), as.character(row$extant_state)] + 1
  }

  sim_df_list[[x]] <- sim_df
})


length(sim_df_list)
z_val_mat <- data.frame(matrix(NA,nrow = 3,ncol = 3))
colnames(z_val_mat) <- c(0,1,2)
rownames(z_val_mat) <- c(0,1,2)
z_val_mat[is.na(z_val_mat)] <- 0
for(col_pos in 1:3){
    for(row_pos in 1:3){
        res <- sapply(sim_df_list, function(x) x[row_pos,col_pos])
        act_val <- partition_table[row_pos,col_pos]
        z <- (act_val - mean(res))/sd(res)
        z_val_mat[row_pos,col_pos] <- z
    }
}
#save Z val matrix
write.csv(z_val_mat,"/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Shaved_stem_Ans_vs_ext_comparison_extant_randomly_Selected_10k_10_TPM_hisg.csv")
#read Z val matrix
z_val_mat <- read.csv("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/data/Shaved_stem_Ans_vs_ext_comparison_extant_randomly_Selected_10k_10_TPM_hisg.csv",row.names = 1)
colnames(z_val_mat) <- c(0,1,2)
col_fun = colorRamp2(c(-20,0,20),colors =  c("#0008ff", "white", "#e40a0a"))
HM_Z <- Heatmap(as.matrix(z_val_mat), name = "Z_score", col =col_fun ,cluster_rows = FALSE,cluster_columns = FALSE, na_col = "gray",
              row_names_gp = gpar(fontsize = 8), show_row_names = TRUE, row_names_side = "left", row_order = as.character(c(0,1,2)),
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
              row_title_rot = 90, column_title_rot = 0, column_title_side = "bottom", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10),
              row_title = "Ancestral_expression_level",column_title = "Extant_expression_level",
              cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", z_val_mat[i, j]), x, y, gp = gpar(fontsize = 10))})
pdf("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Figures/Shaved_stem_Ans_vs_ext_comparison_extant_randomly_Selected_10k_sims_high_exp_INITIAL_THRESHOULD.pdf",height = 10,width = 10)
HM_Z
dev.off() 


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#2B. For the groups that has two paralogs get the number of combinations of expressed-nonexpressed-highly expressed gene combinations in the extant state
high_exp_thresh <- 10

df_TPM <- TPM
df_TPM[df_TPM >= high_exp_thresh] = high_exp_thresh
df_TPM[df_TPM >= 2 & df_TPM < high_exp_thresh ] = 2
df_TPM[df_TPM < 2] = 0

#for ZZtrichome count the presentage of gene that are expressed at high level
#high_exp_genes <- nrow(df_TPM[df_TPM[,c("ZZTrichome")] >= high_exp_thresh,"ZZTrichome",drop =F])/nrow(df_TPM)
#####################################################################################################################################################################
#change paralog dataframe into single gene per row
#Solyc01g110440.5.1
paralogs_df <- paralogs_with_out %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)
head(paralogs_df)
paralogs_df[paralogs_df[,1]== 'Solyc02g038817.1.1',]
paralogs_df[paralogs_df[,2]== 2117 | paralogs_df[,2]== 1261,]
paralogs_df[,1] <- gsub("\\.[^.]*$","",paralogs_df[,1])
paralogs_df <- as.data.frame(paralogs_df, row.names = paralogs_df$in_group)

#chenge ancestral_vals dataframe into single gene per row
ancestral_vals_df <- ancestral_vals %>% 
  mutate(in_group = strsplit(as.character(in_group), ",")) %>% 
  unnest(in_group)
ancestral_vals_df <- as.data.frame(ancestral_vals_df, row.names = ancestral_vals_df$in_group)
ancestral_vals_df <- ancestral_vals_df[,c(2,1,3)]
#ancestral_vals_df[,1] <- gsub("\\.[^.]*$","",ancestral_vals_df[,1])
head(ancestral_vals_df)
#cataegorize the expression level in the ancestral state
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= high_exp_thresh] = high_exp_thresh
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state >= 2 & ancestral_vals_df$ancestral_state < high_exp_thresh ] = 2
ancestral_vals_df$ancestral_state[ancestral_vals_df$ancestral_state < 2] = 0
#ancestral_vals_df[,2] <- gsub("\\.[^.]*$","",ancestral_vals_df[,2])
head(ancestral_vals_df)
#ancestral_vals_df[ancestral_vals_df[,3]== high_exp_thresh,]
#####################################################################################################################################################################
#get extarnt state TPM values for in_group genes in ancestral_vals_df in a new column
colnames(df_TPM)   
tissue <- "ZZTrichome"
#rownames(ancestral_vals_df) <- ancestral_vals_df$in_group
#ancestral_vals_df[ancestral_vals_df[,1] ==  'Solyc02g038817.1.1',]
#ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",rownames(ancestral_vals_df)),tissue]
ancestral_vals_df$extant_state <- df_TPM[gsub("\\.[^.]*$","",ancestral_vals_df$in_group),tissue]
head(ancestral_vals_df)
df_TPM['Solyc01g109320.4',tissue]

ancestral_vals_diff_df <- ancestral_vals_df 
ancestral_vals_diff_df$ancestral_state[ancestral_vals_diff_df$ancestral_state == 2] = "E"
ancestral_vals_diff_df$extant_state[ancestral_vals_diff_df$extant_state == 2] = "E"
ancestral_vals_diff_df$ancestral_state[ancestral_vals_diff_df$ancestral_state == high_exp_thresh] = "H"
ancestral_vals_diff_df$extant_state[ancestral_vals_diff_df$extant_state == high_exp_thresh] = "H"
ancestral_vals_diff_df$ancestral_state[ancestral_vals_diff_df$ancestral_state == 0] = "N"
ancestral_vals_diff_df$extant_state[ancestral_vals_diff_df$extant_state == 0] = "N"
#get the deviation between the ancestral and extant states
#ancestral_vals_diff_df$deviation <- ancestral_vals_diff_df$extant_state - ancestral_vals_diff_df$ancestral_state
head(ancestral_vals_diff_df)
ancestral_vals_diff_df_sub <- ancestral_vals_diff_df[,c("in_group","group","extant_state")]
cat_list <- list()
for (group in unique(ancestral_vals_diff_df_sub$group)) {
    #group <- 210
    sub_df <- ancestral_vals_diff_df_sub[ancestral_vals_diff_df_sub$group == group,]
    
    if (nrow(sub_df) == 2) {
        # combine the extant state of the two genes into one word AND APPEND THAT TO THE LIST
        cat_list[[group]] <- paste(sort(sub_df$extant_state),collapse = "")
        #

    }
}
#DRAW A BAR PLOT OF THE COMBINATIONS
possible_cats <- c("NN","EN","HN","EE","EH","HH")
cat_list_df <- data.frame(table(unlist(cat_list)))
#if var on does not have a possible_cats for a category add 0 to the frequency
for (cat in possible_cats) {
    if (!(cat %in% cat_list_df$Var1)) {
        cat_list_df <- rbind(cat_list_df,data.frame(Var1 = cat,Freq = 0))
    }
}
colnames(cat_list_df) <- c("Combination","Frequency")
#order the dataframe bsed on "NN","HH","EE","HN","EN","EH"
cat_list_df$Combination <- factor(cat_list_df$Combination, levels = c("NN","HH","EE","HN","EN","EH"))

#draw a barplot of the frequency of the combinations in ggplot
pdf("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Ancestral_sates/Figures/At_90th_precentile_extant_combinations.pdf",height = 5,width = 5)
ggplot(cat_list_df, aes(x = Combination, y = Frequency)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Combination") + ylab("Frequency") + ggtitle("25th precentile") +ylim(0,1000)
#chnage y lim

dev.off()
#make