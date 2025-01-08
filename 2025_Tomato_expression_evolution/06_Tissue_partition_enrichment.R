#This is a script to look at is syntanic paralogs that who a partition uniually tend to show a coupled expression for between two tissue types
#####################################################################################################################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(7)
##########################Libraries##############################
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
###########################Data##########################################################################################################################################
TPM <- read.csv("TPM_table_non_T_tissues.txt",row.names=1, header = TRUE, sep = "\t")
colnames(TPM) 
#remove Mature.anthers
TPM <- TPM[,-which(colnames(TPM) == "Mature.anthers")]
dff_TPM <- read.table('Tissue_cluster_groups_TPM.txt',head=T,sep='\t',stringsAsFactors=F)
dff_nont_T <- dff_TPM[match(colnames(TPM) , dff_TPM[,1]),]
#change "Seedling" in group to "Root
dff_nont_T$group[dff_nont_T$group == "Seedling"] <- "Root"
#order rows based of group column
order_groups <- factor(dff_nont_T$group, levels = c("Fruit","Seed","Root","Callus","Floral","Stem","Trichome"))

dff_nont_T <- dff_nont_T[order(order_groups),]
nrow(dff_nont_T)
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

#####################################################################################################################################################################

#table to store log odds ratio
#make empty table using colnames of paralogs_TPM
enrichment_tables_logg_odds <- matrix(NA, nrow = ncol(paralogs_TPM), ncol = ncol(paralogs_TPM))
colnames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
rownames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
#change column order of enrichment_tables to match with the order of dff_nont_T
enrichment_tables_logg_odds <- enrichment_tables_logg_odds[dff_nont_T$tissue,dff_nont_T$tissue]

#table to store p values
enrichment_tables_pvals <- enrichment_tables_logg_odds
#table to store the counts of the paralogs
paralogs_counts <- data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(paralogs_counts) <- c("tissue_1","tissue_2","gene_1_both_tussues","gene_2_both_tussues","gene_1_T1_gene_2_T2","gene_1_T2_gene_2_T1")
for(tissue_1 in unique(colnames(paralogs_TPM))){
  for(tissue_2 in unique(colnames(paralogs_TPM))){
    if (tissue_1 != tissue_2) {
        tissue_1 <- "ZZTrichome"
        tissue_2 <- "Root"
        gene_1_both_tussues <- 0
        gene_2_both_tussues <- 0
        gene_1_T1_gene_2_T2 <- 0
        gene_1_T2_gene_2_T1 <- 0

      #dataframe with two tissue types
      paralogs_TPM_sub_df <- paralogs_TPM[,c(tissue_1,tissue_2)]
      #look at each pair of paralogs
      for(group in unique(paralogs_df$Group)){
        #group <- 204
        sub_partion_df <- paralogs_TPM_sub_df[rownames(paralogs_df[paralogs_df[,2] == group,]),]
        #CHENGING THE EXPRESSED TO NON EXPRESSED TO USING DIFFERENT TPM THRESHOULDS
        #sub_partion_df[sub_partion_df < 10] = 0
        #sub_partion_df[sub_partion_df >= 10] = 1
        sub_partion_df[sub_partion_df >= 2] = 1
        #get every pair of paralogs in the group
        row_combinations <- combn(rownames(sub_partion_df),2, simplify = TRUE)
        for(i in 1:ncol(row_combinations)){
          #i = 1
          #get the partition score for the pair
          sub_partion_df_with_partitions <- sub_partion_df[row_combinations[,i],]
          if(all(colSums(sub_partion_df_with_partitions) == 1)){
            #print the sub_partion_df_with_partitions
            #print(sub_partion_df_with_partitions)
            #get fisher exact scores and do the test
            if (rowSums(sub_partion_df_with_partitions)[1] == 2){
              gene_1_both_tussues <- gene_1_both_tussues + 1
          
          }else if (rowSums(sub_partion_df_with_partitions)[2] == 2){
            gene_2_both_tussues <- gene_2_both_tussues + 1
        }else if (sub_partion_df_with_partitions[1,2] == 1){
          gene_1_T2_gene_2_T1 <- gene_1_T2_gene_2_T1 + 1
        }else if (sub_partion_df_with_partitions[2,2] == 1){
          gene_1_T1_gene_2_T2 <- gene_1_T1_gene_2_T2 + 1
        }
          }
        }      
      }
      #add counts to paralogs_counts
      paralogs_counts <- rbind(paralogs_counts, c(tissue_1,tissue_2,gene_1_both_tussues,gene_2_both_tussues,gene_1_T1_gene_2_T2,gene_1_T2_gene_2_T1))
      #make a contingency table using the counts
      contingency_table <- matrix(c(gene_1_both_tussues,gene_1_T1_gene_2_T2,gene_1_T2_gene_2_T1,gene_2_both_tussues), nrow = 2, ncol = 2)
      #if contingency_table contains 0s fill it with 1s
      if (any(contingency_table == 0)){
        contingency_table <- contingency_table + 1
        #contingency_table[contingency_table == 0] <- 1
      }
      #do the fisher exact test
      fisher_exact_test <- fisher.test(contingency_table)
      #get p value
      p_value <- fisher_exact_test$p.value
      #get odds ratio
      log10_odds_ratio <- log2(fisher_exact_test$estimate)
      #if the p value is less than 0.05 fill enrichment_tables for corresponding tissues
      
      enrichment_tables_logg_odds[tissue_1,tissue_2] <- log10_odds_ratio
      enrichment_tables_pvals[tissue_2,tissue_1] <- p_value
      

    }
  }
} 

#########################################
#modify the enrichment_tables
#save the paralogs_counts
colnames(paralogs_counts) <- c("tissue_1","tissue_2","gene_1_both_tussues","gene_2_both_tussues","gene_1_T1_gene_2_T2","gene_1_T2_gene_2_T1")
write.table(paralogs_counts, file = "tissue_enichment/paralogs_counts_for_enrichment_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)

#save the enrichment_tables_final
write.table(enrichment_tables_pvals, file = "tissue_enichment/enrichment_tables_p_values_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)
write.table(enrichment_tables_logg_odds, file = "tissue_enichment/enrichment_tables_logodds_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)

#FILL NA WITH 0s
#pvalue table
enrichment_tables_pvals_mod <- enrichment_tables_pvals
enrichment_tables_pvals_mod[enrichment_tables_pvals_mod == 0] <- min(enrichment_tables_pvals_mod[enrichment_tables_pvals_mod != 0],na.rm = TRUE)
enrichment_tables_pvals_mod <-  -log10(enrichment_tables_pvals_mod)
enrichment_tables_pvals_mod[is.na(enrichment_tables_pvals_mod)] <- 0

#if there are any 0s fill them with min value


enrichment_tables_logg_odds_mod <- enrichment_tables_logg_odds
enrichment_tables_logg_odds_mod[is.na(enrichment_tables_logg_odds_mod)] <- 0

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
#pvalues
upper_tri <- get_upper_tri(enrichment_tables_pvals_mod)
#log odds ratio
lower_tri <- get_lower_tri(enrichment_tables_logg_odds_mod)
 
#rowSums(sub_partion_df_with_partitions)
#####################################################################################################################################################################
#draw a heatmap of the enrichment_tables matrix using melted comact format
# Melt the correlation matrix
#heatmap for the upper triangle
melted_upper_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap

hm_pvals <- ggplot(data = melted_upper_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#eaee12", high = "#048563", mid = "#00ceb2", 
   midpoint = median(c(-log10(0.05),max(upper_tri,na.rm = TRUE))), limit = c(-log10(0.05),max(upper_tri,na.rm = TRUE)), space = "Lab", na.value = "white",
   name="-log10(pvalues)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
 scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/pvals_enrichment_between_tissue_paires_cheking_for_paired_expression_TPM_10_as_Expressed.pdf", width = 10, height = 10)
hm_pvals
dev.off()

#heatmap for the lower triangle
melted_lower_cormat <- melt(lower_tri, na.rm = TRUE)
# Heatmap
hm_logods <- ggplot(data = melted_lower_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high ="deeppink2", mid = "white", 
   midpoint = 0, limit = c(min(lower_tri, na.rm = TRUE),max(lower_tri,na.rm = TRUE)), space = "Lab", 
   name="log2(odds_ratio)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
  scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/odds_ratio_enrichment_between_tissue_paires_cheking_for_paired_expression_TPM_2_as_Expressed.pdf", width = 10, height = 10)
hm_logods
dev.off()

pdf("tissue_enichment/odds_ratio_and_pvals_enrichment_between_tissue_paires_cheking_for_paired_expression_TPM_2_as_Expressed.pdf", width = 20, height = 10)
grid.arrange(hm_logods,hm_pvals, ncol = 2, nrow = 1)
dev.off()





#####################################################################################################################################################################
#####################################################################################################################################################################
#DO THE ABOVE ANALYSIS USING NON BROADLY EXPRESSED GENES

#table to store log odds ratio
#make empty table using colnames of paralogs_TPM
enrichment_tables_logg_odds <- matrix(NA, nrow = ncol(paralogs_TPM), ncol = ncol(paralogs_TPM))
colnames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
rownames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
#change column order of enrichment_tables to match with the order of dff_nont_T
enrichment_tables_logg_odds <- enrichment_tables_logg_odds[dff_nont_T$tissue,dff_nont_T$tissue]

#table to store p values
enrichment_tables_pvals <- enrichment_tables_logg_odds
#table to store the counts of the paralogs
paralogs_counts <- data.frame(matrix(NA,nrow = 0,ncol = 6))
colnames(paralogs_counts) <- c("tissue_1","tissue_2","gene_1_both_tussues","gene_2_both_tussues","gene_1_T1_gene_2_T2","gene_1_T2_gene_2_T1")
for(tissue_1 in unique(colnames(paralogs_TPM))){
  for(tissue_2 in unique(colnames(paralogs_TPM))){
    if (tissue_1 != tissue_2) {
        #tissue_1 <- "ZZTrichome"
        #tissue_2 <- "Root"
        gene_1_both_tussues <- 0
        gene_2_both_tussues <- 0
        gene_1_T1_gene_2_T2 <- 0
        gene_1_T2_gene_2_T1 <- 0

      #dataframe with two tissue types
      paralogs_TPM_sub_df <- paralogs_TPM[,c(tissue_1,tissue_2)]
      #look at each pair of paralogs
      for(group in unique(paralogs_df$Group)){
        #group <- 204
        sub_partion_df <- paralogs_TPM_sub_df[rownames(paralogs_df[paralogs_df[,2] == group,]),]
        #CHENGING THE EXPRESSED TO NON EXPRESSED TO USING DIFFERENT TPM THRESHOULDS
        #sub_partion_df[sub_partion_df < 10] = 0
        #sub_partion_df[sub_partion_df >= 10] = 1
        sub_partion_df[sub_partion_df >= 2] = 1
        #get every pair of paralogs in the group
        row_combinations <- combn(rownames(sub_partion_df),2, simplify = TRUE)
        for(i in 1:ncol(row_combinations)){
          #i = 1
          #get the partition score for the pair
          sub_partion_df_with_partitions <- sub_partion_df[row_combinations[,i],]
          sub_partion_df_with_partitions_TPM_across <- paralogs_TPM[rownames(sub_partion_df_with_partitions),]
          sub_partion_df_with_partitions_TPM_across[sub_partion_df_with_partitions_TPM_across >= 2] = 1
          if(!any(rowSums(sub_partion_df_with_partitions_TPM_across) > 2)){

            if(all(colSums(sub_partion_df_with_partitions) == 1)){
              #print the sub_partion_df_with_partitions
              #print(sub_partion_df_with_partitions)
              #get fisher exact scores and do the test
              if (rowSums(sub_partion_df_with_partitions)[1] == 2){
                gene_1_both_tussues <- gene_1_both_tussues + 1
            
              }else if (rowSums(sub_partion_df_with_partitions)[2] == 2){
                gene_2_both_tussues <- gene_2_both_tussues + 1
            }else if (sub_partion_df_with_partitions[1,2] == 1){
              gene_1_T2_gene_2_T1 <- gene_1_T2_gene_2_T1 + 1
            }else if (sub_partion_df_with_partitions[2,2] == 1){
              gene_1_T1_gene_2_T2 <- gene_1_T1_gene_2_T2 + 1
            }
            }
          }
        }      
      }
      #add counts to paralogs_counts
      paralogs_counts <- rbind(paralogs_counts, c(tissue_1,tissue_2,gene_1_both_tussues,gene_2_both_tussues,gene_1_T1_gene_2_T2,gene_1_T2_gene_2_T1))
      #make a contingency table using the counts
      contingency_table <- matrix(c(gene_1_both_tussues,gene_1_T1_gene_2_T2,gene_1_T2_gene_2_T1,gene_2_both_tussues), nrow = 2, ncol = 2)
      #if contingency_table contains 0s fill it with 1s
      if (any(contingency_table == 0)){
        #contingency_table[contingency_table == 0] <- 1
        contingency_table <- contingency_table + 1
      }
      #do the fisher exact test
      fisher_exact_test <- fisher.test(contingency_table)
      #get p value
      p_value <- fisher_exact_test$p.value
      #get odds ratio
      log10_odds_ratio <- log10(fisher_exact_test$estimate)
      #if the p value is less than 0.05 fill enrichment_tables for corresponding tissues
      if (p_value < 0.05) {
        enrichment_tables_logg_odds[tissue_1,tissue_2] <- log10_odds_ratio
        enrichment_tables_pvals[tissue_2,tissue_1] <- p_value
      }

    }
  }
} 

#########################################
#save the paralogs_counts
colnames(paralogs_counts) <- c("tissue_1","tissue_2","gene_1_both_tussues","gene_2_both_tussues","gene_1_T1_gene_2_T2","gene_1_T2_gene_2_T1")
write.table(paralogs_counts, file = "tissue_enichment/paralogs_counts_for_enrichment_non_broadly_expressed_2_tissues_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)

#save the enrichment_tables_final
write.table(enrichment_tables_pvals, file = "tissue_enichment/enrichment_tables_p_values_non_broadly_expressed_2_tissues_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)
write.table(enrichment_tables_logg_odds, file = "tissue_enichment/enrichment_tables_logg_odds_non_broadly_expressed_2_tissues_TPM_2_as_Expressed.txt", sep = "\t", quote = FALSE)

#FILL NA WITH 0s
#pvalue table
enrichment_tables_pvals_mod <- enrichment_tables_pvals
enrichment_tables_pvals_mod <-  -log10(enrichment_tables_pvals_mod)
enrichment_tables_pvals_mod[is.na(enrichment_tables_pvals_mod)] <- 0

enrichment_tables_logg_odds_mod <- enrichment_tables_logg_odds
enrichment_tables_logg_odds_mod[is.na(enrichment_tables_logg_odds_mod)] <- 0

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
#pvalues
upper_tri <- get_upper_tri(enrichment_tables_pvals_mod)
#log odds ratio
lower_tri <- get_lower_tri(enrichment_tables_logg_odds_mod)
 
#rowSums(sub_partion_df_with_partitions)
#####################################################################################################################################################################
#draw a heatmap of the enrichment_tables matrix using melted comact format
# Melt the correlation matrix
#heatmap for the upper triangle
melted_upper_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap

hm <- ggplot(data = melted_upper_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#a8a576", high = "#ffa600", mid = "#decf4a", 
   midpoint = median(c(-log10(0.05),max(upper_tri,na.rm = TRUE))), limit = c(-log10(0.05),max(upper_tri,na.rm = TRUE)), space = "Lab", na.value = "white",
   name="-log10(pvalues)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
 scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/pvals_enrichment_between_tissue_paires_cheking_for_paired_expression_excluding_broadlyexpressed_5_tissues_TPM_2_as_Expressed.pdf", width = 10, height = 10)
hm
dev.off()

#heatmap for the lower triangle
melted_lower_cormat <- melt(lower_tri, na.rm = TRUE)
# Heatmap
hm <- ggplot(data = melted_lower_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high ="deeppink2", mid = "white", 
   midpoint = 0, limit = c(min(lower_tri, na.rm = TRUE),max(lower_tri,na.rm = TRUE)), space = "Lab", 
   name="log10(odds_ratio)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
  scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/odds_ratio_enrichment_between_tissue_paires_cheking_for_paired_expression_excluding_broadlyexpressed_5_tissues__TPM_2_as_Expressed.pdf", width = 10, height = 10)
hm
dev.off()



#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
# Re run fisher exact test based off using paralogous counts table without having to run the whole analysis again

#read the paralogs_counts
#counts_table_nm <- "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/tissue_enichment/secondary_analysis/paralogs_counts_for_enrichment_non_broadly_expressed_5_tissues_TPM_2_as_Expressed.txt"
#counts_table_nm <- "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/tissue_enichment/secondary_analysis/paralogs_counts_for_enrichment_non_broadly_expressed_10_tissues_TPM_2_as_Expressed.txt"
counts_table_nm <- "/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/tissue_enichment/paralogs_counts_for_enrichment_non_broadly_expressed_2_tissues_TPM_2_as_Expressed.txt"
paralogs_counts <- read.table(counts_table_nm,header = TRUE,sep = "\t")
head(paralogs_counts)

#table to store log odds ratio
#make empty table using colnames of paralogs_TPM
enrichment_tables_logg_odds <- matrix(NA, nrow = ncol(paralogs_TPM), ncol = ncol(paralogs_TPM))
colnames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
rownames(enrichment_tables_logg_odds) <- colnames(paralogs_TPM)
#change column order of enrichment_tables to match with the order of dff_nont_T
enrichment_tables_logg_odds <- enrichment_tables_logg_odds[dff_nont_T$tissue,dff_nont_T$tissue]

#table to store p values
enrichment_tables_pvals <- enrichment_tables_logg_odds

#function to go though each row of the paralogs_counts and do the fisher exact test
for(i in 1:nrow(paralogs_counts)){
  #i = 1
  tissue_1 <- paralogs_counts[i,1]
  tissue_2 <- paralogs_counts[i,2]
  gene_1_both_tussues <- paralogs_counts[i,3]
  gene_2_both_tussues <- paralogs_counts[i,4]
  gene_1_T1_gene_2_T2 <- paralogs_counts[i,5]
  gene_1_T2_gene_2_T1 <- paralogs_counts[i,6]
  #make a contingency table using the counts
  contingency_table <- matrix(c(gene_1_both_tussues,gene_1_T1_gene_2_T2,gene_1_T2_gene_2_T1,gene_2_both_tussues), nrow = 2, ncol = 2)
  #if contingency_table contains 0s fill it with 1s
  if (any(contingency_table == 0)){
    #contingency_table[contingency_table == 0] <- 1
    contingency_table <- contingency_table + 1
  }
  #do the fisher exact test
  fisher_exact_test <- fisher.test(contingency_table)
  #get p value
  p_value <- fisher_exact_test$p.value
  #get odds ratio
  log2_odds_ratio <- log2(fisher_exact_test$estimate)
  #if the p value is less than 0.05 fill enrichment_tables for corresponding tissues
  
  enrichment_tables_logg_odds[tissue_1,tissue_2] <- log2_odds_ratio
  enrichment_tables_pvals[tissue_2,tissue_1] <- p_value
  
}


#save the enrichment_tables_final
write.table(enrichment_tables_pvals, file = "tissue_enichment/enrichment_tables_p_values_non_broadly_expressed_2_tissues_TPM_2_as_Expressed_phrasing_all.txt", sep = "\t", quote = FALSE)
write.table(enrichment_tables_logg_odds, file = "tissue_enichment/enrichment_tables_logg_odds_non_broadly_expressed_2_tissues_TPM_2_as_Expressed_phrasing_all.txt", sep = "\t", quote = FALSE)

#FILL NA WITH 0s
#pvalue table
enrichment_tables_pvals_mod <- enrichment_tables_pvals
enrichment_tables_pvals_mod <-  -log10(enrichment_tables_pvals_mod)
enrichment_tables_pvals_mod[is.na(enrichment_tables_pvals_mod)] <- 0

enrichment_tables_logg_odds_mod <- enrichment_tables_logg_odds
enrichment_tables_logg_odds_mod[is.na(enrichment_tables_logg_odds_mod)] <- 0

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
#pvalues
upper_tri <- get_upper_tri(enrichment_tables_pvals_mod)
#log odds ratio
lower_tri <- get_lower_tri(enrichment_tables_logg_odds_mod)
 
#rowSums(sub_partion_df_with_partitions)
#####################################################################################################################################################################
#draw a heatmap of the enrichment_tables matrix using melted comact format
# Melt the correlation matrix
#heatmap for the upper triangle
melted_upper_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
#colors values
max_val <- max(upper_tri,na.rm = TRUE)
#max_val <- 30

hm_pvals <- ggplot(data = melted_upper_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#eaee12", high = "#048563", mid = "#00ceb2", 
   midpoint = median(c(-log10(0.05),max_val)), limit = c(-log10(0.05),max_val), space = "Lab", na.value = "#e7e7e7",
   name="-log10(pvalues)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
 scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/pvals_enrichment_between_tissue_paires_cheking_for_paired_expression_excluding_broadlyexpressed_2_tissues_TPM_2_as_Expressed_phrasing_all.pdf", width = 10, height = 10)
hm_pvals
dev.off()

#heatmap for the lower triangle
melted_lower_cormat <- melt(lower_tri, na.rm = TRUE)
# Heatmap
hm_logodds <- ggplot(data = melted_lower_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#000be0", high ="deeppink2", mid = "white", 
   midpoint = 0, limit = c(min(lower_tri, na.rm = TRUE),max(lower_tri,na.rm = TRUE)), space = "Lab", 
   name="log2(odds_ratio)") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() +
  scale_y_discrete(limits = rev, guide = guide_axis(angle = 0))

#draw the heatmap but fliped

#save hm to pdf
pdf("tissue_enichment/odds_ratio_enrichment_between_tissue_paires_cheking_for_paired_expression_excluding_broadlyexpressed_10_tissues__TPM_2_as_Expressed_phrasing_all.pdf", width = 10, height = 10)
hm_logodds
dev.off()

#print both heatmaps ontop of each other

pdf("tissue_enichment/odds_ratio_and_pvals_enrichment_between_tissue_paires_cheking_for_paired_expression_excluding_broadlyexpressed_2_tissues_TPM_2_as_Expressed_phrasing_all.pdf", width = 20, height = 10)
grid.arrange(hm_logodds,hm_pvals, ncol = 2, nrow = 1)
dev.off()
