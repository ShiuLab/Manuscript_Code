####################################################################################
#This script look at the kmer similarities                                        s#
#                                                                                 #
##################################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
####################################################################################
library(cluster)
library("factoextra")
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

####################################################################################
#Data                                                                             #
##################################################################################

min_30 <- read.table("top_features_250_30_min.txt", header = F)
hr_1 <- read.table("top_features_250_1_hr.txt", header = F)
hr_3 <- read.table("top_features_110_3_hr.txt", header = F)
hr_6 <- read.table("top_features_350_6_hr.txt", header = F)
hr_16 <- read.table("top_features_490_16_hr.txt", header = F)
hr_24 <- read.table("top_features_310_24.txt", header = F)

general_kmers <- read.table("general_kmers.txt", header = F)
min_30_specific <- read.table("30_min_specific_kmers.txt", header = F)
hr_1_specifi <-  read.table("1_hr_specific_kmers.txt", header = F)
hr_3_specific <- read.table("3_hr_specific_kmers.txt", header = F)
hr_6_specific <- read.table("6_hr_specific_kmers.txt", header = F)
hr_16_specific <- read.table("16_hr_specific_kmers.txt", header = F)
hr_24_specific <- read.table("24_hr_specific_kmers.txt", header = F)
###########################################################################################

kmer_all <- unique(c(min_30$V1,hr_1$V1,hr_3$V1,hr_6$V1,hr_16$V1,hr_24$V1))
kmer_gen_and_time_specifc <- unique(c(general_kmers$V1,min_30_specific$V1,hr_1_specifi$V1,hr_3_specific$V1,hr_6_specific$V1,hr_16_specific$V1,hr_24_specific$V1))
write.table(as.data.frame(kmer_all),"all_unique_kmaers.txt",quote = F,col.names = F, row.names = F)
write.table(as.data.frame(kmer_gen_and_time_specifc),"general_time_specifc_unique_kmaers.txt",quote = F,col.names = F, row.names = F)
kmer_all_df <- as.data.frame(kmer_all)

file="all_unique_kmaers_no_RC.txt.tamo.dm"
file = "general_time_specifc_unique_kmaers_no_RC.txt.tamo.dm"
kmers <- kmer_all_df
kmers <- as.data.frame(kmer_gen_and_time_specifc)
inp = read.table(file,header=F,sep="\t",na.strings="-",fill=TRUE)
colnames(inp) = kmers$kmer_gen_and_time_specifc
#rownames(inp) = kmers$V1
m = data.matrix(inp)
print(dim(m))
d = as.dist(t(m),upper=FALSE)
hc1 <- agnes(d,diss = TRUE,method="average")
par(cex.lab=0.8, cex.axis=0.1, cex.main=0.7)
plot(hc1, which.plots=2, hang=-1)



pdf('cluster_names.pdf', width = 120)
plot(hc1) 
plot(hc1, which.plots=2, hang=-1)
dev.off()

pdf('gen_time_specifc_cluster_names_cut_tree_0.39.pdf', height=40)
fviz_dend(hc1, cex = 0.6, h = 0.39,rect = T, horiz = T)
dev.off()

'
dend_data <- attr(fviz_dend(hc1, cex = 0.6, h = 0.39,rect = T, horiz = T), "dendrogram")
dend_cuts <- cut(dend_data,  h = 0.23)
fviz_dend(dend_cuts$lower[[2]])
dend_cuts$lower[[2]]
'

kmers$cluster <- cutree(as.hclust(hc1),  h = 0.39)

kmer_orded <- kmers[order(kmers[,2]),]
kmer_orded[kmer_orded[,1]=="GTACCGG.CCGGTAC",]
kmer_orded[kmer_orded[,2]=="39",]
rownames(kmer_orded) <- kmer_orded$kmer_gen_and_time_specifc

kmer_orded$group <- "other"
kmer_orded[general_kmers$V1,3] <- "general_kmers"
kmer_orded[hr_1_specifi$V1,3] <- "1_hr_specific"
kmer_orded[min_30_specific$V1,3] <- "30_min_specific"
kmer_orded[hr_3_specific$V1,3] <- "3_hr_specific"
kmer_orded[hr_6_specific$V1,3] <- "6_hr_specific"
kmer_orded[hr_16_specific$V1,3] <- "16_hr_specific"
kmer_orded[hr_24_specific$V1,3] <- "24_hr_specific"
#kmer_orded[kmer_orded[,3] == "other",]
length(intersect(rownames(kmer_orded[kmer_orded[,3]=="general_kmers",]), general_kmers$V1))
write.table(kmer_orded,"0.39_general_time_Speciifc_kmers_and_similarity_clusters.txt",sep = "\t", quote = F, row.names = F)


kmer_orded <- read.table("0.39_general_time_Speciifc_kmers_and_similarity_clusters.txt",header = T)
kmer_orded[kmer_orded[,3]=="general_kmers",]




kmer_types <-  kmer_orded[,3,drop =F]

colors = structure(c("gray","lightskyblue","lightgoldenrod1", "hotpink1","darkorange","goldenrod","darkviolet","maroon"), names = c("other", "general_kmers", "1_hr_specific","30_min_specific","3_hr_specific","6_hr_specific","16_hr_specific","24_hr_specific"))
kmer_hm <- Heatmap(as.matrix(kmer_types), name = "kmer_type", col = colors, cluster_columns = F,row_dend_reorder = T,
                          row_names_gp = gpar(fontsize = 6), show_row_names = T, row_split = kmer_orded$cluster,
                          row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
                          row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 10))


pdf('cluster_names_cut_tree_hm.pdf', height=120)
kmer_hm
dev.off()




'''
kmers$cluster
order(kmers) <- hc1$order.lab
kmers <- kmers[match(hc1$order.lab,kmers$V1),]

write.csv(kmers,"kmer_cluster_list.txt",sep = "\t")
'''
###########
#making a summery fig
propotion_df <- c()
kmer_groups <-   unique(kmer_orded$group)
#kmer_groups <-  kmer_groups[-2] # only in orther is present
for (cluster in unique(kmer_orded$cluster)) {
  #cluster = 1
  sub_df <- kmer_orded[kmer_orded[,2] ==cluster,]
  for (kmer_type in kmer_groups) {
    if (kmer_type %in% unique(sub_df$group)) {
      #print(kmer_type)
      sub_df_kmer_type <- sub_df[sub_df[,3]==kmer_type,,drop =F]
      prop = (nrow(sub_df_kmer_type)/nrow(sub_df))*100
      propotion_df <- rbind(propotion_df,c(kmer_type,cluster,prop ))

    } else {
      propotion_df <- rbind(propotion_df,c(kmer_type,cluster,0 ))
    }
  }
  
}
 
colnames(propotion_df) <- c("kmer_type", "cluster", "propotion_in_cluster")

propotion_df <- as.data.frame(propotion_df)




#making a boc plot

x_lab <- factor(propotion_df$kmer_type, levels= c("general_kmers", "30_min_specific","1_hr_specific" ,"3_hr_specific" ,"6_hr_specific", "16_hr_specific","24_hr_specific"))
p <- ggplot(propotion_df, aes(x= x_lab, y=as.numeric(propotion_in_cluster))) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(color="black", size=0.4, alpha=0.9) +theme_classic()+
  labs(x="kmer type", y = "% of a kmer present in a cluster") +
  theme(axis.text.x = element_text(angle = 90))
p <- ggplot(propotion_df, aes(x= x_lab, y=as.numeric(propotion_in_cluster))) + 
  geom_boxplot() 

pdf("kmer_type_propotions_in_clusters.pdf", width = 5, height = 5)
p
dev.off()

propotion_df[propotion_df[,1] == "3_hr_specific" & propotion_df[,3] != 0,]
p$data[p$data[,1] == "3_hr_specific" & p$data[,3] != 0,]

#making a heatmaps

sub_df_prop <- propotion_df[propotion_df[,2] ==1,]
tras_df <-t(sub_df_prop)
colnames(tras_df) <- tras_df[1,]
df_df <- as.data.frame(t(tras_df[3,]))
rownames(df_df) <- paste("1","cluster",sep = "_")

kmer_prop_mat <- c()
for (cluster in unique(propotion_df$cluster)) {
  sub_df_prop <- propotion_df[propotion_df[,2] ==cluster,]
  tras_df <-t(sub_df_prop)
  colnames(tras_df) <- tras_df[1,]
  df_df <- as.data.frame(t(tras_df[3,]))
  rownames(df_df) <- paste(cluster,"cluster",sep = "_")
  kmer_prop_mat <- rbind(kmer_prop_mat,df_df)
  
}
kmer_prop_mat <- kmer_prop_mat[,c(1,4,2,7,3,5,6)]
kmer_prop_mat <- as.matrix(kmer_prop_mat)
class(kmer_prop_mat) <- "numeric"
col_fun = colorRamp2(c( 0, 100), c("white", "darkorchid1"))
kmer_hm <- Heatmap(kmer_prop_mat, name = "presentage_similarity", col = col_fun, cluster_columns = F,cluster_rows = F,
                   row_names_gp = gpar(fontsize = 6), show_row_names = T,
                   row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
                   row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 10),column_names_side = "top")
pdf("0.39_thresh_general_time_Specific_kmers_cluster_similarities.pdf", height = 20,width = 10)
kmer_hm
dev.off()
row_order(kmer_hm)


