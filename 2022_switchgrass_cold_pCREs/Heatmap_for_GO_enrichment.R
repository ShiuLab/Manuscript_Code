#setwd('D:\\Switchgrass_allelic_specific\\04_presentation\\Data_from_Thilanka')
rm(list=ls())

library(ComplexHeatmap)
library(circlize)

setwd("~/Documents/02_evolution_of_switchgrass_pCREs/RNA_seq/DE/GO_enrichment/Up_regulation")
d_30 <- read.table("GO_D_30.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h_1 <- read.table("GO_D_1.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h_3 <- read.table("GO_D_3.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h_6 <- read.table("GO_D_6.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h_16 <- read.table("GO_D_16.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
h_24 <- read.table("GO_D_24.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')


all_go <- unique(c(d_30[,1],h_1[,1],h_3[,1],h_6[,1],h_16[,1],h_24[,1]))
res <- matrix(data = NA, nrow=length(all_go),ncol=6)
rownames(res) <- all_go
colnames(res) <- c('30_min','1_hr','3_hr','6_hr','16_hr','24_hr')
for(GO in all_go){
	if(length(d_30[d_30[,1]==GO,2]) > 0) res[GO,'30_min'] <- d_30[d_30[,1]==GO,3]
	if(length(h_1[h_1[,1]==GO,2]) > 0) res[GO,'1_hr'] <- h_1[h_1[,1]==GO,3]
	if(length(h_3[h_3[,1]==GO,2]) > 0) res[GO,'3_hr'] <- h_3[h_3[,1]==GO,3]
	if(length(h_6[h_6[,1]==GO,2]) > 0) res[GO,'6_hr'] <- h_6[h_6[,1]==GO,3]
	if(length(h_16[h_16[,1]==GO,2]) > 0) res[GO,'16_hr'] <- h_16[h_16[,1]==GO,3]
	if(length(h_24[h_24[,1]==GO,2]) > 0) res[GO,'24_hr'] <- h_24[h_24[,1]==GO,3]
	}

#res <- -log10(res)
res[is.na(res)] <- 0
res <- res[nrow(res):1,]
res[res > 10] <- 10
res<-res[order(res[,ncol(res)],decreasing=F),]
#library(pvclust)
#library(gplots)

col_fun = colorRamp2(c(-10,0, 10), c("blue","white", "deeppink2"))
HM <-Heatmap(res, name = "-log10(q)", col = col_fun, row_dend_reorder = T, cluster_columns = F)
pdf('Enriched_GO_down-regulation_clusters_with_go.pdf',height = 10,width = 10)
#lmat <- rbind(4:3,2:1)
#lwid = c(0.8,4)
#lhei = c(0.8,4)
#h <- heatmap.2(res,trace="none",col = colorRampPalette(c("white","red"))(8),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.1,Rowv = FALSE,Colv=FALSE,lmat=lmat,lwid=lwid,lhei=lhei)
HM 
dev.off()

HM <-Heatmap(res, name = "-log10(q)", col = col_fun, row_dend_reorder = T)

HM <- Heatmap(as.matrix(FC_sig), name = "logFC", col = col_fun, row_km = 16,
              row_names_gp = gpar(fontsize = 2), show_row_names = FALSE, 
              row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, column_split = dff$group, 
              row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10))
#####################################################################################################################

remove <- c("metabolic process","mitotic chromosome condensation","DNA replication","activation of protein kinase activity",
            "cytokinin metabolic process","ecognition of pollen","amino acid transmembrane transport",
            "RNA processing","RNA aminoacylation for protein translation","intracellular protein transport",
            "translation")

res_sub <- res[!(row.names(res) %in% remove),]

HM <-Heatmap(res_sub, name = "-log10(q)", col = col_fun, row_dend_reorder = T, cluster_columns = F)
pdf('Sub_set_Enriched_GO_down-regulation_clusters_with_go.pdf',height = 10,width = 10)
#lmat <- rbind(4:3,2:1)
#lwid = c(0.8,4)
#lhei = c(0.8,4)
#h <- heatmap.2(res,trace="none",col = colorRampPalette(c("white","red"))(8),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.1,Rowv = FALSE,Colv=FALSE,lmat=lmat,lwid=lwid,lhei=lhei)
HM 
dev.off()


