#/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Enrichment analysis/all_gene_background/Up_regulation sets/GO_and_Pathway_enrichment_analysis.r
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(gplots)
library(ComplexHeatmap)
library(circlize)
library('GOstats')
library('GSEABase')
######################################################################################################################
GO <- read.csv('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Motif discovary/Choosing positive and negative sets on FDR/without_anther/Enrichments_for_gene_sets/GO/GO_extracted_fromITAG4.0_gene_models_formatted.txt',head=F,sep='\t')
head(GO)
tail(GO)
null_num <- length(unique(GO[,1]))
head(GO)
#functions#
###########

GO_pipeline <- function(Genes_file,short_name_to_save){
	pos <- read.csv(Genes_file,head=F,sep='\t')
	#pos[,1] <- gsub("Pavir.","",pos[,1])
	#pos[,1] <- gsub(".v5.1","",pos[,1])
	pos_go <- unique(merge(GO,pos,by.x=colnames(GO)[1],by.y=colnames(pos)[1]))
	#pos_num <- nrow(pos)
	pos_num <- length(unique(pos_go[,1]))
	#dataframe with gene that are not in pos
	null_go <- GO[!GO[,1] %in% pos_go[,1],]
	#check if any pos_go genes are present in null go
	#null_go[null_go[,1] %in% pos_go[,1],]
	res <- c()
	go <- unique(GO[,2])
	for(i in 1:length(go)){
		Pos <- pos_go[pos_go[,2]==go[i],]
		Null <- null_go[null_go[,2]==go[i],]
		a = nrow(Pos) # genes in cluster and with GO term
		b = pos_num - a # genes in cluster but with no GO term
		cc = nrow(Null) # genes with GO but not in cluster
		d = length(unique(null_go[,1])) - cc# genes with no GO and also not in cluster
		out <- c(go[i],a,b,cc,d)
		if(a+cc != 0) res <- rbind(res,out)
		}
	colnames(res) <- c('GO','Genes_responsive_with_GO','Genes_responsive_with_no_GO','Genes_not_responsive_with_GO','Genes_not_responsive_with_no_GO')
	write.table(res,paste('GO_',short_name_to_save,'.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	enrichment <- c()
	for(i in 1:nrow(res)){
		numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
		exact_test <- fisher.test(numbers, alternative = "two.sided")
        log_odds <-  log10(exact_test$estimate)
        p <- exact_test$p.value
		enrichment <- rbind(enrichment, c(res[i,],log_odds,p))
		}
    #naming last column of enrichment
    colnames(enrichment)[ncol(enrichment)] <- "p_value"
	write.table(enrichment,paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')
	dat <- read.table(paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
    #remove the rows where V7 = 1
    #this is done to remove non informaive p-vals so when p is adjusted the p val distribution does not skew towards 1
    dat <- dat[dat[,7] != 1,]
    #remove any rows where V2 = 0
    #this is to get rid of the inf comming from first class being zero
    dat <- dat[dat[,2] != 0,]
	dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
	dat <- dat[order(dat[,6]),]  
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- cbind(dat,'BP'='')
	dat <- cbind(dat,'CC'='')
	dat <- cbind(dat,'MF'='')
	for(i in 1:nrow(dat)){
		tryCatch( {
			if(!is.null(getGOTerm(dat[i,1])$BP[1])) dat[i,9] <- getGOTerm(dat[i,1])$BP[1]
			if(!is.null(getGOTerm(dat[i,1])$CC[1])) dat[i,10] <- getGOTerm(dat[i,1])$CC[1]
			if(!is.null(getGOTerm(dat[i,1])$MF[1])) dat[i,11] <- getGOTerm(dat[i,1])$MF[1]
			},
			error = function(e) {print(paste("no GO for ",dat[i,1]));NaN},
			finally = {})
		}
	write.table(dat,paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	#subdat <- dat[dat[,8] < 0.05 & dat$BP != '',]
    subdat <- dat[dat[,8] < 0.05,]
    subdat$S <- subdat$BP
    subdat$S[subdat$BP==""] <- subdat$MF[subdat$BP==""]
    subdat$S[subdat$S==""] <- subdat$CC[subdat$S==""]
    #remove any rows where S = ""
    subdat <- subdat[subdat$S != "",]
	if (nrow(subdat) >= 1) {
	write.table(subdat[,c(6,8,12)],paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term_BP.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	df <- subdat[,c('S','V6')]
	rownames(df) <- subdat$S
    df <- df[order(df[,2],decreasing=T),-1,drop=FALSE]
	#change any inf to 1
	df[df == Inf] <- 1
    df <- as.matrix(df)

	col_fun_2 = colorRamp2(c( -1, 0,1), c("#3c00ff","white","#ff0026"))
    kmer_hm_imp <- Heatmap(as.matrix(df), name = "log_odds_ratio", col = col_fun_2,
                           row_names_gp = gpar(fontsize = 6), show_row_names = T,
                           row_gap = unit(3.5, "mm"),column_gap = unit(3.5, "mm"), column_names_gp = gpar(fontsize = 8), border = TRUE, 
                           row_title_rot = 0, column_title_rot = 90, column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 10),column_names_side = "top" )
	pdf(paste("GO_", short_name_to_save,'.fisher.qvalue_GO_term_BP.pdf',sep=''))
    draw(kmer_hm_imp)
    dev.off()
	}else{
  print("No_enrichment")
}
}

GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_all_tissues.txt','both_no_exp')
GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_expressed_in_all_tissues.txt','both_exp')
GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/only_one_paralog_expressed_in_all_tissues.txt','only_one_exp')
GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/one_paralog_expressed_in_one_tissue_while_other_non.txt','one_exp_one_not')
GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_al_46_different_tissues.txt','not_exp_in_all_46_tissues')
GO_pipeline('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_tissues_contain_trichome_but_not_in_23.txt','exp_trichome_tissue')
	
	pos <- read.csv('/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/both_paralogs_not_expressed_in_tissues_contain_trichome_but_not_in_23.txt',head=F,sep='\t')
	short_name_to_save <- 'exp_trichome_tissue'
	#####################################################################################################################
#combine all the results
NDNNNN <- read.table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Go_enrichment/GO_both_exp.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
NNDNNN <- read.table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Go_enrichment/GO_only_one_exp.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')
NNNNDN <- read.table("/Users/thilanka_ranaweera/Documents/prj_01_S_lycopersicum_trichome/Synteny/Expression_fr_paralogs/Go_enrichment/GO_both_no_exp.fisher.qvalue_GO_term_BP.txt",head=T,sep='\t')


all_go <- unique(c(NDNNNN[,3],NNDNNN[,3],NNNNDN[,3]))
res <- matrix(data = NA, nrow=length(all_go),ncol=3)
rownames(res) <- all_go
colnames(res) <- c('NDNNNN','NNDNNN','NNNNDN')
for(GO in all_go){
    #if not GO == ""
    if (GO == "") next
	if(length(NDNNNN[NDNNNN[,3]==GO,1]) > 0) res[GO,'NDNNNN'] <- NDNNNN[NDNNNN[,3]==GO,1]
	if(length(NNDNNN[NNDNNN[,3]==GO,1]) > 0) res[GO,'NNDNNN'] <- NNDNNN[NNDNNN[,3]==GO,1]
	if(length(NNNNDN[NNNNDN[,3]==GO,1]) > 0) res[GO,'NNNNDN'] <- NNNNDN[NNNNDN[,3]==GO,1]
		
	}

#res <- -log10(res)
res[is.na(res)] <- 0
res <- res[nrow(res):1,]
res[res > 10] <- 10
res<-res[order(res[,ncol(res)],decreasing=F),]
colnames(res) <- c("both_exp","only_one_exp","no_exp")
#library(pvclust)
#library(gplots)

col_fun = colorRamp2(c(-2,0, 2), c("#258ef7","white", "#f900e4"))
HM <-Heatmap(res, name = "log(odds_ratio)", col = col_fun, row_dend_reorder = T, cluster_columns = F)
pdf('./Enriched_GO_divergence_types.pdf',height = 25,width = 5)
#lmat <- rbind(4:3,2:1)
#lwid = c(0.8,4)
#lhei = c(0.8,4)
#h <- heatmap.2(res,trace="none",col = colorRampPalette(c("white","red"))(8),dendrogram='none',cex.axis=0.4,notecex=0.4,cexRow=0.1,cexCol=0.1,Rowv = FALSE,Colv=FALSE,lmat=lmat,lwid=lwid,lhei=lhei)
HM 
dev.off()

#for the presentation
col_fun = colorRamp2(c(-2,0, 2), c("blue","white", "deeppink2"))
HM <-Heatmap(res, name = "log(odds_ratio)", col = col_fun, row_dend_reorder = T, cluster_columns = F)
pdf('Enriched_GO_up-regulation_FOR_PRESENTATION.pdf',height = 5,width = 5)
HM 
dev.off()