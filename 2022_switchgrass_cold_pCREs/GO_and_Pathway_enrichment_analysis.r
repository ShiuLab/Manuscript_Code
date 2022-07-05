
rm(list=ls())
setwd("~/Documents/02_evolution_of_switchgrass_pCREs/RNA_seq/DE/GO_enrichment/Up_regulation")
 
Enrichment <- function(k,n,C,G){
	return((k / C) / (n/ G))
	}
library('GOstats')
library('GSEABase')
library(pvclust)
library(gplots)
GO <- read.table('Switchgrass_GO_annotation.txt',head=F,sep='\t')
# null <- rownames(read.table('D:\\Transformer_switchgrass\\02_RNA-seq_results\\Swithchgrass_genes_FC_nochange_FC0.5_FDR0.1_unique.txt',head=T,sep='\t'))
# null <- gsub("Pavir.","",null)
# null <- gsub(".v5.1","",null)
# null_go <- unique(merge(GO,t(t(null)),by.x='V1',by.y='V1'))
# null_num <- length(null)
null_go <- unique(GO)
null_num <- 80278

GO_pipeline <- function(Genes_file,short_name_to_save){
	pos <- read.table(Genes_file,head=F,sep='\t')
	pos[,1] <- gsub("Pavir.","",pos[,1])
	pos[,1] <- gsub(".v5.1","",pos[,1])
	pos_go <- unique(merge(GO,pos,by.x='V1',by.y='V1'))
	pos_num <- nrow(pos)
	res <- c()
	go <- unique(GO[,2])
	for(i in 1:length(go)){
		Pos <- pos_go[pos_go[,2]==go[i],]
		Null <- null_go[null_go[,2]==go[i],]
		a = nrow(Pos) # genes in cluster and with GO term
		b = pos_num - a # genes in cluster but with no GO term
		cc = nrow(Null) # genes with GO but not in cluster
		d = null_num - cc# genes with no GO and also not in cluster
		out <- c(go[i],a,b,cc,d)
		if(a+cc != 0) res <- rbind(res,out)
		}
	colnames(res) <- c('GO','Genes_responsive_with_GO','Genes_responsive_with_no_GO','Genes_not_responsive_with_GO','Genes_not_responsive_with_no_GO')
	write.table(res,paste('GO_',short_name_to_save,'.txt',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	enrichment <- c()
	for(i in 1:nrow(res)){
		numbers <- matrix(as.numeric(res[i,2:5]),nrow = 2)
		p <- fisher.test(numbers, alternative = "two.sided")[[1]][1]
		a = as.numeric(res[i,2])
		b = as.numeric(res[i,3])
		cc = as.numeric(res[i,4])
		d = as.numeric(res[i,5])
		if(Enrichment(a, a+b, a+cc, a+b+cc+d) >= 1) direction = '+' else direction = '-'
		enrichment <- rbind(enrichment, c(res[i,],direction,p))
		}
	write.table(enrichment,paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),row.names=F,col.names=F,quote=F,sep='\t')

	dat <- read.table(paste('GO_',short_name_to_save,'.fisher.pvalue',sep=''),head=F,sep='\t',stringsAsFactors=F)
	dat <- cbind(dat,p.adjust(dat[,7], method = "BH"))
	dat <- dat[order(dat[,8]),]  
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
	
	subdat <- dat[dat[,8] < 0.05 & dat$BP != '',]


	if (nrow(subdat) >= 2) {
	  print(nrow(subdat))
	  subdat <- cbind(subdat,'logP'=0)
	  for(i in 1:nrow(subdat)){
	    if(subdat[i,6] == '-')subdat[i,ncol(subdat)] = log10(subdat[i,8])
	    if(subdat[i,6] == '+')subdat[i,ncol(subdat)] = -log10(subdat[i,8])
	  }
	  subdat[subdat$logP < -10,12] <- -10
	  subdat[subdat$logP > 10,12] <- 10
	  subdat$color <- 10 #added thilanka
	  subdat[1,ncol(subdat)] = -10 #added thilanka
	  df <- as.matrix(subdat[,12:13]) #added thilanka
	  #df <- as.matrix(cbind(subdat$logP,c(10,-10)))
	  rownames(df) <- subdat$BP
	  df <- df[order(df[,1],decreasing=T),]
	pdf(paste("GO_", short_name_to_save,'.fisher.qvalue_GO_term_BP.pdf',sep=''))
	heatmap.2(df,
		trace="none", # remove trace lines
		col = colorRampPalette(c("blue","white","red"))(21), # color palette
		dendrogram="none", # remove dendrogram
		cexRow=0.8,# row label font size
		labCol=F, # remove column labels
		Rowv=F, # turn off row clustering
		Colv=F, # turn off column clustering
		lmat=rbind(4:3,2:1),
		lwid=c(0.8,4), # dimensions of display array cell widths
		lhei=c(0.8,4)) # dimensions of display cell heights
	dev.off()
	write.table(subdat[,c(9,8,12)],paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term_BP.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	}else if (nrow(subdat) == 1) {
	  subdat <- cbind(subdat,'logP'=0)
	  for(i in 1:nrow(subdat)){
	    if(subdat[i,6] == '-')subdat[i,ncol(subdat)] = log10(subdat[i,8])
	    if(subdat[i,6] == '+')subdat[i,ncol(subdat)] = -log10(subdat[i,8])
	  }
	  subdat[subdat$logP < -10,12] <- -10
	  subdat[subdat$logP > 10,12] <- 10
	  subdat$color <- 10 #added thilanka
	  subdat[1,ncol(subdat)] = -10 #added thilanka
	  write.table(subdat[,c(9,8,12)],paste('GO_',short_name_to_save,'.fisher.qvalue_GO_term_BP.txt',sep=''),row.names=F,col.names=T,quote=F,sep='\t')
	}else{
  print("No_enrichment")
}
	}

GO_pipeline("upregulated_genes_30min.txt","D_30")
GO_pipeline("upregulated_genes_1hr.txt","D_1")
GO_pipeline("upregulated_genes_3hr.txt","D_3")
GO_pipeline("upregulated_genes_6hr.txt","D_6")
GO_pipeline("upregulated_genes_16hr.txt","D_16")
GO_pipeline("upregulated_genes_24hr.txt","D_24")








up_list <- list.files(pattern = "*_cluster.txt",all.files = F, path = "../")
#up_list <- up_list[-c(seq(1,14,1))]
for (i in up_list) {
  
  go_name <- gsub("_cluster.txt","",i)
  print(i)
  GO_pipeline(paste("..",i,sep = "/"),go_name)
}


gsub("_cluster.txt","",up_list[1])


go_name
#GO_pipeline('downregulated_genes_16hr.txt','test')


GO_pipeline("../NDNNUN_cluster.txt","NDNNUN") # no enrichment
#Pathway_pipeline("downregulated_genes_16hr.txt","16_D")

GO_pipeline("../NDNNUU_cluster.txt","NDNNUU") # no enrichment
#Pathway_pipeline("downregulated_genes_1hr.txt","1_D")

GO_pipeline("downregulated_genes_24hr.txt","24_D")
#Pathway_pipeline("downregulated_genes_24hr.txt","24_D")

GO_pipeline("downregulated_genes_30min.txt","30_m_D")
#Pathway_pipeline("downregulated_genes_30min.txt","30_m_D")

GO_pipeline("downregulated_genes_3hr.txt","3_D")
#Pathway_pipeline("downregulated_genes_3hr.txt","3_D")

GO_pipeline("downregulated_genes_6hr.txt","6_D")
#Pathway_pipeline("downregulated_genes_6hr.txt","6_D")
###
GO_pipeline("upregulated_genes_16hr.txt","16_U")
#Pathway_pipeline("upregulated_genes_16hr.txt","16_U")

GO_pipeline("upregulated_genes_1hr.txt","1_U")
#Pathway_pipeline("upregulated_genes_1hr.txt","1_U")

GO_pipeline("upregulated_genes_24hr.txt","24_U")
#Pathway_pipeline("upregulated_genes_24hr.txt","24_U")

GO_pipeline("upregulated_genes_30min.txt","30_m_U")
#Pathway_pipeline("upregulated_genes_30min.txt","30_m_U")

GO_pipeline("upregulated_genes_3hr.txt","3_U")
#Pathway_pipeline("upregulated_genes_3hr.txt","3_U")

GO_pipeline("upregulated_genes_6hr.txt","6_U")
#Pathway_pipeline("upregulated_genes_6hr.txt","6_U")

######################################################################################################################

