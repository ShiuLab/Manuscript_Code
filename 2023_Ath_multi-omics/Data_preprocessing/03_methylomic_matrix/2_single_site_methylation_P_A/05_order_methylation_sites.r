dat <- read.table('Methylation_sites_listgenome_wide_618_accessions.txt',head=F,sep='\t')
pos1 <- regexpr('_C',dat[,1])
tem  <- substr(dat[,1],1,pos1-1)
pos2 <- regexpr('_',tem)
chr <- substr(tem,1,pos2-1)
pos <- substr(tem,pos2+1,nchar(tem))
dat <- cbind(dat,chr,pos)
dat <- dat[order(dat[,2],as.numeric(as.character(dat[,3]))),]
write.table(dat,'Methylation_sites_listgenome_wide_618_accessions_ordered.txt',row.names=F,col.names=F,quote=F,sep='\t')
