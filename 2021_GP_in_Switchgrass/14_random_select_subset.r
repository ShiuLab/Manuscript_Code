library(rrBLUP)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
geno_file = args[1]
start = as.numeric(args[2])
stop = as.numeric(args[3])
step = as.numeric(args[4])
total_number <- as.numeric(args[5]) + 1
for(i in seq(start,stop,step)){
	geno <- fread(geno_file,select=c(1,sample(seq(2,total_number),i,replace = FALSE)))
	fwrite(geno,paste('geno_',i,'.csv',sep=''),sep = ",",quote=FALSE)
	}




