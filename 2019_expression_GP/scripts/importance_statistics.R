library('data.table')

# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

file_name <- args[1]
type <- args[2]

if(type == 'imp'){
  d <- fread(args[1], sep='\t', header=F)
  names(d) <- c('ID', 'imp')
  quantile <- ecdf(d$imp)
  d$perc <- quantile(d$imp)
  d$zscore <- (d$imp - mean(d$imp))/ sd(d$imp)
  d$rank <- NA
  d$rank[order(d$imp)] <- nrow(d):1
  print(head(d))
  write.csv(d, paste(file_name, '_stats.csv', sep=''), quote=F, row.names=F)
}

if(type == 'coef'){
  d <- fread(args[1], sep=' ', header=T)
  d2 <- as.data.frame(colMeans(abs(d[,5:length(d)])))
  names(d2) <- c('imp')
  quantile <- ecdf(d2$imp)
  d2$perc <- quantile(d2$imp)
  d2$zscore <- (d2$imp - mean(d2$imp))/ sd(d2$imp)
  d2$rank <- NA
  d2$rank[order(d2$imp)] <- nrow(d2):1
  print(head(d2))
  write.csv(d2, paste(file_name, '_stats.csv', sep=''), quote=F, row.names=T)
}



