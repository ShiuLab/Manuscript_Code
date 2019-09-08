#####################################
# Make phenotype predictions using rrBLUP
#
# Arguments: [1] id (i.e. wheat_599_CIMMYT)
#            [2] JobNum (use PSB JobArray)
#            [3] Traits to apply to (default = all)
#            [4] Output directory
#
#
# Written by: Christina Azodi
# Original: 4.26.17
# Modified: 
#####################################
library(rrBLUP)


# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args)==0) {
  stop("Need at least one argument (ID)", call.=FALSE)
} 

id = args[1]
jobNum = as.numeric(args[2])
trait = args[3]
save_dir = args[4]

## load the phenotypes and PCs
setwd(paste("/mnt/home/azodichr/03_GenomicSelection/", id, '/', sep=''))
Y <- read.csv('01_Data/pheno.csv', row.names=1)
X <- read.csv('01_Data/geno.csv', row.names=1)
cvs <- read.csv('01_Data/CVFs.csv', row.names=1)

# Make sure Y is in the same order as X:
Y <- Y[match(rownames(X), rownames(Y)),]

# Make the relationship matrix from the markers
M=tcrossprod(scale(X))  # centered and scaled XX'
M=M/mean(diag(M))
rownames(M) <- 1:nrow(X)

# Make output directory
setwd(save_dir)
dir.create(id)
setwd(id)
dir.create('03_rrBLUP')
setwd('03_rrBLUP')


# For trait in Y dataframe
for(i in 1:length(Y)){
  names(Y)[i]
  dir.create(paste('trait_',names(Y)[i], sep=''))
  setwd(paste('trait_',names(Y)[i], sep=''))
  dir.create('output')
  y=Y[, names(Y)[i]]
  CV.fold= paste('cv_', toString(jobNum-1), sep='')
  #CV.fold = 'cv_2'
  
  if(CV.fold =='cv_0'){
    # fit model to the entire data set and save model
    df <- data.frame(y=y,gid=1:nrow(X)) # Set up dataframe with traits and genotype labels (same order as in A1) 
    rrblup <- kin.blup(df,K=M,geno="gid",pheno='y') #optional parameters: fixed effects, gaussian kernel, covariates
    write.csv(cbind(y, rrBLUP_pred = rrblup$g), file='full_pred.csv')
    save(rrblup,file='full_model.RData')
  }
  
  else{
    # fit the model using the training sets designated by CV.
    tst=cvs[,CV.fold]
    yhat <- data.frame(cbind(y, yhat = 0))
    yhat$yhat <- as.numeric(yhat$yhat)
    row.names(yhat) <- row.names(Y)
    
    for(i in 1:5){
      # Make training (TRN) and testing (TST) dfs
      test <- which(tst==i)
      yNA <- y
      yNA[test] <- NA # Mask yields for validation set
      df <- data.frame(y=yNA,gid=1:nrow(X)) # Set up dataframe with traits and genotype labels (same order as in A1) 
      
      # Build rrBLUP model and save yhat for the masked values
      rrblup <- kin.blup(df,K=M,geno="gid",pheno='y') #optional parameters: fixed effects, gaussian kernel, covariates
      yhat$yhat[test] <- rrblup$g[test]
      
      #save(rrblup,file=paste(CV.fold,'.',i,'.RData', sep=''))
    }
    
    write.table(yhat, paste('output/',CV.fold,'.csv', sep=''), sep=',', row.names=FALSE, col.names=TRUE)
    accuracy <- cor(yhat$y, yhat$yhat)
    end.time <- Sys.time()
    time.taken <- difftime(end.time, start.time, units='sec')
    df_out <- data.frame(CV.fold, accuracy, time.taken)
    write.table(df_out, 'accuracy.csv', append=TRUE, sep=',', row.names=FALSE, col.names=FALSE)  
  }
  setwd('../')
}

print('Complete')
