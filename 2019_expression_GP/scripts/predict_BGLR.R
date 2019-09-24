#####################################
# Make phenotype predictions using BGLR
#
# Arguments: [1]   X_file
#            [2]   Y_file
#            [3]   trait (col_name or all)
#            [4]   CVFs_file
#            [5]   cvJobNum (cv_##)
#            [6]   BGLR_model (BL, BRR, BayesA, BayesB)
#            [7]   optional: save directory
#
#
# Written by: Christina Azodi
# Original: 4.26.17
# Modified: 
#####################################
library(BGLR)

# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args) < 6) {
  stop("Need 6 arguments: X_file Y_file trait CVFs_file cvJobNum BGLR_model [optional: save directory]", call.=FALSE)
} else if (length(args) < 7) {
  # default output file
  args[6] <- "/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/"
}

#cvs_file <- '/Volumes/ShiuLab/17_GP_SNP_Exp/maize/01_Data/CVFs.csv'
#X_file <- '/Volumes/ShiuLab/17_GP_SNP_Exp/maize/01_Data/kin.csv'
#Y_file <- '/Volumes/ShiuLab/17_GP_SNP_Exp/maize/01_Data/pheno.csv'
#trait <- 'HT'
#BGLR_model <- 'BL'
#jobNum <- 1
#setwd('/Volumes/ShiuLab/17_GP_SNP_Exp/maize/')

X_file <- args[1]
Y_file <- args[2]
trait <- args[3]
cvs_file <- args[4]
jobNum <- as.numeric(args[5])
BGLR_model <- args[6]
save_dir <- args[7]
df0 <- 5

## load the phenotypes and PCs
X <- read.csv(X_file, row.names=1)
Y <- read.csv(Y_file, row.names=1)
cvs <- read.csv(cvs_file, row.names=1)

# Make sure Y is in the same order as X:
Y <- Y[match(rownames(X), rownames(Y)),]

if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

# Make output directory
setwd(save_dir)
X_name <- tail(unlist(strsplit(X_file, '/')), n=1)
X_name <- gsub('.csv','', X_name)
X_name <- gsub('.txt','', X_name)

for(i in 1:length(Y)){
  #i <- 1
  save_name <- paste(BGLR_model, X_name, names(Y)[i], 'yhat.csv', sep='_')
  save_name_coef <- paste(BGLR_model, X_name, names(Y)[i], 'coef.csv', sep='_')
  y=Y[, names(Y)[i]]
  CV.fold= paste('cv_', toString(jobNum), sep='')
  print(CV.fold)
  ETA=list(list(X=X,model=BGLR_model)) 
  
  # load the folds
  tst=cvs[,CV.fold]
  yhat <- data.frame(cbind(y, yhat = 0))
  yhat$yhat <- as.numeric(yhat$yhat)
  row.names(yhat) <- row.names(Y)
  coef <- c()
  for(j in 1:5){
    print(paste('fold =',j))
    test <- which(tst==j)
    yNA <- y
    yNA[test] <- NA # Mask yields for validation set
    fm=BGLR(y=yNA,ETA=ETA,verbose=FALSE,nIter=12000,burnIn=2000)
    yhat$yhat[test] <- fm$yHat[test]
    coef_temp <- fm$ETA[[1]]$b
    coef <- cbind(coef, coef_temp)
  }
  
  unlink('*.dat')
  
  # Calculate accuracy 
  accuracy <- cor(yhat$y, yhat$yhat)
  
  # Save mean coefficients
  coef_df <- as.data.frame(coef)
  coef_df$mean <- rowMeans(coef_df)
  coef_to_save <- t(coef_df)
  coef_to_save <- subset(coef_to_save, row.names(coef_to_save) != 'coef_temp')
  coef_to_save2 <- data.frame(BGLR_model, X_name, names(Y)[i], CV.fold, coef_to_save)
  colnames(coef_to_save2) <- c('model', 'x_file', 'y', 'cv_Num', unlist(colnames(coef_to_save)))
  write.table(coef_to_save2, save_name_coef, append=T, row.names=F, quote=F, col.names=!file.exists(save_name_coef))  
  
  # save  results files
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units='sec')
  df_out <- data.frame(BGLR_model, X_name, names(Y)[i], CV.fold, accuracy, time.taken)
  colnames(df_out) <- c('model', 'x_file', 'y', 'cv_Num', 'accuracy_PCC', 'run_time')
  write.table(df_out, 'accuracy.txt', append=T, row.names=F, quote=F, col.names=!file.exists('accuracy.txt'))  
  
  # Save predicted values to save_name
  job_ID <- paste(BGLR_model, X_name, names(Y)[i], CV.fold, sep='_')
  colnames(yhat) <- c('y', job_ID)
  yhat_t <- t(yhat)
  yhat_t <- yhat_t[-c(1), ,drop=FALSE]
  write.table(yhat_t, save_name, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(save_name))
}


print('Complete')
