#####################################
# Make phenotype predictions using rrBLUP
#
# Arguments: [1]   X_file
#            [2]   Y_file
#            [3]   trait (col_name or all)
#            [4]   feat to use
#            [5]   CVFs_file
#            [6]   cvJobNum (cv_##)
#            [7]   optional: save directory
#
#
# Written by: Christina Azodi
# Original: 4.26.17
# Modified: 7.18.2019 to add feature selection
#####################################
library(rrBLUP)
Sys.sleep(sample(1:6, 1))
# Removes all existing variables from the workspace
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

start.time <- Sys.time()

# Read in arguments with 5 as the default PCs to include
if (length(args) < 6) {
  stop("Need 7 arguments: X_file Y_file trait CVFs_file cvJobNum [optional: save directory]", call.=FALSE)
} else if (length(args) < 7) {
  # default output file
  args[7] <- "/mnt/research/ShiuLab/17_GP_SNP_Exp/maize/"
}

X_file <- args[1]
Y_file <- args[2]
trait <- args[3]
feat_file <- args[4]
cvs_file <- args[5]
jobNum <- as.numeric(args[6])
save_dir <- args[7]

print('Loading in data...')
## load the phenotypes and PCs
X <- read.csv(X_file, row.names=1)
Y <- read.csv(Y_file, row.names=1)
cvs <- read.csv(cvs_file, row.names=1)

# Make sure Y is in the same order as X:
Y <- Y[match(rownames(X), rownames(Y)),]


# Subset X if feat_file is not all
if (feat_file != 'all'){
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character')
  print(setdiff(FEAT, colnames(X)))
  X <- X[FEAT]
}

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
  save_name_yhat <- paste('rrBLUP', X_name, names(Y)[i], 'yhat.csv', sep='_')
  save_name_coef <- paste('rrBLUP', X_name, names(Y)[i], 'coef.csv', sep='_')
  y=Y[, names(Y)[i]]
  CV.fold= paste('cv_', toString(jobNum), sep='')
  print(paste('Running ', CV.fold, sep=''))
        
  # load the folds
  tst=cvs[,CV.fold]
  yhat <- data.frame(cbind(y, yhat = 0))
  yhat$yhat <- as.numeric(yhat$yhat)
  row.names(yhat) <- row.names(Y)
  coef <- c()
  for(j in 1:5){
    print(paste('fold =',j))
    test <- which(tst==j)
    train <- which(tst!=j)
    
    test_x <- as.matrix(X[test,])
    test_y <- y[test]
    train_x <- X[train,]
    train_y <- y[train]
    
    trained_model <- mixed.solve(train_y, Z=train_x, K=NULL, SE=FALSE, return.Hinv=FALSE)
    coef_temp <- trained_model$u
    coef <- cbind(coef, coef_temp)
    
    effect_size <- as.matrix(coef_temp)
    pred_y <- (test_x %*% effect_size) 
    pred_y <- pred_y[,1] + trained_model$beta
    
    yhat$yhat[test] <- pred_y

  }
  
  # Calculate accuracy 
  accuracy <- cor(yhat$y, yhat$yhat)
  
  # Save mean coefficients
  coef_df <- as.data.frame(coef)
  coef_df$mean <- rowMeans(coef_df)
  coef_to_save <- t(coef_df)
  coef_to_save <- subset(coef_to_save, row.names(coef_to_save) != 'coef_temp')
  coef_to_save2 <- data.frame('rrBLUP', X_name, names(Y)[i], CV.fold, coef_to_save)
  colnames(coef_to_save2) <- c('model', 'x_file', 'y', 'cv_Num', unlist(colnames(coef_to_save)))
  write.table(coef_to_save2, save_name_coef, append=T, row.names=F, quote=F, col.names=!file.exists(save_name_coef))  
  
  # save  results files
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units='sec')
  df_out <- data.frame('rrBLUP', X_name, names(Y)[i], feat_file, CV.fold, accuracy, time.taken)
  colnames(df_out) <- c('model', 'x_file', 'y','feat', 'cv_Num', 'accuracy_PCC', 'run_time')
  write.table(df_out, 'accuracy.txt', append=T, row.names=F, quote=F, col.names=!file.exists('accuracy.txt'))  
  
  # Save predicted values to save_name
  job_ID <- paste('rrBLUP', X_name, names(Y)[i], CV.fold, sep='_')
  colnames(yhat) <- c('y', job_ID)
  yhat_t <- t(yhat)
  yhat_t <- yhat_t[-c(1), ,drop=FALSE]
  write.table(yhat_t, save_name_yhat, append=T, sep=',', row.names=T, quote=F, col.names=!file.exists(save_name_yhat))
}

unlink('*.dat')
print('Complete')
