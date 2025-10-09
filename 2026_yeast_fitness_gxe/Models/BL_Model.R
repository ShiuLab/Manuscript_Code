##################################################################################################
# Description: Bayesian LASSO for predicting phenotype within a cross-validation (cv) scheme
#
# Arguments:
#   [1] X_file:     genotype matrix
#   [2] Y_file:     phenotype matrix
#   [3] feat_file   selected features or "all" for all features in X_file
#   [4] test_file   file with samples in test set
#   [5] trait:      name of trait (column in phenotype matrix) or "all" for all traits in Y_file
#   [6] cvf_file:   cross-validation scheme file (specifies which sample is part of each cv fold)
#   [7] fold:       cross-validation fold number
#   [8] number:     number of cross-validation repetitions
#   [9] save_name:  output file save name
#   [10] save_dir:  directory to save output file
# 
# Output:
#   [1] Fitted BL model object as .RDS file
#   [2] Feature coefficients file
#   [3] BGLR marker effects file (ETA_1_b.bin)
#   [4] BGLR lambda parameters file (ETA_1_lambda.dat)
#   [5] BGLR mu file (mu.dat)
#   [6] BGLR varE file (varE.dat)
#   [7] R-squared values for each cross-validation repetition and of test set
#   [8] Predicted values for each cross-validation repetition and of test set
#
# Original code by Christina Azodi: 
# https://github.com/ShiuLab/Manuscript_Code/blob/master/2019_expression_GP/scripts/predict_BGLR.R
#
# Modified by: Kenia Segura Abá
##################################################################################################

# Load necessary packages
library(BGLR)
library(data.table)
library(dplyr)

set.seed(23) # for reproducibility

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 10) {
    stop("Need 10 arguments: X_file Y_file feat_file test_file trait cvf_file fold number save_name save_dir", call.=FALSE)
} else {
    X_file <- args[1]
    Y_file <- args[2]
    feat_file <- args[3]
    test_file <- args[4]
    trait <- args[5]
    cvf_file <- args[6]
    fold <- as.numeric(args[7])
    number <- as.numeric(args[8])
    save_name <- args[9]
    save_dir <- args[10]
}

# Read in data
print("Reading in data...")
if (feat_file != "all") {
    print("Pulling features to use...")
    FEAT <- scan(feat_file, what="character") # determine which features are included
    X <- fread(X_file, select=c("ID", FEAT), showProgress=TRUE) # subset genotype data 
    X <- as.matrix(X, rownames=1, colnames=1)
    # X[1:5,1:5]
} else {
    X <- fread(X_file, showProgress=TRUE)
    X <- as.matrix(X, rownames=1, colnames=1)
    # X[1:5,1:5]
}

Y <- read.csv(Y_file, row.names=1)
Test <- scan(test_file, what="character")
cvs <- read.csv(cvf_file, row.names=1)

# Process cross-validation file
cvs_all <- merge(Y, cvs, by="row.names", all.x=TRUE) # merge Y_file and cvs_file
rownames(cvs_all) <- cvs_all$Row.names # set row names to sample name 
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)] # keep only cvs columns
cvs_all[is.na(cvs_all)] = 0 # samples in test file will be "NA", set to 0

# Make sure X and Y have the same order of rows as cvs_all
X <- X[rownames(cvs_all),]
Y <- Y[rownames(cvs_all),]

# Trait or traits to be modelled
if (trait == "all") {
    print("Modeling all traits...")
} else {
    Y <- Y[trait]
}

# Collect results
setwd(save_dir) # set output directory as working directory
file <- "RESULTS_BL.txt"
if (!file.exists(file)) {
    cat("Date", "RunTime", "Trait", "ID", "Alg", "NumInstances", "FeatureNum",
        "CVfold", "CV_rep", "MSE_val", "MSE_val_sd", "MSE_val_se", "r2_val",
        "r2_val_sd", "r2_val_se", "PCC_val", "PCC_val_sd","PCC_val_se",
        "MSE_test", "MSE_test_sd", "MSE_test_se", "r2_test", "r2_test_sd",
        "r2_test_se", "PCC_test", "PCC_test_sd", "PCC_test_se\n", file=file,
        append=FALSE, sep="\t")
} else {message("RESULTS_BL.txt exists") }

# Evaluation metrics
mse <- function(preds, actual){ return(mean((actual-preds)^2)) } # mean squared error
se <- function(vector){ return(sd(vector)/length(vector)) } # standard error
r2_score <- function(preds, actual) {
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss))
} # same as sklearn's r2_score function to compute coefficient of determination (R^2)

`%notin%` <- Negate(`%in%`)

# Run model
for (i in 1:length(Y)) { # loop through selected trait(s)
    print(sprintf("Modeling trait %s...", names(Y)[i]))
    Coef <- c() # feature effects
    pred_val <- c() # predicted value of validation
    pred_test <- c() # predicted value of test set
    Start <- Sys.time() # start time
    for (k in 1:number) { # cross-validation repetitions
        print(sprintf("CV repetition number %i", k))
        tst <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
        Coeff <- c() # estimated feature effects
        yhat_val <- data.frame(Y[i], yhat=0, row.names=row.names(Y)) # dataframe to hold predicted values
        yhat_test <- data.frame(Y[Test,i], row.names=Test)
        colnames(yhat_test) <- colnames(Y[i])
        for (j in 1:fold) { # cross-validion fold number
            print(sprintf("CV fold number %i", j))
            validation <- which(tst==j) # validation set for this fold
            training <- which(tst!=j & tst!=0) # training set is all other data excluding test set
            test <- which(tst==0) # testing set
            yNA <- Y[,i] # label (dependent variable)
            yNA[validation] <- NA # mask validation sample values
            yNA[test] <- NA # mask test sample values
            
            # Build Bayesian LASSO model
            start <- Sys.time()
            ETA <- list(list(X=X, model="BL")) #, saveEffects=TRUE)) # input genotype data and model type
            model <- BGLR(y=as.matrix(yNA), ETA=ETA, verbose=FALSE, nIter=32000, burnIn=3200) #, saveAt = paste("BL_", save_name, "_", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), "_", sep="")) # takes about 13 minutes
            print("Saving model...")
            saveRDS(model, file=paste("BL_", save_name, "_", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), ".RDS", sep=""))
            end <- Sys.time() # end time
            print(sprintf("Model elapsed time: %f", end-start))
            
            # Extract results from model
            Coeff <- rbind(Coeff, model$ETA[[1]]$b) # feature effects
            yhat_val$yhat[validation] <- model$yHat[validation] # predicted labels for validation set
            yhat_test[paste("fold",j,sep="_")] <- model$yHat[test] # collect predicted labels
        }
        Coef <- rbind(Coef, colMeans(Coeff)) # mean feature coefficients
        pred_val <- cbind(pred_val, yhat_val[which(rownames(yhat_val) %notin% Test),]$yhat) # predicted values of validation set
        pred_test <- cbind(pred_test, rowMeans(yhat_test[1:fold+1])) # predicted values of test set
    }
    # save average feature coefficients
    write.csv(Coef, paste("Coef_", save_name, "_", names(Y)[i], ".csv", sep=""),
        row.names=FALSE, quote=FALSE)
    
    # save predicted values
    write.csv(pred_val, paste("Predict_value_cv_", save_name, "_", names(Y)[i],
        ".csv", sep=""), row.names=FALSE, quote=FALSE)
    write.csv(pred_test, paste("Predict_value_test_", save_name, "_", names(Y)[i],
        ".csv", sep=""), row.names=FALSE, quote=FALSE)
    
    # save model performances
    print("Writing model results to file...")
    RunTime <- Sys.time()-Start
    ID <- paste(names(Y)[i], save_name, sep="_")
    NumInstances <- nrow(X)-length(Test)
    mse_val <- apply(pred_val, 2, mse, Y[which(rownames(Y) %notin% Test),i])
    r2_val <- apply(pred_val, 2, r2_score, Y[which(rownames(Y) %notin% Test),i])
    pcc_val <- apply(pred_val, 2, cor, Y[which(rownames(Y) %notin% Test),i])
    mse_test <- apply(pred_test, 2, mse, Y[which(rownames(Y) %in% Test),i])
    r2_test <- apply(pred_test, 2, r2_score, Y[which(rownames(Y) %in% Test),i])
    pcc_test <- apply(pred_test, 2, cor, Y[which(rownames(Y) %in% Test),i])
    cat(paste(Sys.time(), RunTime, names(Y)[i], ID, "Bayesian LASSO",
        NumInstances, ncol(X), fold, number, mean(mse_val), sd(mse_val), se(mse_val),
        mean(r2_val), sd(r2_val), se(r2_val), mean(pcc_val), sd(pcc_val),
        se(pcc_val), mean(mse_test), sd(mse_test), se(mse_test), mean(r2_test),
        sd(r2_test), se(r2_test), mean(pcc_test), sd(pcc_test), se(pcc_test),
        sep="\t"), file=file, append=T, sep="\n")
}

