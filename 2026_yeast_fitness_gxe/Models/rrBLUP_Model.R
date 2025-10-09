# Description: rrBLUP regression script
# Original code written by Peipei Wang (https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_GP_in_Switchgrass/13_rrBLUP_training_test_split_fread_predict_values.r)
# Modified by: Kenia Segura Abá

library(rrBLUP)
library(data.table)
set.seed(42)
args = commandArgs(trailingOnly=TRUE)
X_file <- args[1] # your genetic matrix, e.g., geno.csv
Y_file <- args[2] # your phenotypic matrix, e.g., pheno.csv
feat_file <- args[3] # selected features or "all" for all the markers in the genetic matrix
trait <- args[4] # the column name of your target trait, or "all" for all the traits in the pheno matrix
test_file <- args[5] # file with individuals in test set
cv <- as.numeric(args[6]) # the fold number of the cross-validation scheme
number <- as.numeric(args[7]) # how many times your want to repeat the cross-validation scheme
cvs_file <- args[8] # the CVs file
save_name <- args[9]
scale <- args[10] # whether to scale feature table or not

# Coefficient of determination (R^2) function and standard error function
r2_score <- function(preds, actual) {
	# This function is comparable to sklearn's r2_score function
	# It computes the coefficient of determination (R^2)
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss)) # return R^2 value
}

se <- function(x) { return(sd(x)/sqrt(length(x))) }

# Read in label and test set files
Y <- read.csv(Y_file, row.names=1) 
Test <- scan(test_file, what='character')

# if file is larger than 10Mb, using fread to read the file
if(file.size(X_file) > 10*1024*1024){
	# Subset X if feat_file is not all
	if (feat_file != 'all'){
		print('Pulling features to use...')
		FEAT <- scan(feat_file, what='character')
		X <- as.matrix(fread(X_file,select=c('ID',FEAT)), rownames=1)
	} else{
		X <- as.matrix(fread(X_file),rownames=1)
	}
}else{
	X <- read.csv(X_file, row.names=1) 
	# Subset X if feat_file is not all
	if (feat_file != 'all'){
		print('Pulling features to use...')
		FEAT <- scan(feat_file, what='character')
		X <- X[FEAT]
	}
}

if (scale == 'y') {
	print('Scaling features...')
	X <- scale(X)
}

# Cross-validation fold assignments
cvs <- read.csv(cvs_file, row.names=1)
cvs_all <- merge(Y,cvs,by="row.names",all.x=TRUE)
rownames(cvs_all) <- cvs_all$Row.names
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)]
cvs_all[is.na(cvs_all)] = 0

# make sure X and Y have the same order of rows as cvs_all
X <- X[rownames(cvs_all),]
Y <- Y[rownames(cvs_all),]


if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

# File to save results to
file <- "RESULTS_rrblup.txt"
if (!file.exists(file)) {
    cat("Date", "Trait", "ID", "NumInstances", "NumFeatures", "CV-Fold", "NumRepetitions", "R2_val", "R2_val_sd", "R2_val_se",  "PCC_val", "PCC_val_sd", "PCC_val_se",  "R2_test", "R2_test_sd", "R2_test_se", "PCC_test", "PCC_test_sd", "PCC_test_se\n", file=file, append=FALSE, sep="\t")
} else {
	message("RESULTS_rrblup.txt exists in: ")
	print(getwd())
}

for(i in 1:length(Y)){
	print(names(Y)[i])
	Coef <- c() # Model marker effects per trait
	Error <- c() # Model residual error term (Ve) per trait
	Beta <- c() # Model fixed effects (β) per trait
	pred_val <- c() # Store predicted values for validation set
	pred_test <- c() # Store predicted values for test set
	for(k in 1:number){ # Repeat the CV scheme "number" times
		print(k)
		tst = cvs_all[,k]
		Coeff <- c()
		Errors <- c() # Model residual error term (Ve) per j-cv repetition
		Betas <- c() # Model fixed effects (β) per j-cv repetition
		y_test <- c()
		yhat <- data.frame(cbind(Y, yhat = 0))
		yhat$yhat <- as.numeric(yhat$yhat)
		row.names(yhat) <- row.names(Y)
		for(j in 1:cv){ # Each CV fold
			validation <- which(tst==j)
			training <- which(tst!=j & tst!=0)
			test <- which(tst==0)
			yNA <- Y[,i]
			yNA[validation] <- NA # Mask yields for validation set
			yNA[test] <- NA # Mask yields for test set
			# Build rrBLUP model and save yhat for the masked values
			# predict marker effects
			coeff <- mixed.solve(y=Y[training,i], Z=X[training,], K=NULL, SE=FALSE, return.Hinv=FALSE)
			Coeff <- rbind(Coeff,coeff$u)
			Errors <- rbind(Errors, coeff$Ve) # Model residual error term (Ve) per cv fold
			Betas <- rbind(Betas, coeff$beta) # Model fixed effects term (β) per cv fold
			effect_size <- as.matrix(coeff$u)
			# predict breeding 
			yhat$yhat[validation] <- (as.matrix(X[validation,]) %*% effect_size)[,1] + c(coeff$beta)
			yhat$yhat[test] <- (as.matrix(X[test,]) %*% effect_size)[,1] + c(coeff$beta)
			y_test <- cbind(y_test,yhat$yhat)
		}
		y_test <- cbind(y_test,rowMeans(y_test))
		Coef <- rbind(Coef,colMeans(Coeff))
		Error <- rbind(Error, colMeans(Errors)) # Model average residual errors across folds
		Beta <- rbind(Beta, colMeans(Betas)) # Model average fixed effects across folds
		pred_val <- cbind(pred_val,yhat$yhat[which(tst!=0)])
		pred_test <- cbind(pred_test,y_test[which(tst==0),ncol(y_test)])
	}
	write.csv(Coef,paste('Coef_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Coefficients
	write.csv(Error,paste('Error_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average residual errors across folds
	write.csv(Beta,paste('Beta_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average fixed effects across folds

	colnames(pred_val) <- paste(names(Y)[i],'_',1:number,sep='')
	rownames(pred_val) <- rownames(X)[-test]
	colnames(pred_test) <- paste(names(Y)[i],'_',1:number,sep='')
	rownames(pred_test) <- rownames(X)[test]
	write.csv(pred_val, paste('Predict_value_cv_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)
	write.csv(pred_test, paste('Predict_value_test_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)

	# Calculate model performance and save results to file
	pcc.cv <- apply(pred_val, 2, function(x) cor(x,Y[rownames(pred_val),names(Y)[i]], method="pearson"))
	pcc.test <- apply(pred_test, 2, function(x) cor(x,Y[rownames(pred_test),names(Y)[i]], method="pearson"))
	r2.cv <- apply(pred_val, 2, function(x) r2_score(x,Y[rownames(pred_val),names(Y)[i]]))
	r2.test <- apply(pred_test, 2, function(x) r2_score(x,Y[rownames(pred_test),names(Y)[i]]))
	write.csv(pcc.cv, paste('PCC_cv_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)
	write.csv(pcc.test, paste('PCC_test_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)
	write.csv(r2.cv, paste('R2_cv_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)
	write.csv(r2.test, paste('R2_test_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=TRUE,quote=FALSE)

	print("Saving results...")
	cat(as.character(Sys.time()), "\t", file=file, append=TRUE, sep="")
	cat(names(Y)[i], "\t", file=file, append=TRUE, sep="")
	cat(save_name, "\t", file=file, append=TRUE, sep="")
	cat(dim(X)[1], "\t", file=file, append=TRUE, sep="")
	cat(dim(X)[2], "\t", file=file, append=TRUE, sep="")
	cat(cv, "\t", file=file, append=TRUE, sep="")
	cat(number, "\t", file=file, append=TRUE, sep="")
	cat(mean(r2.cv), "\t", file=file, append=TRUE, sep="")
	cat(sd(r2.cv), "\t", file=file, append=TRUE, sep="")
	cat(se(r2.cv), "\t", file=file, append=TRUE, sep="")
	cat(mean(pcc.cv), "\t", file=file, append=TRUE, sep="")
	cat(sd(pcc.cv), "\t", file=file, append=TRUE, sep="")
	cat(se(pcc.cv), "\t", file=file, append=TRUE, sep="")
	cat(mean(r2.test), "\t", file=file, append=TRUE, sep="")
	cat(sd(r2.test), "\t", file=file, append=TRUE, sep="")
	cat(se(r2.test), "\t", file=file, append=TRUE, sep="")
	cat(mean(pcc.test), "\t", file=file, append=TRUE, sep="")
	cat(sd(pcc.test), "\t", file=file, append=TRUE, sep="")
	cat(se(pcc.test), "\n", file=file, append=TRUE, sep="")
}

