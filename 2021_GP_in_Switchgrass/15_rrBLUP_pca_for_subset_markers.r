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

Y <- read.csv(Y_file, row.names=1) 
Test <- scan(test_file, what='character')

# Subset X if feat_file is not all
if (feat_file != 'all'){
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character')
  X <- fread(X_file,select=c('ID',FEAT))
  fwrite(X,paste('geno',feat_file,'.csv',sep=''),sep = ",",quote=FALSE)
  feat_method <- tail(unlist(strsplit(feat_file, '/')), n=1)
  X <- read.csv(paste('geno',feat_file,'.csv',sep=''), row.names=1) 
} else{
	X <- as.matrix(fread(X_file),rownames=1)
	}

# get rid of markers with no variation among individuals
X_tem <- c()
col_names <- c()
for(i in 1:ncol(X)){
	if(length(unique(X[,i:i])) > 1){
		X_tem <- cbind(X_tem,as.matrix(X[,i:i]))
		col_names <- c(col_names,colnames(X)[i])
		}
	else{
		print(i)
		}
	if(i%%1000 == 0){
		print(i)
		}
	}
rownames(X_tem) <- rownames(X)
colnames(X_tem) <- col_names
X = X_tem	

#### Calculate the G matrix
X=scale(X) # Centers (subtract the column means) and Scales (dividing the centered columns by their stdev)
G=tcrossprod(X) # Take the cross product X transpose
G=G/mean(diag(G))

#### Calculate eigenvalues
EVD=eigen(G)
rownames(EVD$vectors)=rownames(G)
save(EVD,file='EVD.RData')
write.csv(EVD$vectors,paste("PCA_matrix_",X_file,".csv",sep=""),row.names=T,quote=F)
write.csv(EVD$vectors[,1:5],paste("PCA5_",X_file,".csv",sep=""),row.names=T,quote=F)

X = EVD$vectors[,1:5]

# make sure X and Y have the same order of rows
X <- X[rownames(Y),]
cvs <- read.csv(cvs_file, row.names=1)
cvs_all <- merge(Y,cvs,by="row.names",all.x=TRUE)
rownames(cvs_all) <- cvs_all$Row.names
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)]
cvs_all[is.na(cvs_all)] = 0

if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

R2_cv <- c()
R2_test <- c()
Predict_validation <- c()
Predict_test <- c()
for(i in 1:length(Y)){
	print(names(Y)[i])
	Accuracy_CV <- c()
	Accuracy_test <- c()
	Coef <- c()
	for(k in 1:number){
		print(k)
		tst = cvs_all[,k]
		Coeff <- c()
		y_test <- c()
		yhat <- data.frame(cbind(Y, yhat = 0))
		yhat$yhat <- as.numeric(yhat$yhat)
		row.names(yhat) <- row.names(Y)
		for(j in 1:cv){
			validation <- which(tst==j)
			training <- which(tst!=j & tst!=0)
			test <- which(tst==0)
			yNA <- Y[,i]
			yNA[validation] <- NA # Mask yields for validation set
			yNA[test] <- NA # Mask yields for test set
			# Build rrBLUP model and save yhat for the masked values
			# predict marker effects
			coeff <- mixed.solve(y=Y[training,i], Z=X[training,], K=NULL, method='ML', SE=FALSE, return.Hinv=FALSE)
			Coeff <- rbind(Coeff,coeff$u)
			effect_size <- as.matrix(coeff$u)
			# predict breeding 
			# rrblup <- mixed.solve(y=yNA, K=A.mat(X))
			# yhat$yhat[validation] <- rrblup$u[validation]
			yhat$yhat[validation] <- (as.matrix(X[validation,]) %*% effect_size)[,1] + c(coeff$beta)
			# yhat$yhat[test] <- rrblup$u[test]
			yhat$yhat[test] <- (as.matrix(X[test,]) %*% effect_size)[,1] + c(coeff$beta)
			y_test <- cbind(y_test,yhat$yhat)
			}
		accuracy_CV <- cor(yhat[which(tst!=0),i], yhat$yhat[which(tst!=0)])
		Accuracy_CV <- c(Accuracy_CV,accuracy_CV^2)
		y_test <- cbind(y_test,rowMeans(y_test))
		accuracy_test <- cor(yhat[which(tst==0),i], y_test[which(tst==0),ncol(y_test)])
		Accuracy_test <- c(Accuracy_test,accuracy_test^2)
		Coef <- rbind(Coef,colMeans(Coeff))
		pred_val <- cbind(pred_val,yhat$yhat[which(tst!=0)])
		pred_test <- cbind(pred_test,y_test[which(tst==0),ncol(y_test)])
		}
	R2_cv <- cbind(R2_cv,Accuracy_CV)
	R2_test <- cbind(R2_test,Accuracy_test)
	write.csv(Coef,paste('Coef_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=F,quote=F)
	colnames(pred_val) <- paste(names(Y)[i],'_',1:number,sep='')
	Predict_validation <- cbind(Predict_validation,pred_val)
	colnames(pred_test) <- paste(names(Y)[i],'_',1:number,sep='')
	Predict_test <- cbind(Predict_test,pred_test)
	}
colnames(R2_cv) <- names(Y)
write.csv(R2_cv,paste('R2_cv_results_',save_name,'.csv',sep=''),row.names=F,quote=F)
colnames(R2_test) <- names(Y)
write.csv(R2_test,paste('R2_test_results_',save_name,'.csv',sep=''),row.names=F,quote=F)
rownames(Predict_validation) <- rownames(X)[which(tst!=0)]
write.csv(Predict_validation,paste('Predict_value_cv_',save_name,'.csv',sep=''),row.names=T,quote=F)
rownames(Predict_test) <- rownames(X)[test]
write.csv(Predict_test,paste('Predict_value_test_',save_name,'.csv',sep=''),row.names=T,quote=F)