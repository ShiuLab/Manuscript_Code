library(rrBLUP)
set.seed(42)
args = commandArgs(trailingOnly=TRUE)
X_file <- args[1] # your genetic matrix, e.g., geno.csv
Y_file <- args[2] # your phenotypic matrix, e.g., pheno.csv
feat_file <- args[3] # selected features or "all" for all the markers in the genetic matrix
trait <- args[4] # the column name of your target trait, or "all" for all the traits in the pheno matrix
cv <- as.numeric(args[5]) # the fold number of the cross-validation scheme
number <- as.numeric(args[6]) # how many times your want to repeat the cross-validation scheme
cvs_file <- args[7] # the CVs file
save_name <- args[8] # a short name for your file to be saved

Y <- read.csv(Y_file, row.names=1) 
X <- read.csv(X_file, row.names=1) 
# make X and Y have the same order of rows
X <- X[rownames(Y),]
cvs <- read.csv(cvs_file, row.names=1)
# Make the relationship matrix from the markers
M=tcrossprod(scale(X))  # centered and scaled XX'
M=M/mean(diag(M))
rownames(M) <- 1:nrow(X)

# Subset X if feat_file is not all
if (feat_file != 'all'){
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character')
  X <- X[FEAT]
  feat_method <- tail(unlist(strsplit(feat_file, '/')), n=1)
}

if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

R2 <- c()
for(i in 1:length(Y)){
	print(names(Y)[i])
	Accuracy <- c()
	Coef <- c()
	for(k in 1:number){
		print(k)
		tst = cvs[,k]
		Coeff <- c()
		yhat <- data.frame(cbind(Y, yhat = 0))
		yhat$yhat <- as.numeric(yhat$yhat)
		row.names(yhat) <- row.names(Y)
		for(j in 1:cv){
			validation <- which(tst==j)
			training <- which(tst!=j)
			yNA <- Y[,i]
			yNA[validation] <- NA # Mask yields for validation set
			# Build rrBLUP model and save yhat for the masked values
			# predict marker effects
			coeff <- mixed.solve(y=Y[training,i], Z=X[training,], K=NULL, method='ML', SE=FALSE, return.Hinv=FALSE)
			Coeff <- rbind(Coeff,coeff$u)
			# predict breeding 
			rrblup <- mixed.solve(y=yNA, K=A.mat(X))
			yhat$yhat[validation] <- rrblup$u[validation]
			# Note that the breeding trait values can be calculated as below. The final accuracy is almost the same as code in above.
			# yhat$yhat[validation] <- (X[validation,] %*% coeff$u)[,1] + coeff$beta
			}
		accuracy <- cor(yhat[,i], yhat$yhat)
		Accuracy <- c(Accuracy,accuracy^2)
		Coef <- rbind(Coef,colMeans(Coeff))
		}
	R2 <- cbind(R2,Accuracy)
	write.csv(Coef,paste('Coef_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=F,quote=F)
	}
colnames(R2) <- names(Y)
write.csv(R2,paste('R2_results_',save_name,'.csv',sep=''),row.names=F,quote=F)
