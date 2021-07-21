library(rrBLUP)
library(data.table)
set.seed(42)
args = commandArgs(TRUE)
X_file <- args[1] # your genetic matrix, e.g., geno.csv
Y_file <- args[2] # your phenotypic matrix, e.g., pheno.csv
feat_file <- args[3] # selected features or "all" for all the markers in the genetic matrix
trait <- args[4] # the column name of your target trait, or "all" for all the traits in the pheno matrix
cv <- as.numeric(args[5]) # the fold number of the cross-validation scheme
number <- as.numeric(args[6]) # how many times your want to repeat the cross-validation scheme
cvs_file <- args[7] # the CVs file
save_name <- args[8] # a short name for your file to be saved

Y <- read.csv(Y_file, row.names=1) 
cvs <- read.csv(cvs_file, row.names=1)

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

# make sure X and Y have the same order of rows as cvs
X <- X[rownames(cvs),]
Y <- Y[rownames(cvs),]

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
			# coeff <- mixed.solve(y=Y[training,i], Z=X[training,], K=NULL, SE=FALSE, return.Hinv=FALSE)
			# Coeff <- rbind(Coeff,coeff$u)
			# # predict breeding 
			# rrblup <- mixed.solve(y=yNA, K=A.mat(X))
			# yhat$yhat[validation] <- rrblup$u[validation]
			xTraining <- X[training,]
			xvalidation <- X[validation,]
			yTraining <- Y[training,i]
			yvalidation <- Y[validation,i]
			train_set <- data.frame(cbind(yTraining, xTraining))
			validation_set <- data.frame(cbind(yvalidation, xvalidation))
			fm2 <- lm(yTraining ~ ., data=train_set)
			pred <- predict(fm2, newdata=(validation_set), interval="confidence")
			results <- data.frame(pred[,'fit'])
			yhat$yhat[validation] <- results$pred....fit..
			Coeff <- rbind(Coeff,fm2$coefficients[2:length(fm2$coefficients)])
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
