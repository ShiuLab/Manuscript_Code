install.packages("flexmix")
#####################################################################################################################################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##########################Libraries##############################
library("ggplot2")
library("dplyr")
library("reshape2")
library("flexmix")
#####################################################################################################################################################################
####Looking at gene/gene ks values
paralogs_vsls <- read.table("ks_ka_tomato_paralogs_MCScanx_default_onlylast5col")
#getting rid of the negative values
paralogs_filtrd <- paralogs_vsls[paralogs_vsls[,5] >= 0,]
#now need to get rid of the duplicated paralog combination
paralogs_filtrd$names<- paste(paralogs_filtrd[,1],paralogs_filtrd[,2],sep = "_")
paralogs_filtrd_2 <- unique(paralogs_filtrd, by = "names")
#getting row names
rownames(paralogs_filtrd_2) <- paralogs_filtrd_2$names
#extract KS column
Ks_df <- paralogs_filtrd_2[,5,drop =F]
#Visualizing the data
Ks_df_fin <- Ks_df
Ks_df_fin[Ks_df_fin[,1] >5 ,1] <- 5
ks_dit <- Ks_df_fin %>%
  ggplot(aes(x = V5)) + geom_histogram(bins = 1000) 
pdf("Ks_ditribution.pdf",width=40,height=10)
ks_dit
dev.off()
#####################################################################################################################################################################
#####################################################################################################################################################################
####Looking at block median ks values
paralogs_med <- read.table("meadian_ks_for_syntanic_block.txt",sep = "\t",header = T,row.names = 1)

ks_block_dit <- paralogs_med %>%
  ggplot(aes(x = median_KS)) + geom_histogram(bins = 50) 
pdf("Ks_for_blocks_ditribution.pdf",width=40,height=10)
ks_block_dit

dev.off()
#####################################################################################################################################################################
# Modeling median Ks for syntanic blocks with three ditributions
TPM2 <- paralogs_med
TPM3 <- paralogs_med$median_KS
break_divs <- seq(0,max(TPM3),by=0.05)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
#TPM3.count <- log10(TPM3.count+1) 
bin <- seq(0,max(TPM3),by=0.05)
bin <- bin[-1]
plot(bin,TPM3.count, ylim = c(0,100), xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
#abline(v = 0.68)
#abline(v = 1.05)
#abline(v = 1.35)
abline(v = 0.7698904) # 5% of the  second normal distribution
abline(v = 0.9479036) #5% of the  third normal distribution
df_TPM3 <- data.frame(bin,TPM3.count)
# Model Fitting
# Power distribution model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))
## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height
fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)) + C2*exp(-(bin-mean2)**2/(2*sigma2**2)) + C3*exp(-(bin-mean3)**2/(2*sigma3**2)),
                 data = df_TPM3,
                 start=list(A1=-0.5,C1=1000,C2 =1000,C3 = 1000,mean1=0.68,sigma1=0.1,mean2=1.1,sigma2=0.05,mean3=1.35,sigma3=0.4),algorithm = "port")   ###residual sum-of-squares: 0.8483
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_C2 <- fit_value[3]
fit_C3 <- fit_value[4]
fit_mean1 <- fit_value[5]
fit_sigma1 <- fit_value[6]
fit_mean2 <- fit_value[7]
fit_sigma2 <- fit_value[8]
fit_mean3 <- fit_value[9]
fit_sigma3 <- fit_value[10]
#pdf("keep_TPM_at_0_and_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,100), xlab="median Ks", ylab="Frequncy", xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
curve(fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)), col = 5, add = TRUE)
curve(fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)), col = 6, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) + fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)) + fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)) , col = 2, add = TRUE)

#getting the quntiles
paralogs_med_dist_1 <- paralogs_med[paralogs_med[,1] > 0.45 & paralogs_med[,1] < 0.92 ,]
quantile(paralogs_med_dist_1, probs = c(0.025,0.975)) 

#final fig
pdf("median_Ks_for_syntanic_blocks.pdf",width = 10,height = 10)
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,100), xlab="median Ks", ylab="Frequncy", xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
x = bin
#hist(paralogs_med$median_KS, breaks = 50)
curve(x^fit_A1, col = 3, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
curve(fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)), col = 5, add = TRUE)
curve(fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)), col = 6, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) + fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)) + fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)) , col = 2, add = TRUE)
abline(v = 0.5201583)
#abline(v = 0.9062312 )
abline(v = 0.7698904) # 5% of the  second normal distribution
#abline(v = 0.9479036) #5% of the  third normal distribution
dev.off()

#function to find presentile
pres <- function(mu,sig,z) {
  return(mu + sig*z)
}
pres(0.95913,0.11539,-1.64)
pres(1.35013,0.24526,-1.64)
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

# Modeling median Ks for syntanic blocks with two ditributions
TPM2 <- paralogs_med
TPM3 <- paralogs_med$median_KS
break_divs <- seq(0,max(TPM3),by=0.08)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
#TPM3.count <- log10(TPM3.count+1) 
bin <- seq(0,max(TPM3),by=0.08)
bin <- bin[-1]
plot(bin,TPM3.count, ylim = c(0,150), xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
abline(v = 0.68)
#abline(v = 1.05)
abline(v = 1.35)
#abline(v = 0.7698904) # 5% of the  second normal distribution
#abline(v = 0.9479036) #5% of the  third normal distribution
df_TPM3 <- data.frame(bin,TPM3.count)
# Model Fitting
# Power distribution model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))
## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height

fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)) + C2*exp(-(bin-mean2)**2/(2*sigma2**2)) ,
                 data = df_TPM3,
                 start=list(A1=-0.5,C1=100,C2 =120,mean1=0.6,sigma1=0.1,mean2=1.6,sigma2=0.1),algorithm = "port",control = nls.control(maxiter = 1000))   ###residual sum-of-squares: 0.8483
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_C2 <- fit_value[3]

fit_mean1 <- fit_value[4]
fit_sigma1 <- fit_value[5]
fit_mean2 <- fit_value[6]
fit_sigma2 <- fit_value[7]

#pdf("keep_TPM_at_0_and_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,150), xlab="median Ks", ylab="Frequncy", xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
curve(fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)), col = 5, add = TRUE)
#curve(fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)), col = 6, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) + fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)) , col = 2, add = TRUE)

#getting the quntiles
paralogs_med_dist_1 <- paralogs_med[paralogs_med[,1] > 0.45 & paralogs_med[,1] < 0.92 ,]
quantile(paralogs_med_dist_1, probs = c(0.025,0.975)) 




#function to find presentile
pres <- function(mu,sig,z) {
  return(mu + sig*z)
}
#95th precentile of the first distribution
pres(fit_mean1,fit_sigma1,1.645)
#5th precentile of the second distribution
pres(fit_mean1,fit_sigma1,-1.645)
#final fig
pdf("median_Ks_for_syntanic_blocks_with_two_distributions.pdf",width = 10,height = 10)
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,150), xlab="median Ks", ylab="Frequncy", xaxt="n")
axis(1, at = seq(0, 3, by = 0.1), las=2)
x = bin
#hist(paralogs_med$median_KS, breaks = 50)
curve(x^fit_A1, col = 3, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
curve(fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)), col = 5, add = TRUE)
#curve(fit_C3*exp(-(x-fit_mean3)**2/(2*fit_sigma3**2)), col = 6, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) + fit_C2*exp(-(x-fit_mean2)**2/(2*fit_sigma2**2)) , col = 2, add = TRUE)
abline(v = pres(fit_mean1,fit_sigma1,1.645)) #95th of the  first normal distribution
#abline(v = 0.9062312 )
abline(v = pres(fit_mean1,fit_sigma1,-1.645)) # 5% of the  second normal distribution

dev.off()


Z_of_val <- function(mu,sig,x) {
  #return(mu + sig*z)
  return((x - mu)/sig)
}

Z_of_val(fit_mean2,fit_sigma2,0.89355 )
# Corresponding presentile of the second distribution is 9.1699


