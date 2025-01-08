setwd("~/Documents/prj_01_S_lycopersicum_trichome/calculating TPM")
dat <- read.table('TMP_for_all_tissues.txt',head=T,stringsAsFactors=F,sep='\t')
TPM <- c()
for(i in 2:ncol(dat)){
	TPM <- c(TPM,dat[,i])
	}
##########################################
###                                    ###
### keep TPM=0, and don't log y-axis   ###
###                                    ###
##########################################
TPM2 <- log10(TPM+1)
TPM3 <- TPM2[TPM2<=2.5]
break_divs <- seq(0,2.5,by=0.01)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
bin <- seq(0.01,2.5,by=0.01)
plot(bin,TPM3.count)
df_TPM3 <- data.frame(bin,TPM3.count)
# Model Fitting
# Exponential model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))

?nls## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height
fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)),data = df_TPM3,start=list(A1=-0.5,C1=6000,mean1=1.25,sigma1=1),algorithm='port')   ####residual sum-of-squares: 3.283e+09
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_mean1 <- fit_value[3]
fit_sigma1 <- fit_value[4]
pdf("keep_TPM_at_0_not_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,max(df_TPM3$TPM3.count)))
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 2, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
#lines(df_TPM3$bin, predict(fit_mixed), col = 2)
dev.off()

##########################################
###                                    ###
### remove TPM=0, and don't log y-axis ###
###                                    ###
##########################################
TPM2 <- log10(TPM+1)
TPM3 <- TPM2[TPM2<=2.5 & TPM2>0]
break_divs <- seq(0,2.5,by=0.01)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
bin <- seq(0.01,2.5,by=0.01)
plot(bin,TPM3.count)
df_TPM3 <- data.frame(bin,TPM3.count)
# Model Fitting
# Power distribution model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))
## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height
fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)),data = df_TPM3,start=list(A1=-0.5,C1=6000,mean1=1.25,sigma1=1),algorithm='port')   ####residual sum-of-squares: 908830311
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_mean1 <- fit_value[3]
fit_sigma1 <- fit_value[4]
pdf("Remove_TPM_at_0_not_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,max(df_TPM3$TPM3.count)))
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 2, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
#lines(df_TPM3$bin, predict(fit_mixed), col = 2)
dev.off()

####################################
###                              ###
### remove TPM=0, and log y-axis ###
###                              ###
####################################
TPM2 <- log10(TPM+1)
TPM3 <- TPM2[TPM2<=2.5 & TPM2>0]
break_divs <- seq(0,2.5,by=0.01)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
TPM3.count <- log10(TPM3.count)
bin <- seq(0.01,2.5,by=0.01)
plot(bin,TPM3.count)
df_TPM3 <- data.frame(bin,TPM3.count)
# Model Fitting
# Power distribution model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))
## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height
fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)),data = df_TPM3,start=list(A1=-0.5,C1=6000,mean1=1.25,sigma1=1),algorithm='port')   ###residual sum-of-squares: 1.86
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_mean1 <- fit_value[3]
fit_sigma1 <- fit_value[4]
pdf("Remove_TPM_at_0_and_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,5))
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 2, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
#lines(df_TPM3$bin, predict(fit_mixed), col = 2)
dev.off()

##################################
###                            ###
### keep TPM=0, and log y-axis ###
###                            ###
##################################
##This is the final model that I used
TPM2 <- log10(TPM+1)
TPM3 <- TPM2[TPM2<=2.5]
break_divs <- seq(0,2.5,by=0.01)
TPM3.cut <- cut(TPM3,break_divs,right=FALSE)
TPM3.count <- as.numeric(table(TPM3.cut))
TPM3.count <- log10(TPM3.count)
bin <- seq(0.01,2.5,by=0.01)
plot(bin,TPM3.count)
df_TPM3 <- data.frame(bin,TPM3.count) 

# Model Fitting
# Power distribution model: y = bin^a, here a is the power, and the parameter we need to figure out
# Gaussian Model or normal distribution model: y = C0*exp(-(bin-mean0)**2/(2*sigma0**2))
## C0 is the height of the peak; mean0 is the x axis value of the peak; sigma0 is half of the x axis distance between two point at half of the peak height
fit_mixed <- nls(TPM3.count ~ bin^A1 + C1*exp(-(bin-mean1)**2/(2*sigma1**2)),data = df_TPM3,start=list(A1=-0.5,C1=6000,mean1=1.25,sigma1=1),algorithm='port')   ###residual sum-of-squares: 0.8483
# Get model parameters 
fit_value <- fit_mixed$m$getAllPars()
fit_A1 <- fit_value[1]
fit_C1 <- fit_value[2]
fit_mean1 <- fit_value[3]
fit_sigma1 <- fit_value[4]
pdf("keep_TPM_at_0_and_log_y-axis.pdf")
plot(df_TPM3$bin,df_TPM3$TPM3.count,ylim=c(0,5))
x = bin
curve(x^fit_A1, col = 3, add = TRUE)
curve(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 4, add = TRUE)
curve(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)), col = 2, add = TRUE)
abline( v=0.4751045)
#lines(df_TPM3$bin, predict(fit_mixed), col = 2)
dev.off()

# lookig at summery stats
confint(fit_mixed)
summary(fit_mixed)
plot(x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) ~ df_TPM3$bin,ylim=c(0,5))
lines(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) ~ df_TPM3$bin)
lines(x^fit_A1  ~ df_TPM3$bin)
abline( v=0.4751045)

#extract and draw a histogram of the TPM values that are less than 0.4751045

#make a dataframe with x^fit_A1 + fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) and corresponding df_TPM3$bin values
plot(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)) ~ df_TPM3$bin)
df_sim_TPM <- data.frame(fit_C1*exp(-(x-fit_mean1)**2/(2*fit_sigma1**2)),df_TPM3$bin)
head(df_sim_TPM)
bin_vec <- c()
for(i in 1:nrow(df_sim_TPM)){
	bin_vec <- c(bin_vec,rep(df_sim_TPM[i,2],df_sim_TPM[i,1]*100))
	}
hist(bin_vec,breaks=100,main="Histogram of TPM values that are less than 0.4751045",xlab="TPM",ylab="Frequency")
#get the precentile = 0.4751045 of the bin_vec

precentile <- quantile(bin_vec,0.5)
#make a dataframe channg 0-100 quntile and get corrosponding value
df_precentile <- data.frame(0:100,quantile(bin_vec,0:100/100))

log10(10+1)


TPM_tresh_at_10 <- 10^(0.47) -1

TPM_tresh_at_25 <- 10^(0.69) -1

TPM_tresh_at_50 <- 10^(1.26) -1

TPM_tesh_at_75 <- 10^(1.83) -1

TPM_tesh_at_90 <- 10^(2.21) -1
