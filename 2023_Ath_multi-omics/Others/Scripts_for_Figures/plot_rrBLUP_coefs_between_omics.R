###################################################
# Density plots of Genomic Prediction Analysis of
# flowering time to compare top features (genomic, 
# transcriptomic, methylomic [gene body]) by rank
# of rrBLUP feature coefficients.
#
# Written by: Kenia Segura Abá
####################################################

library("dplyr")
library("ggplot2")
library("ggExtra")
library("data.table")
library("MASS")
library("viridis")

# Read in feature coefficients for genomic, transcriptomic, and methylomic data
print("Reading in files...")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Data/") # Path to files
T <- fread("Rank_Coef_Expression_FT10_mean.csv")
M <- fread("Rank_Coef_GB_methy_FT10_mean.csv")
G <- fread("Rank_Coef_SNP_FT10_mean.csv")
G2 <- fread("Rank_max_per_gene_Coef_SNP_FT10_mean.csv") # SNP w/max value represents the gene

# Rename G2 column names
G2 <- rename(G2, abs_coef = max_abs_coef, coef = max_coef)
G2$rank <- seq(1, nrow(G2)) # Add rank column

# Split methylomic data by type of gene body methylation
mCG <- M[which(M$type=="mCG"),]
mCHH <- M[which(M$type=="mCHH"),]
mCHG <- M[which(M$type=="mCHG"),]

# Add feature normalized rank (1-0 descending importance)
T$rank_norm <- 1 - ((T$rank - 1)/(nrow(T) - 1))
M$rank_norm <- 1 - ((M$rank - 1)/(nrow(M) - 1))
G$rank_norm <- 1 - ((G$rank - 1)/(nrow(G) - 1))
G2$rank_norm <- 1 - ((G2$rank - 1)/(nrow(G2) - 1))
mCG$rank_norm <- 1 - ((mCG$rank - 1)/(nrow(mCG) - 1))
mCHH$rank_norm <- 1 - ((mCHH$rank - 1)/(nrow(mCHH) - 1))
mCHG$rank_norm <- 1 - ((mCHG$rank - 1)/(nrow(mCHG) - 1))

# Plot feature coefficient distribution
T.quant <- quantile(T$coef, c(0.95,0.05,0.99,0.01))
g <- ggplot(T, aes(coef)) + geom_density() + theme_bw() + xlab("Expression Feature Coefficients") +
            geom_vline(xintercept=T.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=T.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=T.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=T.quant[4], color="blue", linetype="dashed")
G.quant <- quantile(G$coef, c(0.95,0.05,0.99,0.01))
g2 <- ggplot(G, aes(coef)) + geom_density() + theme_bw() + xlab("SNP Feature Coefficients") +
            geom_vline(xintercept=G.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=G.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=G.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=G.quant[4], color="blue", linetype="dashed")
M.quant <- quantile(M$coef, c(0.95,0.05,0.99,0.01))
p <- ggplot(M, aes(coef)) + geom_density() + theme_bw() + xlab("GB Methylation Feature Coefficients") +
            geom_vline(xintercept=M.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=M.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=M.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=M.quant[4], color="blue", linetype="dashed")
mCG.quant <- quantile(mCG$coef, c(0.95,0.05,0.99,0.01))
l <- ggplot(mCG, aes(coef)) + geom_density() + theme_bw() + xlab("GB Methylation CG Feature Coefficients") +
            geom_vline(xintercept=mCG.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCG.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCG.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCG.quant[4], color="blue", linetype="dashed")
mCHH.quant <- quantile(mCHH$coef, c(0.95,0.05,0.99,0.01))
o <- ggplot(mCHH, aes(coef)) + geom_density() + theme_bw() + xlab("GB Methylation CHH Feature Coefficients") +
            geom_vline(xintercept=mCHH.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCHH.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCHH.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCHH.quant[4], color="blue", linetype="dashed")
mCHG.quant <- quantile(mCHG$coef, c(0.95,0.05,0.99,0.01))
t <- ggplot(mCHG, aes(coef)) + geom_density() + theme_bw() + xlab("GB Methylation CHG Feature Coefficients") +
            geom_vline(xintercept=mCHG.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCHG.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCHG.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCHG.quant[4], color="blue", linetype="dashed")
G2.quant <- quantile(G2$coef, c(0.95,0.05,0.99,0.01))
s <- ggplot(G2, aes(coef)) + geom_density() + theme_bw() + xlab("SNP Feature Coefficients") +
            geom_vline(xintercept=G2.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=G2.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=G2.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=G2.quant[4], color="blue", linetype="dashed")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Figures/Based_On_Coef/rrBLUP")
ggsave("Density_dist_Coef_Expression_FT10_mean.pdf", plot=g, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_Coef_SNP_FT10_mean.pdf", plot=g2, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_Coef_GB_methy_FT10_mean.pdf", plot=p, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_Coef_GB_methyCG_FT10_mean.pdf", plot=l, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_Coef_GB_methyCHH_FT10_mean.pdf", plot=o, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_Coef_GB_methyCHG_FT10_mean.pdf", plot=t, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_max_per_gene_Coef_SNP_FT10_mean.pdf", plot=s, width=5, height=5, device="pdf", useDingbats=FALSE)

# Merge omics features
print("Merging omics datasets...")
G2vT <- merge(T, G2, by="genes") # G2 on y axis, T on x axis
G2vM <- merge(M, G2, by="genes")
G2vmCG <- merge(mCG, G2, by="genes")
G2vmCHH <- merge(mCHH, G2, by="genes")
G2vmCHG <- merge(mCHG, G2, by="genes")
TvM <- merge(M, T, by="genes")
TvmCG <- merge(mCG, T, by="genes")
TvmCHH <- merge(mCHH, T, by="genes")
TvmCHG <- merge(mCHG, T, by="genes")

# Sort and add feature rank column based on absolute coefficient
print("Ranking features based on absolute value of coefficient...")
G2vT <- G2vT[order(G2vT$abs_coef.x, decreasing=TRUE),] ; G2vT$rank.x <- seq(1, nrow(G2vT))
G2vM <- G2vM[order(G2vM$abs_coef.x, decreasing=TRUE),] ; G2vM$rank.x <- seq(1, nrow(G2vM))
G2vmCG <- G2vmCG[order(G2vmCG$abs_coef.x, decreasing=TRUE),] ; G2vmCG$rank.x <- seq(1, nrow(G2vmCG))
G2vmCHH <- G2vmCHH[order(G2vmCHH$abs_coef.x, decreasing=TRUE),] ; G2vmCHH$rank.x <- seq(1, nrow(G2vmCHH))
G2vmCHG <- G2vmCHG[order(G2vmCHG$abs_coef.x, decreasing=TRUE),] ; G2vmCHG$rank.x <- seq(1, nrow(G2vmCHG))
TvM <- TvM[order(TvM$abs_coef.x, decreasing=TRUE),] ; TvM$rank.x <- seq(1, nrow(TvM))
TvmCG <- TvmCG[order(TvmCG$abs_coef.x, decreasing=TRUE),] ; TvmCG$rank.x <- seq(1, nrow(TvmCG))
TvmCHH <- TvmCHH[order(TvmCHH$abs_coef.x, decreasing=TRUE),] ; TvmCHH$rank.x <- seq(1, nrow(TvmCHH))
TvmCHG <- TvmCHG[order(TvmCHG$abs_coef.x, decreasing=TRUE),] ; TvmCHG$rank.x <- seq(1, nrow(TvmCHG))

G2vT <- G2vT[order(G2vT$abs_coef.y, decreasing=TRUE),] ; G2vT$rank.y <- seq(1, nrow(G2vT))
G2vM <- G2vM[order(G2vM$abs_coef.y, decreasing=TRUE),] ; G2vM$rank.y <- seq(1, nrow(G2vM))
G2vmCG <- G2vmCG[order(G2vmCG$abs_coef.y, decreasing=TRUE),] ; G2vmCG$rank.y <- seq(1, nrow(G2vmCG))
G2vmCHH <- G2vmCHH[order(G2vmCHH$abs_coef.y, decreasing=TRUE),] ; G2vmCHH$rank.y <- seq(1, nrow(G2vmCHH))
G2vmCHG <- G2vmCHG[order(G2vmCHG$abs_coef.y, decreasing=TRUE),] ; G2vmCHG$rank.y <- seq(1, nrow(G2vmCHG))
TvM <- TvM[order(TvM$abs_coef.y, decreasing=TRUE),] ; TvM$rank.y <- seq(1, nrow(TvM))
TvmCG <- TvmCG[order(TvmCG$abs_coef.y, decreasing=TRUE),] ; TvmCG$rank.y <- seq(1, nrow(TvmCG))
TvmCHH <- TvmCHH[order(TvmCHH$abs_coef.y, decreasing=TRUE),] ; TvmCHH$rank.y <- seq(1, nrow(TvmCHH))
TvmCHG <- TvmCHG[order(TvmCHG$abs_coef.y, decreasing=TRUE),] ; TvmCHG$rank.y <- seq(1, nrow(TvmCHG))

# Add feature normalised rank (1-0 descending importance)
print("Normalizing feature rank...")
G2vT$rank_norm.x <- 1 - ((G2vT$rank.x - 1)/(nrow(G2vT) - 1))
G2vM$rank_norm.x <- 1 - ((G2vM$rank.x - 1)/(nrow(G2vM) - 1))
G2vmCG$rank_norm.x <- 1 - ((G2vmCG$rank.x - 1)/(nrow(G2vmCG) - 1))
G2vmCHH$rank_norm.x <- 1 - ((G2vmCHH$rank.x - 1)/(nrow(G2vmCHH) - 1))
G2vmCHG$rank_norm.x <- 1 - ((G2vmCHG$rank.x - 1)/(nrow(G2vmCHG) - 1))
TvM$rank_norm.x <- 1 - ((TvM$rank.x - 1)/(nrow(TvM) - 1))
TvmCG$rank_norm.x <- 1 - ((TvmCG$rank.x - 1)/(nrow(TvmCG) - 1))
TvmCHH$rank_norm.x <- 1 - ((TvmCHH$rank.x - 1)/(nrow(TvmCHH) - 1))
TvmCHG$rank_norm.x <- 1 - ((TvmCHG$rank.x - 1)/(nrow(TvmCHG) - 1))

G2vT$rank_norm.y <- 1 - ((G2vT$rank.y - 1)/(nrow(G2vT) - 1))
G2vM$rank_norm.y <- 1 - ((G2vM$rank.y - 1)/(nrow(G2vM) - 1))
G2vmCG$rank_norm.y <- 1 - ((G2vmCG$rank.y - 1)/(nrow(G2vmCG) - 1))
G2vmCHH$rank_norm.y <- 1 - ((G2vmCHH$rank.y - 1)/(nrow(G2vmCHH) - 1))
G2vmCHG$rank_norm.y <- 1 - ((G2vmCHG$rank.y - 1)/(nrow(G2vmCHG) - 1))
TvM$rank_norm.y <- 1 - ((TvM$rank.y - 1)/(nrow(TvM) - 1))
TvmCG$rank_norm.y <- 1 - ((TvmCG$rank.y - 1)/(nrow(TvmCG) - 1))
TvmCHH$rank_norm.y <- 1 - ((TvmCHH$rank.y - 1)/(nrow(TvmCHH) - 1))
TvmCHG$rank_norm.y <- 1 - ((TvmCHG$rank.y - 1)/(nrow(TvmCHG) - 1))

# Generate density plots (Adapted from Christina Azodi on GitHub ShiuLab/Manuscript_Code/2019_expression_GP/scripts/plot_impCOR_G.R)
get_density <- function(x, y, ...){ # density of points
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}

print("Calculating point density based on absolute value of feature coefficients...")
G2vT$density_abs_coef <- get_density(G2vT$abs_coef.x, G2vT$abs_coef.y, n=100)
G2vM$density_abs_coef <- get_density(G2vM$abs_coef.x, G2vM$abs_coef.y, n=100)
G2vmCG$density_abs_coef <- get_density(G2vmCG$abs_coef.x, G2vmCG$abs_coef.y, n=100)
G2vmCHH$density_abs_coef <- get_density(G2vmCHH$abs_coef.x, G2vmCHH$abs_coef.y, n=100)
G2vmCHG$density_abs_coef <- get_density(G2vmCHG$abs_coef.x, G2vmCHG$abs_coef.y, n=100)
TvM$density_abs_coef <- get_density(TvM$abs_coef.x, TvM$abs_coef.y, n=100)
TvmCG$density_abs_coef <- get_density(TvmCG$abs_coef.x, TvmCG$abs_coef.y, n=100)
TvmCHH$density_abs_coef <- get_density(TvmCHH$abs_coef.x, TvmCHH$abs_coef.y, n=100)
TvmCHG$density_abs_coef <- get_density(TvmCHG$abs_coef.x, TvmCHG$abs_coef.y, n=100)

print("Calculating point density based on original value of feature coefficients...")
G2vT$density_coef <- get_density(G2vT$coef.x, G2vT$coef.y, n=100)
G2vM$density_coef <- get_density(G2vM$coef.x, G2vM$coef.y, n=100)
G2vmCG$density_coef <- get_density(G2vmCG$coef.x, G2vmCG$coef.y, n=100)
G2vmCHH$density_coef <- get_density(G2vmCHH$coef.x, G2vmCHH$coef.y, n=100)
G2vmCHG$density_coef <- get_density(G2vmCHG$coef.x, G2vmCHG$coef.y, n=100)
TvM$density_coef <- get_density(TvM$coef.x, TvM$coef.y, n=100) 
TvmCG$density_coef <- get_density(TvmCG$coef.x, TvmCG$coef.y, n=100)
TvmCHH$density_coef <- get_density(TvmCHH$coef.x, TvmCHH$coef.y, n=100)
TvmCHG$density_coef <- get_density(TvmCHG$coef.x, TvmCHG$coef.y, n=100)

# Functions for plotting density scatterplot
theme_set(theme_bw(base_size=16)) # set plot theme
#Rho Unicode: U+1D746, UTF-8: F0 9D 9D 86
plot_one_percentile <- function(data, x, y, density, correlation, px, py, xlab, ylab, save_name){
    # Plot dashed lines for one percentile cutoff per omic type (2 total lines)
    p <- ggplot(data) + geom_point(aes(x=x, y=y, color=density), shape=16, alpha=0.75, size=1) + 
            scale_color_viridis() + theme_bw() + 
            geom_hline(yintercept=py, color="red", linetype="dashed") + 
            geom_vline(xintercept=px, color="red", linetype="dashed") +
            ggtitle(paste("rho = ", round(correlation, 3))) + xlab(xlab) + ylab(ylab) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p1 <- ggMarginal(p, type="density")
    return (p1)
}

plot_two_percentile <- function(data, x, y, density, correlation, px, py, px2, py2, xlab, ylab){
    # Plot dashed lines for two percentile cutoffs per omic type (4 total lines)
    p <- ggplot(data) + geom_point(aes(x=x, y=y, color=density), shape=16, alpha=0.75, size=1) + 
            scale_color_viridis() + theme_bw() +
            geom_hline(yintercept=py, color="red", linetype="dashed") + 
            geom_hline(yintercept=py2, color="blue", linetype="dashed") +
            geom_vline(xintercept=px, color="red", linetype="dashed") +
            geom_vline(xintercept=px2, color="blue", linetype="dashed") +
            ggtitle(paste("rho = ", round(correlation, 3))) + xlab(xlab) + ylab(ylab) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p1 <- ggMarginal(p, type="density")
    return (p1)
}

plot_density_coef <- function(data, x, y, density, px, py, px2, py2, xlab, ylab, save_name, abs=TRUE, nper=2){
    # data = merged dataframe between omics types
    # px, py, px2, py2 = corresponding percentile values for individual omics types (x, y)
    # save_name = plot file name
    if (abs==TRUE){ # based on absolute value of feature coefficient or importance score
        print("Calculating correlation between omics pair absolute value... ") ; print(deparse(substitute(data)))
        correlation <- cor(x, y, use="complete.obs", method ="spearman")
        print("Plotting... ")
        if (nper==2){
            p <- plot_one_percentile(data, x, y, density, correlation, px, py, xlab, ylab)
            ggsave(paste("Abs_", save_name, sep=""), plot=p, width=5, height=5, device="pdf", useDingbats=FALSE)
        } else if (nper==4){
            p2 <- plot_two_percentile(data, x, y, density, correlation, px, py, px2, py2, xlab, ylab)
            ggsave(paste("Abs_four_quant_", save_name, sep="", plot=p2, width=5, height=5, device="pdf", useDingbats=FALSE))
        }
    } else { # based on original value of feature coefficient or importance score
        print("Calculating correlation between omics pair original value... ") ; print(deparse(substitute(data)))
        correlation <- cor(x, y, use="complete.obs", method ="spearman")
        print("Plotting... ")
        if (nper==2){
            l <- plot_one_percentile(data, x, y, density, correlation, px, py, xlab, ylab)
            ggsave(paste(save_name, sep=""), plot=l, width=5, height=5, device="pdf", useDingbats=FALSE)
        } else if (nper==4){
            l2 <- plot_two_percentile(data, x, y, density, correlation, px, py, px2, py2, xlab, ylab)
            ggsave(paste("Four_quant_", save_name, sep=""), plot=l2, width=5, height=5, device="pdf", useDingbats=FALSE)
        }
    }
}

plot_density_rank <- function(data, x, y, density, px, py, px2, py2, xlab, ylab, save_name, norm=TRUE, nper=2){
    if (norm==TRUE){ # based on normalized rank of features
        print("Calculating correlation between omics pair normalized rank... ") ; print(deparse(substitute(data)))
        correlation <- cor(x, y, use="complete.obs", method ="spearman")
        print("Plotting... ") ; print(save_name)
        if (nper==2){
            o <- plot_one_percentile(data, x, y, density, correlation, px, py, xlab, ylab)
            ggsave(paste("Rank_norm_", save_name, sep=""), plot=o, width=5, height=5, device="pdf", useDingbats=FALSE)
        } else if (nper==4){
            o2 <- plot_two_percentile(data, x, y, density, correlation, px, py, px2, py2, xlab, ylab)
            ggsave(paste("Rank_norm_four_quant_", save_name, sep=""), plot=o2, width=5, height=5, device="pdf", useDingbats=FALSE)
        }
    } else { # based on rank of features
        print("Calculating correlation between omics pair rank... ") ; print(deparse(substitute(data)))
        correlation <- cor(x, y, use="complete.obs", method ="spearman")
        print("Plotting... ") ; print(save_name)
        if (nper==2){
            t <- plot_one_percentile(data, x, y, density, correlation, px, py, xlab, ylab)
            #t + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
            ggsave(paste("Rank_", save_name, sep=""), plot=t, width=5, height=5, device="pdf", useDingbats=FALSE)
        } else if (nper==4){
            t2 <- plot_two_percentile(data, x, y, density, correlation, px, py, px2, py2, xlab, ylab)
            #t2 + scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10")
            ggsave(paste("Rank_four_quant_", save_name, sep=""), plot=t2, width=5, height=5, device="pdf", useDingbats=FALSE)
        }
    }
}

################################### Absolute Value of Coef ##############################################
print("Generating density plots based on absolute value of feature coefficients...")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Figures/Based_On_Coef/rrBLUP")
print("...for 95th percentile")
T.per95 <- quantile(G2vT$abs_coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vT$abs_coef.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(G2vT$abs_coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vT$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vT, G2vT$abs_coef.x, G2vT$abs_coef.y, G2vT$density_abs_coef, T.per95, G2.per95, T.per05, G2.per05, "Gene Expression Feature Abs Coef ", "SNP Feature Abs Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
M.per95 <- quantile(G2vM$abs_coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vM$abs_coef.y, 0.95, na.rm=TRUE)
M.per05 <- quantile(G2vM$abs_coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vM$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vM, G2vM$abs_coef.x, G2vM$abs_coef.y, G2vM$density_abs_coef, M.per95, G2.per95, M.per05, G2.per05, "GB Methylation Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCG.per95 <- quantile(G2vmCG$abs_coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCG$abs_coef.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(G2vmCG$abs_coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCG$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCG, G2vmCG$abs_coef.x, G2vmCG$abs_coef.y, G2vmCG$density_abs_coef, mCG.per95, G2.per95, mCG.per05, G2.per05, "GB Methlyation CG Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCHH.per95 <- quantile(G2vmCHH$abs_coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHH$abs_coef.y, 0.95, na.rm=TRUE)
mCHH.per05 <- quantile(G2vmCHH$abs_coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHH$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCHH, G2vmCHH$abs_coef.x, G2vmCHH$abs_coef.y, G2vmCHH$density_abs_coef, mCHH.per95, G2.per95, mCHH.per05, G2.per05, "GB Methylation CHH Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCHG.per95 <- quantile(G2vmCHG$abs_coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHG$abs_coef.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(G2vmCHG$abs_coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHG$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCHG, G2vmCHG$abs_coef.x, G2vmCHG$abs_coef.y, G2vmCHG$density_abs_coef, mCHG.per95, G2.per95, mCHG.per05, G2.per05, "GB Methylation CHG Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
T.per95 <- quantile(TvM$abs_coef.x, 0.95, na.rm=TRUE); M.per95 <- quantile(TvM$abs_coef.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(TvM$abs_coef.x, 0.05, na.rm=TRUE); M.per05 <- quantile(TvM$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvM, TvM$abs_coef.x, TvM$abs_coef.y, TvM$density_abs_coef, M.per95, T.per95, M.per05, T.per05, "GB Methylation Feature Abs Coef ", "Gene Expression Feature Abs Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCG.per95 <- quantile(TvmCG$abs_coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCG$abs_coef.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(TvmCG$abs_coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCG$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCG, TvmCG$abs_coef.x, TvmCG$abs_coef.y, TvmCG$density_abs_coef, mCG.per95, T.per95, mCG.per05, T.per05, "GB Methlyation CG Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCHH.per95 <- quantile(TvmCHH$abs_coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHH$abs_coef.y, 0.95, na.rm=TRUE)
mCHH.per05 <- quantile(TvmCHH$abs_coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHH$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCHH, TvmCHH$abs_coef.x, TvmCHH$abs_coef.y, TvmCHH$density_abs_coef, mCHH.per95, T.per95, mCHH.per05, T.per05, "GB Methylation CHH Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
mCHG.per95 <- quantile(TvmCHG$abs_coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHG$abs_coef.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(TvmCHG$abs_coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHG$abs_coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCHG, TvmCHG$abs_coef.x, TvmCHG$abs_coef.y, TvmCHG$density_abs_coef, mCHG.per95, T.per95, mCHG.per05, T.per05, "GB Methylation CHG Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=TRUE, nper=2)
print("...for 99th percentile")
T.per99 <- quantile(G2vT$abs_coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vT$abs_coef.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(G2vT$abs_coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vT$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vT, G2vT$abs_coef.x, G2vT$abs_coef.y, G2vT$density_abs_coef, T.per99, G2.per99, T.per01, G2.per01, "Gene Expression Feature Abs Coef ", "SNP Feature Abs Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
M.per99 <- quantile(G2vM$abs_coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vM$abs_coef.y, 0.99, na.rm=TRUE)
M.per01 <- quantile(G2vM$abs_coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vM$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vM, G2vM$abs_coef.x, G2vM$abs_coef.y, G2vM$density_abs_coef, M.per99, G2.per99, M.per01, G2.per01, "GB Methylation Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCG.per99 <- quantile(G2vmCG$abs_coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCG$abs_coef.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(G2vmCG$abs_coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCG$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCG, G2vmCG$abs_coef.x, G2vmCG$abs_coef.y, G2vmCG$density_abs_coef, mCG.per99, G2.per99, mCG.per01, G2.per01, "GB Methlyation CG Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCHH.per99 <- quantile(G2vmCHH$abs_coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHH$abs_coef.y, 0.99, na.rm=TRUE)
mCHH.per01 <- quantile(G2vmCHH$abs_coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHH$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCHH, G2vmCHH$abs_coef.x, G2vmCHH$abs_coef.y, G2vmCHH$density_abs_coef, mCHH.per99, G2.per99, mCHH.per01, G2.per01, "GB Methylation CHH Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCHG.per99 <- quantile(G2vmCHG$abs_coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHG$abs_coef.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(G2vmCHG$abs_coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHG$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCHG, G2vmCHG$abs_coef.x, G2vmCHG$abs_coef.y, G2vmCHG$density_abs_coef, mCHG.per99, G2.per99, mCHG.per01, G2.per01, "GB Methylation CHG Feature Abs Coef ", "SNP Feature Abs Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
T.per99 <- quantile(TvM$abs_coef.x, 0.99, na.rm=TRUE); M.per99 <- quantile(TvM$abs_coef.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(TvM$abs_coef.x, 0.01, na.rm=TRUE); M.per01 <- quantile(TvM$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvM, TvM$abs_coef.x, TvM$abs_coef.y, TvM$density_abs_coef, M.per99, T.per99, M.per01, T.per01, "GB Methylation Feature Abs Coef ", "Gene Expression Feature Abs Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCG.per99 <- quantile(TvmCG$abs_coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCG$abs_coef.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(TvmCG$abs_coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCG$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCG, TvmCG$abs_coef.x, TvmCG$abs_coef.y, TvmCG$density_abs_coef, mCG.per99, T.per99, mCG.per01, T.per01, "GB Methlyation CG Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCHH.per99 <- quantile(TvmCHH$abs_coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHH$abs_coef.y, 0.99, na.rm=TRUE)
mCHH.per01 <- quantile(TvmCHH$abs_coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHH$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCHH, TvmCHH$abs_coef.x, TvmCHH$abs_coef.y, TvmCHH$density_abs_coef, mCHH.per99, T.per99, mCHH.per01, T.per01, "GB Methylation CHH Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)
mCHG.per99 <- quantile(TvmCHG$abs_coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHG$abs_coef.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(TvmCHG$abs_coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHG$abs_coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCHG, TvmCHG$abs_coef.x, TvmCHG$abs_coef.y, TvmCHG$density_abs_coef, mCHG.per99, T.per99, mCHG.per01, T.per01, "GB Methylation CHG Feature Abs Coef ",  "Gene Expression Feature Abs Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=TRUE, nper=2)

################################### Original Coef ##############################################
print("Generating density plots based on original value of feature coefficients...")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Figures/Based_On_Coef/rrBLUP")
print("...for 95th percentile")
T.per95 <- quantile(G2vT$coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vT$coef.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(G2vT$coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vT$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vT, G2vT$coef.x, G2vT$coef.y, G2vT$density_coef, T.per95, G2.per95, T.per05, G2.per05, "Gene Expression Feature Coef ", "SNP Feature Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vT, G2vT$coef.x, G2vT$coef.y, G2vT$density_coef, T.per95, G2.per95, T.per05, G2.per05, "Gene Expression Feature Coef ", "SNP Feature Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
M.per95 <- quantile(G2vM$coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vM$coef.y, 0.95, na.rm=TRUE)
M.per05 <- quantile(G2vM$coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vM$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vM, G2vM$coef.x, G2vM$coef.y, G2vM$density_coef, M.per95, G2.per95, M.per05, G2.per05, "GB Methylation Feature Coef ", "SNP Feature Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vM, G2vM$coef.x, G2vM$coef.y, G2vM$density_coef, M.per95, G2.per95, M.per05, G2.per05, "GB Methylation Feature Coef ", "SNP Feature Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCG.per95 <- quantile(G2vmCG$coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCG$coef.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(G2vmCG$coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCG$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCG, G2vmCG$coef.x, G2vmCG$coef.y, G2vmCG$density_coef, mCG.per95, G2.per95, mCG.per05, G2.per05, "GB Methlyation CG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCG, G2vmCG$coef.x, G2vmCG$coef.y, G2vmCG$density_coef, mCG.per95, G2.per95, mCG.per05, G2.per05, "GB Methlyation CG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCHH.per95 <- quantile(G2vmCHH$coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHH$coef.y, 0.95, na.rm=TRUE)
mCHH.per05 <- quantile(G2vmCHH$coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHH$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCHH, G2vmCHH$coef.x, G2vmCHH$coef.y, G2vmCHH$density_coef, mCHH.per95, G2.per95, mCHH.per05, G2.per05, "GB Methylation CHH Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCHH, G2vmCHH$coef.x, G2vmCHH$coef.y, G2vmCHH$density_coef, mCHH.per95, G2.per95, mCHH.per05, G2.per05, "GB Methylation CHH Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCHG.per95 <- quantile(G2vmCHG$coef.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHG$coef.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(G2vmCHG$coef.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHG$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(G2vmCHG, G2vmCHG$coef.x, G2vmCHG$coef.y, G2vmCHG$density_coef, mCHG.per95, G2.per95, mCHG.per05, G2.per05, "GB Methylation CHG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCHG, G2vmCHG$coef.x, G2vmCHG$coef.y, G2vmCHG$density_coef, mCHG.per95, G2.per95, mCHG.per05, G2.per05, "GB Methylation CHG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
T.per95 <- quantile(TvM$coef.x, 0.95, na.rm=TRUE); M.per95 <- quantile(TvM$coef.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(TvM$coef.x, 0.05, na.rm=TRUE); M.per05 <- quantile(TvM$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvM, TvM$coef.x, TvM$coef.y, TvM$density_coef, M.per95, T.per95, M.per05, T.per05, "GB Methylation Feature Coef ", "Gene Expression Feature Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(TvM, TvM$coef.x, TvM$coef.y, TvM$density_coef, M.per95, T.per95, M.per05, T.per05, "GB Methylation Feature Coef ", "Gene Expression Feature Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCG.per95 <- quantile(TvmCG$coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCG$coef.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(TvmCG$coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCG$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCG, TvmCG$coef.x, TvmCG$coef.y, TvmCG$density_coef, mCG.per95, T.per95, mCG.per05, T.per05, "GB Methlyation CG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCG, TvmCG$coef.x, TvmCG$coef.y, TvmCG$density_coef, mCG.per95, T.per95, mCG.per05, T.per05, "GB Methlyation CG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCHH.per95 <- quantile(TvmCHH$coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHH$coef.y, 0.95, na.rm=TRUE)
mCHH.per05 <- quantile(TvmCHH$coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHH$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCHH, TvmCHH$coef.x, TvmCHH$coef.y, TvmCHH$density_coef, mCHH.per95, T.per95, mCHH.per05, T.per05, "GB Methylation CHH Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCHH, TvmCHH$coef.x, TvmCHH$coef.y, TvmCHH$density_coef, mCHH.per95, T.per95, mCHH.per05, T.per05, "GB Methylation CHH Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
mCHG.per95 <- quantile(TvmCHG$coef.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHG$coef.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(TvmCHG$coef.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHG$coef.y, 0.05, na.rm=TRUE)
plot_density_coef(TvmCHG, TvmCHG$coef.x, TvmCHG$coef.y, TvmCHG$density_coef, mCHG.per95, T.per95, mCHG.per05, T.per05, "GB Methylation CHG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCHG, TvmCHG$coef.x, TvmCHG$coef.y, TvmCHG$density_coef, mCHG.per95, T.per95, mCHG.per05, T.per05, "GB Methylation CHG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_95.pdf", abs=FALSE, nper=4)
print("...for 99th percentile")
T.per99 <- quantile(G2vT$coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vT$coef.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(G2vT$coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vT$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vT, G2vT$coef.x, G2vT$coef.y, G2vT$density_coef, T.per99, G2.per99, T.per01, G2.per01, "Gene Expression Feature Coef ", "SNP Feature Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vT, G2vT$coef.x, G2vT$coef.y, G2vT$density_coef, T.per99, G2.per99, T.per01, G2.per01, "Gene Expression Feature Coef ", "SNP Feature Coef ",  "SNPvExpression_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
M.per99 <- quantile(G2vM$coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vM$coef.y, 0.99, na.rm=TRUE)
M.per01 <- quantile(G2vM$coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vM$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vM, G2vM$coef.x, G2vM$coef.y, G2vM$density_coef, M.per99, G2.per99, M.per01, G2.per01, "GB Methylation Feature Coef ", "SNP Feature Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vM, G2vM$coef.x, G2vM$coef.y, G2vM$density_coef, M.per99, G2.per99, M.per01, G2.per01, "GB Methylation Feature Coef ", "SNP Feature Coef ", "SNPvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCG.per99 <- quantile(G2vmCG$coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCG$coef.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(G2vmCG$coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCG$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCG, G2vmCG$coef.x, G2vmCG$coef.y, G2vmCG$density_coef, mCG.per99, G2.per99, mCG.per01, G2.per01, "GB Methlyation CG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCG, G2vmCG$coef.x, G2vmCG$coef.y, G2vmCG$density_coef, mCG.per99, G2.per99, mCG.per01, G2.per01, "GB Methlyation CG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCHH.per99 <- quantile(G2vmCHH$coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHH$coef.y, 0.99, na.rm=TRUE)
mCHH.per01 <- quantile(G2vmCHH$coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHH$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCHH, G2vmCHH$coef.x, G2vmCHH$coef.y, G2vmCHH$density_coef, mCHH.per99, G2.per99, mCHH.per01, G2.per01, "GB Methylation CHH Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCHH, G2vmCHH$coef.x, G2vmCHH$coef.y, G2vmCHH$density_coef, mCHH.per99, G2.per99, mCHH.per01, G2.per01, "GB Methylation CHH Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCHG.per99 <- quantile(G2vmCHG$coef.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHG$coef.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(G2vmCHG$coef.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHG$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(G2vmCHG, G2vmCHG$coef.x, G2vmCHG$coef.y, G2vmCHG$density_coef, mCHG.per99, G2.per99, mCHG.per01, G2.per01, "GB Methylation CHG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(G2vmCHG, G2vmCHG$coef.x, G2vmCHG$coef.y, G2vmCHG$density_coef, mCHG.per99, G2.per99, mCHG.per01, G2.per01, "GB Methylation CHG Feature Coef ", "SNP Feature Coef ", "SNPvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
T.per99 <- quantile(TvM$coef.x, 0.99, na.rm=TRUE); M.per99 <- quantile(TvM$coef.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(TvM$coef.x, 0.01, na.rm=TRUE); M.per01 <- quantile(TvM$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvM, TvM$coef.x, TvM$coef.y, TvM$density_coef, M.per99, T.per99, M.per01, T.per01, "GB Methylation Feature Coef ", "Gene Expression Feature Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(TvM, TvM$coef.x, TvM$coef.y, TvM$density_coef, M.per99, T.per99, M.per01, T.per01, "GB Methylation Feature Coef ", "Gene Expression Feature Coef", "ExpressionvGB_methy_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCG.per99 <- quantile(TvmCG$coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCG$coef.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(TvmCG$coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCG$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCG, TvmCG$coef.x, TvmCG$coef.y, TvmCG$density_coef, mCG.per99, T.per99, mCG.per01, T.per01, "GB Methlyation CG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCG, TvmCG$coef.x, TvmCG$coef.y, TvmCG$density_coef, mCG.per99, T.per99, mCG.per01, T.per01, "GB Methlyation CG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCHH.per99 <- quantile(TvmCHH$coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHH$coef.y, 0.99, na.rm=TRUE)
mCHH.per01 <- quantile(TvmCHH$coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHH$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCHH, TvmCHH$coef.x, TvmCHH$coef.y, TvmCHH$density_coef, mCHH.per99, T.per99, mCHH.per01, T.per01, "GB Methylation CHH Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCHH, TvmCHH$coef.x, TvmCHH$coef.y, TvmCHH$density_coef, mCHH.per99, T.per99, mCHH.per01, T.per01, "GB Methylation CHH Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHH_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
mCHG.per99 <- quantile(TvmCHG$coef.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHG$coef.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(TvmCHG$coef.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHG$coef.y, 0.01, na.rm=TRUE)
plot_density_coef(TvmCHG, TvmCHG$coef.x, TvmCHG$coef.y, TvmCHG$density_coef, mCHG.per99, T.per99, mCHG.per01, T.per01, "GB Methylation CHG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=2)
plot_density_coef(TvmCHG, TvmCHG$coef.x, TvmCHG$coef.y, TvmCHG$density_coef, mCHG.per99, T.per99, mCHG.per01, T.per01, "GB Methylation CHG Feature Coef ",  "Gene Expression Feature Coef ", "ExpressionvGB_methyCHG_Coef_FT10_mean_density_plot_99.pdf", abs=FALSE, nper=4)
