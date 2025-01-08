###################################################
# Density plots of Genomic Prediction Analysis of
# flowering time to compare top features (genomic, 
# transcriptomic, methylomic [gene body]) by rank
# of RF feature importance scores.
#
# Written by: Kenia Segura Abá
####################################################

library("dplyr")
library("ggplot2")
library("ggExtra")
library("data.table")
library("MASS")
library("viridis")

# Read in feature importance scores for genomic, transcriptomic, and methylomic data
print("Reading in files...")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Data/") # Path to files
T <- fread("Rank_RF_Expression_FT10_mean_imp.csv")
M <- fread("Rank_RF_GB_methy_FT10_mean_imp.csv")
G <- fread("Rank_RF_SNP_FT10_mean_imp.csv")
G2 <- fread("Rank_max_per_gene_RF_SNP_FT10_mean_imp.csv") # SNP w/max value represents the gene

# Rename G2 column names
G2 <- rename(G2, mean_imp = abs_imp) # importance scores are positive values

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

# Plot feature importance score distribution
T.quant <- quantile(T$mean_imp, c(0.95,0.05,0.99,0.01))
g <- ggplot(T, aes(mean_imp)) + geom_density() + theme_bw() + xlab("Expression Feature Importance Scores") +
            geom_vline(xintercept=T.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=T.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=T.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=T.quant[4], color="blue", linetype="dashed")
G.quant <- quantile(G$mean_imp, c(0.95,0.05,0.99,0.01))
g2 <- ggplot(G, aes(mean_imp)) + geom_density() + theme_bw() + xlab("SNP Feature Importance Scores") +
            geom_vline(xintercept=G.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=G.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=G.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=G.quant[4], color="blue", linetype="dashed")
M.quant <- quantile(M$mean_imp, c(0.95,0.05,0.99,0.01))
p <- ggplot(M, aes(mean_imp)) + geom_density() + theme_bw() + xlab("GB Methylation Feature Importance Scores") +
            geom_vline(xintercept=M.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=M.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=M.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=M.quant[4], color="blue", linetype="dashed")
mCG.quant <- quantile(mCG$mean_imp, c(0.95,0.05,0.99,0.01))
l <- ggplot(mCG, aes(mean_imp)) + geom_density() + theme_bw() + xlab("GB Methylation CG Feature Importance Scores") +
            geom_vline(xintercept=mCG.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCG.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCG.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCG.quant[4], color="blue", linetype="dashed")
mCHH.quant <- quantile(mCHH$mean_imp, c(0.95,0.05,0.99,0.01))
o <- ggplot(mCHH, aes(mean_imp)) + geom_density() + theme_bw() + xlab("GB Methylation CHH Feature Importance Scores") +
            geom_vline(xintercept=mCHH.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCHH.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCHH.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCHH.quant[4], color="blue", linetype="dashed")
mCHG.quant <- quantile(mCHG$mean_imp, c(0.95,0.05,0.99,0.01))
t <- ggplot(mCHG, aes(mean_imp)) + geom_density() + theme_bw() + xlab("GB Methylation CHG Feature Importance Scores") +
            geom_vline(xintercept=mCHG.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=mCHG.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=mCHG.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=mCHG.quant[4], color="blue", linetype="dashed")
G2.quant <- quantile(G2$mean_imp, c(0.95,0.05,0.99,0.01))
s <- ggplot(G2, aes(mean_imp)) + geom_density() + theme_bw() + xlab("SNP Feature Importance Scores") +
            geom_vline(xintercept=G2.quant[1], color="red", linetype="dashed") + 
            geom_vline(xintercept=G2.quant[2], color="red", linetype="dashed") +
            geom_vline(xintercept=G2.quant[3], color="blue", linetype="dashed") +
            geom_vline(xintercept=G2.quant[4], color="blue", linetype="dashed")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Figures/Based_On_Coef/RF")
ggsave("Density_dist_RF_Expression_FT10_mean_imp.pdf", plot=g, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_RF_SNP_FT10_mean_imp.pdf", plot=g2, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_RF_GB_methy_FT10_mean_imp.pdf", plot=p, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_RF_GB_methyCG_FT10_mean_imp.pdf", plot=l, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_RF_GB_methyCHH_FT10_mean_imp.pdf", plot=o, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_RF_GB_methyCHG_FT10_mean_imp.pdf", plot=t, width=5, height=5, device="pdf", useDingbats=FALSE)
ggsave("Density_dist_max_per_gene_RF_SNP_FT10_mean_imp.pdf", plot=s, width=5, height=5, device="pdf", useDingbats=FALSE)

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

# Sort and add feature rank column based on importance scores
print("Ranking features based on importance scores...")
G2vT <- G2vT[order(G2vT$mean_imp.x, decreasing=TRUE),] ; G2vT$rank.x <- seq(1, nrow(G2vT))
G2vM <- G2vM[order(G2vM$mean_imp.x, decreasing=TRUE),] ; G2vM$rank.x <- seq(1, nrow(G2vM))
G2vmCG <- G2vmCG[order(G2vmCG$mean_imp.x, decreasing=TRUE),] ; G2vmCG$rank.x <- seq(1, nrow(G2vmCG))
G2vmCHH <- G2vmCHH[order(G2vmCHH$mean_imp.x, decreasing=TRUE),] ; G2vmCHH$rank.x <- seq(1, nrow(G2vmCHH))
G2vmCHG <- G2vmCHG[order(G2vmCHG$mean_imp.x, decreasing=TRUE),] ; G2vmCHG$rank.x <- seq(1, nrow(G2vmCHG))
TvM <- TvM[order(TvM$mean_imp.x, decreasing=TRUE),] ; TvM$rank.x <- seq(1, nrow(TvM))
TvmCG <- TvmCG[order(TvmCG$mean_imp.x, decreasing=TRUE),] ; TvmCG$rank.x <- seq(1, nrow(TvmCG))
TvmCHH <- TvmCHH[order(TvmCHH$mean_imp.x, decreasing=TRUE),] ; TvmCHH$rank.x <- seq(1, nrow(TvmCHH))
TvmCHG <- TvmCHG[order(TvmCHG$mean_imp.x, decreasing=TRUE),] ; TvmCHG$rank.x <- seq(1, nrow(TvmCHG))

G2vT <- G2vT[order(G2vT$mean_imp.y, decreasing=TRUE),] ; G2vT$rank.y <- seq(1, nrow(G2vT))
G2vM <- G2vM[order(G2vM$mean_imp.y, decreasing=TRUE),] ; G2vM$rank.y <- seq(1, nrow(G2vM))
G2vmCG <- G2vmCG[order(G2vmCG$mean_imp.y, decreasing=TRUE),] ; G2vmCG$rank.y <- seq(1, nrow(G2vmCG))
G2vmCHH <- G2vmCHH[order(G2vmCHH$mean_imp.y, decreasing=TRUE),] ; G2vmCHH$rank.y <- seq(1, nrow(G2vmCHH))
G2vmCHG <- G2vmCHG[order(G2vmCHG$mean_imp.y, decreasing=TRUE),] ; G2vmCHG$rank.y <- seq(1, nrow(G2vmCHG))
TvM <- TvM[order(TvM$mean_imp.y, decreasing=TRUE),] ; TvM$rank.y <- seq(1, nrow(TvM))
TvmCG <- TvmCG[order(TvmCG$mean_imp.y, decreasing=TRUE),] ; TvmCG$rank.y <- seq(1, nrow(TvmCG))
TvmCHH <- TvmCHH[order(TvmCHH$mean_imp.y, decreasing=TRUE),] ; TvmCHH$rank.y <- seq(1, nrow(TvmCHH))
TvmCHG <- TvmCHG[order(TvmCHG$mean_imp.y, decreasing=TRUE),] ; TvmCHG$rank.y <- seq(1, nrow(TvmCHG))

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

print("Calculating point density based on importance scores...")
G2vT$density_mean_imp <- get_density(G2vT$mean_imp.x, G2vT$mean_imp.y, n=100)
G2vM$density_mean_imp <- get_density(G2vM$mean_imp.x, G2vM$mean_imp.y, n=100)
G2vmCG$density_mean_imp <- get_density(G2vmCG$mean_imp.x, G2vmCG$mean_imp.y, n=100)
#G2vmCHH$density_mean_imp <- get_density(G2vmCHH$mean_imp.x, G2vmCHH$mean_imp.y, n=100) # at least two quantiles are zero, so cannot get density
G2vmCHG$density_mean_imp <- get_density(G2vmCHG$mean_imp.x, G2vmCHG$mean_imp.y, n=100)
TvM$density_mean_imp <- get_density(TvM$mean_imp.x, TvM$mean_imp.y, n=100) 
TvmCG$density_mean_imp <- get_density(TvmCG$mean_imp.x, TvmCG$mean_imp.y, n=100)
#TvmCHH$density_mean_imp <- get_density(TvmCHH$mean_imp.x, TvmCHH$mean_imp.y, n=100)
TvmCHG$density_mean_imp <- get_density(TvmCHG$mean_imp.x, TvmCHG$mean_imp.y, n=100)

print("Calculating point density based on feature normalized rank...") # based on the normalized rank
G2vT$density_rank_norm <- get_density(G2vT$rank_norm.x, G2vT$rank_norm.y, n=100)
G2vM$density_rank_norm <- get_density(G2vM$rank_norm.x, G2vM$rank_norm.y, n=100)
G2vmCG$density_rank_norm <- get_density(G2vmCG$rank_norm.x, G2vmCG$rank_norm.y, n=100)
#G2vmCHH$density_rank_norm <- get_density(G2vmCHH$rank_norm.x, G2vmCHH$rank_norm.y, n=100)
G2vmCHG$density_rank_norm <- get_density(G2vmCHG$rank_norm.x, G2vmCHG$rank_norm.y, n=100)
TvM$density_rank_norm <- get_density(TvM$rank_norm.x, TvM$rank_norm.y, n=100)
TvmCG$density_rank_norm <- get_density(TvmCG$rank_norm.x, TvmCG$rank_norm.y, n=100)
#TvmCHH$density_rank_norm <- get_density(TvmCHH$rank_norm.x, TvmCHH$rank_norm.y, n=100)
TvmCHG$density_rank_norm <- get_density(TvmCHG$rank_norm.x, TvmCHG$rank_norm.y, n=100)

print("Calculating point density based on feature rank...") # based on original rank
G2vT$density_rank <- get_density(G2vT$rank.x, G2vT$rank.y, n=100)
G2vM$density_rank <- get_density(G2vM$rank.x, G2vM$rank.y, n=100)
G2vmCG$density_rank <- get_density(G2vmCG$rank.x, G2vmCG$rank.y, n=100)
#G2vmCHH$density_rank <- get_density(G2vmCHH$rank.x, G2vmCHH$rank.y, n=100)
G2vmCHG$density_rank <- get_density(G2vmCHG$rank.x, G2vmCHG$rank.y, n=100)
TvM$density_rank <- get_density(TvM$rank.x, TvM$rank.y, n=100)
TvmCG$density_rank <- get_density(TvmCG$rank.x, TvmCG$rank.y, n=100)
#TvmCHH$density_rank <- get_density(TvmCHH$rank.x, TvmCHH$rank.y, n=100)
TvmCHG$density_rank <- get_density(TvmCHG$rank.x, TvmCHG$rank.y, n=100)

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
            geom_hline(yintercept=py, color="blue", linetype="dashed") + 
            geom_hline(yintercept=py2, color="red", linetype="dashed") +
            geom_vline(xintercept=px, color="blue", linetype="dashed") +
            geom_vline(xintercept=px2, color="red", linetype="dashed") +
            ggtitle(paste("rho = ", round(correlation, 3))) + xlab(xlab) + ylab(ylab) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p1 <- ggMarginal(p, type="density")
    return (p1)
}

plot_density_imp <- function(data, x, y, density, px, py, px2, py2, xlab, ylab, save_name, nper=2){
    # data = merged dataframe between omics types
    # px, py, px2, py2 = corresponding percentile values for individual omics types (x, y)
    # save_name = plot file name
    # based on original value of feature coefficient or importance score
    print("Calculating correlation between omics pair feature importance scores... ") ; print(deparse(substitute(data)))
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

################################### Original Importance Scores ##############################################
print("Generating density plots based on importance scores...")
setwd("/mnt/home/seguraab/Shiu_Lab/Collabs/Multi_Omic/Figures/Based_On_Coef/RF")
T.per95 <- quantile(G2vT$mean_imp.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vT$mean_imp.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(G2vT$mean_imp.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vT$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(G2vT, G2vT$mean_imp.x, G2vT$mean_imp.y, G2vT$density_mean_imp, T.per95, G2.per95, T.per05, G2.per05, "Gene Expression Feature Importance Scores ", "SNP Feature Importance Scores ",  "SNPvExpression_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(G2vT, G2vT$mean_imp.x, G2vT$mean_imp.y, G2vT$density_mean_imp, T.per95, G2.per95, T.per05, G2.per05, "Gene Expression Feature Importance Scores ", "SNP Feature Importance Scores ",  "SNPvExpression_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
M.per95 <- quantile(G2vM$mean_imp.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vM$mean_imp.y, 0.95, na.rm=TRUE)
M.per05 <- quantile(G2vM$mean_imp.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vM$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(G2vM, G2vM$mean_imp.x, G2vM$mean_imp.y, G2vM$density_mean_imp, M.per95, G2.per95, M.per05, G2.per05, "GB Methylation Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methy_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(G2vM, G2vM$mean_imp.x, G2vM$mean_imp.y, G2vM$density_mean_imp, M.per95, G2.per95, M.per05, G2.per05, "GB Methylation Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methy_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
mCG.per95 <- quantile(G2vmCG$mean_imp.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCG$mean_imp.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(G2vmCG$mean_imp.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCG$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(G2vmCG, G2vmCG$mean_imp.x, G2vmCG$mean_imp.y, G2vmCG$density_mean_imp, mCG.per95, G2.per95, mCG.per05, G2.per05, "GB Methlyation CG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCG_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(G2vmCG, G2vmCG$mean_imp.x, G2vmCG$mean_imp.y, G2vmCG$density_mean_imp, mCG.per95, G2.per95, mCG.per05, G2.per05, "GB Methlyation CG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCG_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
#mCHH.per95 <- quantile(G2vmCHH$mean_imp.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHH$mean_imp.y, 0.95, na.rm=TRUE)
#mCHH.per05 <- quantile(G2vmCHH$mean_imp.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHH$mean_imp.y, 0.05, na.rm=TRUE)
#plot_density_imp(G2vmCHH, G2vmCHH$mean_imp.x, G2vmCHH$mean_imp.y, G2vmCHH$density_mean_imp, mCHH.per95, G2.per95, mCHH.per05, G2.per05, "GB Methylation CHH Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHH_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
#plot_density_imp(G2vmCHH, G2vmCHH$mean_imp.x, G2vmCHH$mean_imp.y, G2vmCHH$density_mean_imp, mCHH.per95, G2.per95, mCHH.per05, G2.per05, "GB Methylation CHH Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHH_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
mCHG.per95 <- quantile(G2vmCHG$mean_imp.x, 0.95, na.rm=TRUE); G2.per95 <- quantile(G2vmCHG$mean_imp.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(G2vmCHG$mean_imp.x, 0.05, na.rm=TRUE); G2.per05 <- quantile(G2vmCHG$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(G2vmCHG, G2vmCHG$mean_imp.x, G2vmCHG$mean_imp.y, G2vmCHG$density_mean_imp, mCHG.per95, G2.per95, mCHG.per05, G2.per05, "GB Methylation CHG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHG_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(G2vmCHG, G2vmCHG$mean_imp.x, G2vmCHG$mean_imp.y, G2vmCHG$density_mean_imp, mCHG.per95, G2.per95, mCHG.per05, G2.per05, "GB Methylation CHG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHG_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
T.per95 <- quantile(TvM$mean_imp.x, 0.95, na.rm=TRUE); M.per95 <- quantile(TvM$mean_imp.y, 0.95, na.rm=TRUE)
T.per05 <- quantile(TvM$mean_imp.x, 0.05, na.rm=TRUE); M.per05 <- quantile(TvM$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(TvM, TvM$mean_imp.x, TvM$mean_imp.y, TvM$density_mean_imp, M.per95, T.per95, M.per05, T.per05, "GB Methylation Feature Importance Scores ", "Gene Expression Feature Importance Scores", "ExpressionvGB_methy_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(TvM, TvM$mean_imp.x, TvM$mean_imp.y, TvM$density_mean_imp, M.per95, T.per95, M.per05, T.per05, "GB Methylation Feature Importance Scores ", "Gene Expression Feature Importance Scores", "ExpressionvGB_methy_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
mCG.per95 <- quantile(TvmCG$mean_imp.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCG$mean_imp.y, 0.95, na.rm=TRUE)
mCG.per05 <- quantile(TvmCG$mean_imp.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCG$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(TvmCG, TvmCG$mean_imp.x, TvmCG$mean_imp.y, TvmCG$density_mean_imp, mCG.per95, T.per95, mCG.per05, T.per05, "GB Methlyation CG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCG_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(TvmCG, TvmCG$mean_imp.x, TvmCG$mean_imp.y, TvmCG$density_mean_imp, mCG.per95, T.per95, mCG.per05, T.per05, "GB Methlyation CG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCG_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
#mCHH.per95 <- quantile(TvmCHH$mean_imp.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHH$mean_imp.y, 0.95, na.rm=TRUE)
#mCHH.per05 <- quantile(TvmCHH$mean_imp.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHH$mean_imp.y, 0.05, na.rm=TRUE)
#plot_density_imp(TvmCHH, TvmCHH$mean_imp.x, TvmCHH$mean_imp.y, TvmCHH$density_mean_imp, mCHH.per95, T.per95, mCHH.per05, T.per05, "GB Methylation CHH Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHH_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
#plot_density_imp(TvmCHH, TvmCHH$mean_imp.x, TvmCHH$mean_imp.y, TvmCHH$density_mean_imp, mCHH.per95, T.per95, mCHH.per05, T.per05, "GB Methylation CHH Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHH_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
mCHG.per95 <- quantile(TvmCHG$mean_imp.x, 0.95, na.rm=TRUE); T.per95 <- quantile(TvmCHG$mean_imp.y, 0.95, na.rm=TRUE)
mCHG.per05 <- quantile(TvmCHG$mean_imp.x, 0.05, na.rm=TRUE); T.per05 <- quantile(TvmCHG$mean_imp.y, 0.05, na.rm=TRUE)
plot_density_imp(TvmCHG, TvmCHG$mean_imp.x, TvmCHG$mean_imp.y, TvmCHG$density_mean_imp, mCHG.per95, T.per95, mCHG.per05, T.per05, "GB Methylation CHG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHG_RF_FT10_mean_imp_density_plot_95.pdf", nper=2)
plot_density_imp(TvmCHG, TvmCHG$mean_imp.x, TvmCHG$mean_imp.y, TvmCHG$density_mean_imp, mCHG.per95, T.per95, mCHG.per05, T.per05, "GB Methylation CHG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHG_RF_FT10_mean_imp_density_plot_95.pdf", nper=4)
print("...for 99th percentile")
T.per99 <- quantile(G2vT$mean_imp.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vT$mean_imp.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(G2vT$mean_imp.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vT$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(G2vT, G2vT$mean_imp.x, G2vT$mean_imp.y, G2vT$density_mean_imp, T.per99, G2.per99, T.per01, G2.per01, "Gene Expression Feature Importance Scores ", "SNP Feature Importance Scores ",  "SNPvExpression_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(G2vT, G2vT$mean_imp.x, G2vT$mean_imp.y, G2vT$density_mean_imp, T.per99, G2.per99, T.per01, G2.per01, "Gene Expression Feature Importance Scores ", "SNP Feature Importance Scores ",  "SNPvExpression_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
M.per99 <- quantile(G2vM$mean_imp.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vM$mean_imp.y, 0.99, na.rm=TRUE)
M.per01 <- quantile(G2vM$mean_imp.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vM$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(G2vM, G2vM$mean_imp.x, G2vM$mean_imp.y, G2vM$density_mean_imp, M.per99, G2.per99, M.per01, G2.per01, "GB Methylation Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methy_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(G2vM, G2vM$mean_imp.x, G2vM$mean_imp.y, G2vM$density_mean_imp, M.per99, G2.per99, M.per01, G2.per01, "GB Methylation Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methy_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
mCG.per99 <- quantile(G2vmCG$mean_imp.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCG$mean_imp.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(G2vmCG$mean_imp.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCG$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(G2vmCG, G2vmCG$mean_imp.x, G2vmCG$mean_imp.y, G2vmCG$density_mean_imp, mCG.per99, G2.per99, mCG.per01, G2.per01, "GB Methlyation CG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCG_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(G2vmCG, G2vmCG$mean_imp.x, G2vmCG$mean_imp.y, G2vmCG$density_mean_imp, mCG.per99, G2.per99, mCG.per01, G2.per01, "GB Methlyation CG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCG_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
#mCHH.per99 <- quantile(G2vmCHH$mean_imp.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHH$mean_imp.y, 0.99, na.rm=TRUE)
#mCHH.per01 <- quantile(G2vmCHH$mean_imp.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHH$mean_imp.y, 0.01, na.rm=TRUE)
#plot_density_imp(G2vmCHH, G2vmCHH$mean_imp.x, G2vmCHH$mean_imp.y, G2vmCHH$density_mean_imp, mCHH.per99, G2.per99, mCHH.per01, G2.per01, "GB Methylation CHH Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHH_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
#plot_density_imp(G2vmCHH, G2vmCHH$mean_imp.x, G2vmCHH$mean_imp.y, G2vmCHH$density_mean_imp, mCHH.per99, G2.per99, mCHH.per01, G2.per01, "GB Methylation CHH Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHH_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
mCHG.per99 <- quantile(G2vmCHG$mean_imp.x, 0.99, na.rm=TRUE); G2.per99 <- quantile(G2vmCHG$mean_imp.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(G2vmCHG$mean_imp.x, 0.01, na.rm=TRUE); G2.per01 <- quantile(G2vmCHG$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(G2vmCHG, G2vmCHG$mean_imp.x, G2vmCHG$mean_imp.y, G2vmCHG$density_mean_imp, mCHG.per99, G2.per99, mCHG.per01, G2.per01, "GB Methylation CHG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHG_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(G2vmCHG, G2vmCHG$mean_imp.x, G2vmCHG$mean_imp.y, G2vmCHG$density_mean_imp, mCHG.per99, G2.per99, mCHG.per01, G2.per01, "GB Methylation CHG Feature Importance Scores ", "SNP Feature Importance Scores ", "SNPvGB_methyCHG_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
T.per99 <- quantile(TvM$mean_imp.x, 0.99, na.rm=TRUE); M.per99 <- quantile(TvM$mean_imp.y, 0.99, na.rm=TRUE)
T.per01 <- quantile(TvM$mean_imp.x, 0.01, na.rm=TRUE); M.per01 <- quantile(TvM$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(TvM, TvM$mean_imp.x, TvM$mean_imp.y, TvM$density_mean_imp, M.per99, T.per99, M.per01, T.per01, "GB Methylation Feature Importance Scores ", "Gene Expression Feature Importance Scores", "ExpressionvGB_methy_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(TvM, TvM$mean_imp.x, TvM$mean_imp.y, TvM$density_mean_imp, M.per99, T.per99, M.per01, T.per01, "GB Methylation Feature Importance Scores ", "Gene Expression Feature Importance Scores", "ExpressionvGB_methy_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
mCG.per99 <- quantile(TvmCG$mean_imp.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCG$mean_imp.y, 0.99, na.rm=TRUE)
mCG.per01 <- quantile(TvmCG$mean_imp.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCG$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(TvmCG, TvmCG$mean_imp.x, TvmCG$mean_imp.y, TvmCG$density_mean_imp, mCG.per99, T.per99, mCG.per01, T.per01, "GB Methlyation CG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCG_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(TvmCG, TvmCG$mean_imp.x, TvmCG$mean_imp.y, TvmCG$density_mean_imp, mCG.per99, T.per99, mCG.per01, T.per01, "GB Methlyation CG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCG_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
#mCHH.per99 <- quantile(TvmCHH$mean_imp.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHH$mean_imp.y, 0.99, na.rm=TRUE)
#mCHH.per01 <- quantile(TvmCHH$mean_imp.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHH$mean_imp.y, 0.01, na.rm=TRUE)
#plot_density_imp(TvmCHH, TvmCHH$mean_imp.x, TvmCHH$mean_imp.y, TvmCHH$density_mean_imp, mCHH.per99, T.per99, mCHH.per01, T.per01, "GB Methylation CHH Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHH_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
#plot_density_imp(TvmCHH, TvmCHH$mean_imp.x, TvmCHH$mean_imp.y, TvmCHH$density_mean_imp, mCHH.per99, T.per99, mCHH.per01, T.per01, "GB Methylation CHH Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHH_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
mCHG.per99 <- quantile(TvmCHG$mean_imp.x, 0.99, na.rm=TRUE); T.per99 <- quantile(TvmCHG$mean_imp.y, 0.99, na.rm=TRUE)
mCHG.per01 <- quantile(TvmCHG$mean_imp.x, 0.01, na.rm=TRUE); T.per01 <- quantile(TvmCHG$mean_imp.y, 0.01, na.rm=TRUE)
plot_density_imp(TvmCHG, TvmCHG$mean_imp.x, TvmCHG$mean_imp.y, TvmCHG$density_mean_imp, mCHG.per99, T.per99, mCHG.per01, T.per01, "GB Methylation CHG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHG_RF_FT10_mean_imp_density_plot_99.pdf", nper=2)
plot_density_imp(TvmCHG, TvmCHG$mean_imp.x, TvmCHG$mean_imp.y, TvmCHG$density_mean_imp, mCHG.per99, T.per99, mCHG.per01, T.per01, "GB Methylation CHG Feature Importance Scores ",  "Gene Expression Feature Importance Scores ", "ExpressionvGB_methyCHG_RF_FT10_mean_imp_density_plot_99.pdf", nper=4)
