################################################################################
# FIGURE 1
################################################################################
rm(list=ls())
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

setwd("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

# Isolate growth condition labels
cond <- c("YPACETATE", "YPD14", "YPD40", "YPD42", "YPD6AU", "YPDANISO10", 
        "YPDANISO20", "YPDANISO50", "YPDBENOMYL200", "YPDBENOMYL500", "YPDCAFEIN40", 
        "YPDCAFEIN50", "YPDCHX05", "YPDCHX1", "YPDCUSO410MM", "YPDDMSO", "YPDETOH", 
        "YPDFLUCONAZOLE", "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPDHU", "YPDKCL2M", 
        "YPDLICL250MM", "YPDMV", "YPDNACL15M", "YPDNACL1M", "YPDNYSTATIN", "YPDSDS", 
        "YPDSODIUMMETAARSENITE", "YPETHANOL", "YPGALACTOSE", "YPRIBOSE", "YPGLYCEROL", 
        "YPXYLOSE", "YPSORBITOL")
new_cond <- c("YP Acetate 2%", "YPD 14ºC", "YPD 40ºC", "YPD 42ºC", "YPD 6-Azauracile 600 µg/ml",
        "YPD Anisomycin 10 µg/ml", "YPD Anisomycin 20 µg/ml", "YPD Anisomycin 50 µg/ml",
        "YPD Benomyl 200 µg/ml", "YPD Benomyl 500 µg/ml", "YPD Caffeine 40 mM", "YPD Caffeine 50 mM",
        "YPD Cycloheximide 0.5 µg/ml", "YPD Cycloheximide 1 µg/ml", "YPD CuSO4 10 mM", "YPD DMSO 6%",
        "YPD Ethanol 15%", "YPD Fluconazole 20 µg/ml", "YPD Formamide 4%", "YPD Formamide 5%",
        "YPD Hydroxyurea 30 mg/ml", "YPD KCL 2 M", "YPD LiCl 250 mM", "YPD Methylviologen 20 mM",
        "YPD NaCl 1.5 M", "YPD NaCl 1 M", "YPD Nystatin 10 µg/ml", "YPD SDS 0.2%", 
        "YPD Sodium metaarsenite 2.5 mM", "YP Ethanol 2%", "YP Galactose 2%", "YP Ribose 2%",
        "YP Glycerol 2%", "YP Xylose 2%", "YP Sorbitol 2%")
conds <- as.data.frame(cbind(cond, new_cond))

pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # fitness data
pCorEnvs <- cor(as.matrix(pheno), method="pearson")
write.csv(pCorEnvs, "Data/Peter_2018/pheno_corr_envs.csv", quote=F, row.names=T)
pCorIso <- cor(as.matrix(t(pheno)), method="pearson")
write.csv(pCorIso, "Data/Peter_2018/pheno_corr_isolates.csv", quote=F, row.names=T)

pheno <- reshape2::melt(pheno) # pivot longer
pheno <- left_join(conds, pheno, by=c("cond"="variable")) # add condition labels

################################## FIGURE 1A ##################################
# Heatmap of pCorEnvs
rdbu_r <- rev(brewer.pal(n=9, "RdBu")) # reversed color palette
col_fun = colorRamp2(seq(-1,1,.25), rdbu_r) #same as cm
#cm = ColorMapping(name="PCC", colors=rdbu_r, levels=seq(-1,1,.2)) # color mapping
pdf("Scripts/Data_Vis/Section_1/pheno_corr_envs_v4.pdf", height=9, width=9)
hm3 <- ComplexHeatmap::Heatmap(as.matrix(pCorEnvs),
        col=col_fun, show_row_dend=F, show_column_dend=F,
        heatmap_legend_param=list(title="PCC", at=seq(-1,1,.25), color_bar="continuous"))
draw(hm3)
dev.off()

# clusters from pCorEnvs heatmap
hr <- hclust(dist(as.matrix(pCorEnvs), method="euclidean"), method="complete")
cluster_df <- as.data.frame(colnames(pCorEnvs))
for (k in 10:20){
        clusters <- as.data.frame(cutree(hr, k=k))
        cluster_df <- merge(cluster_df, clusters, by.x="colnames(pCorEnvs)", by.y="row.names")
}
colnames(cluster_df) <- c("Envs", paste("k=",10:20, sep="")) # k=19 looks the best
write.csv(cluster_df, "Scripts/Data_Vis/Section_1/pheno_corr_envs_v3_clusters.csv", quote=F)

# for k=19, calculate the minimum PCC for the environments in each cluster
for (grp in 1:19){
        envs <- cluster_df[which(cluster_df["k=19"]==grp), "Envs"]
        if (length(envs) > 1) print(paste(grp, ";", min(pCorEnvs[envs, envs]))) # minimum PCC
}
# [1] "2 ; 0.975362803318604"
# [1] "3 ; 0.518246302099321"
# [1] "4 ; 0.58645803304682"
# [1] "7 ; 0.784820502498895"
# [1] "8 ; 0.485669214866689"
# [1] "9 ; 0.58347711984573"
# [1] "17 ; 0.730897529069503"
# [1] "18 ; 0.752072512382595"

# Binned pCorEnv heatmap
pCorEnvs_binned <- apply(pCorEnvs, 2, cut, seq(-1,1,.2)) # bins on each column
rownames(pCorEnvs_binned) <- rownames(pCorEnvs) # set rownames
pCorEnvs_binned_melt <- reshape2::melt(pCorEnvs_binned) # pivot longer
pCorEnvs_melt <- reshape2::melt(as.matrix(pCorEnvs)) # pivot cor data longer
pCorEnvs_binned <- merge(pCorEnvs_melt, pCorEnvs_binned_melt, by=c("Var1", "Var2")) # add bin info to cor data
pCorEnvs_binned_av <- pCorEnvs_binned %>% group_by(Var1, Var2, value.y) %>% summarize(bin_mean = mean(value.x))  # bin average
pCorEnvs_binned_av2 <- reshape2::dcast(pCorEnvs_binned_av, formula=Var1~Var2, value.var="bin_mean") # pivot binned cor data wider
rownames(pCorEnvs_binned_av2) <- pCorEnvs_binned_av2$Var1 # set rownames
pCorEnvs_binned_av2 <- pCorEnvs_binned_av2[-1] # drop Var1

# draw heatmap
rdbu_r <- rev(brewer.pal(n=11, "RdBu")) # reversed color palette
col_fun = colorRamp2(seq(-1,1,.2), rdbu_r) #same as cm
pdf("Scripts/Data_Vis/pheno_corr_envs_binned.pdf", height=9, width=9)
hm3 <- ComplexHeatmap::Heatmap(as.matrix(pCorEnvs_binned_av2),
        col=col_fun, heatmap_legend_param=list(title="PCC", at=seq(-1,1,.2),
                color_bar="discrete"))
draw(hm3)
dev.off()
# the binned and original heatmaps look the same

######## Relate environment clusters to the narrow-sense heritabilities ########
str(hm3)
data <- hm3@matrix
env_order <- c("YPDCHX05", "YPDCHX1", "YPDANISO50", "YPDANISO10", "YPDANISO20",
             "YPDDMSO", "YPDMV", "YPDSDS", "YPD40", "YPD42", "YPDKCL2M",
             "YPDCAFEIN40", "YPDCAFEIN50", "YPDBENOMYL200", "YPDBENOMYL500",
             "YPDETOH", "YPDNYSTATIN", "YPACETATE", "YPXYLOSE", "YPRIBOSE",
             "YPSORBITOL", "YPGLYCEROL", "YPETHANOL", "YPGALACTOSE",
             "YPDLICL250MM", "YPDNACL15M", "YPDNACL1M", "YPDFORMAMIDE4",
             "YPDFORMAMIDE5", "YPDHU", "YPD14", "YPDFLUCONAZOLE",
             "YPDSODIUMMETAARSENITE", "YPD6AU", "YPDCUSO410MM")
clustered_data <- data[env_order, env_order]

h2 <- read.csv("Data/Peter_2018/Heritability_h2_H2_sommer.csv", row.names=1) # heritability data
rownames(h2) <- h2$Condition
h2 <- h2[, c("h2", "h2_SE")]
h2 <- h2[env_order, ] # reorder to match heatmap

# plot the heritabilities according to the env order of Figure 1A
rdylbu <- rev(brewer.pal(n=11, "RdYlBu")) # color palette
col_fun = colorRamp2(seq(0,1,.1), rdylbu)
h2_matrix <- as.matrix(h2$h2)
rownames(h2_matrix) <- rownames(h2)
pdf("Scripts/Data_Vis/Section_1/heritiabilities_ordered_by_Fig1A.pdf")
ht_h2 <- Heatmap(h2_matrix, name = "h²",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(title = "Heritability"),
  column_title = NULL,
  row_names_gp = gpar(fontsize = 8)
)

draw(ht_h2)
dev.off()

########################## ADDITIONAL VISUALIZATIONS ##########################
# Fitness distributions histogram for each column
ggplot(pheno, aes(x=value)) + theme_bw(8) + geom_histogram(binwidth=0.01) + xlab("Fitness") +
        theme(axis.text.x=element_text(color="black",  size=7)) +
        theme(axis.text.y=element_text(color="black", size=7)) +
        facet_wrap(~new_cond, scales="free_x")
ggsave("Scripts/Data_Vis/Section_3/fitness_distributions.pdf", width=8.5, height=11, device="pdf", useDingbats=FALSE)

# Violin plot of fitness distributions in each environment
ggplot(pheno, aes(x=reorder(new_cond, -value), y=value)) + theme_bw(8) +
        geom_violin(color="#5B9BD5") + geom_boxplot(width=0.1) + ylab("Fitness") +
        theme(axis.text.x=element_text(color="black",  size=9, face="bold", angle=55, hjust=1, )) +
        theme(axis.text.y=element_text(color="black", size=9, face="bold"))
ggsave("Scripts/Data_Vis/fitness_violinplot.pdf", width=14, height=4, device="pdf", useDingbats=FALSE)

# Fitness variance
pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # fitness data
stats <- do.call(cbind, lapply(pheno, summary))
stats <- rbind(stats, apply(pheno, 2, var))
rownames(stats)[7] <- "Var."
write.csv(t(stats), "Scripts/Data_Vis/Section_3/fitness_stats.csv", quote=F)

# Fitness correlations among replicates (ensure data quality)
reps <- read.csv("Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h.csv")
reps <- reps[!is.na(reps$Growth.ratio),]
length(unique(reps$Strain)) # only have replicate info for 382 strains :(
stats <- reps %>% group_by(Condition, Strain) %>% 
        dplyr::summarize(min=min(Growth.ratio), Q1=quantile(Growth.ratio, 0.25), 
                median=median(Growth.ratio), mean=mean(Growth.ratio), 
                sd=sd(Growth.ratio), Q3=quantile(Growth.ratio, 0.75), 
                max=max(Growth.ratio))
IDs <- read.csv("Data/Peter_2018/IDcorrespondance.txt", sep="\t")
stats <- right_join(IDs, stats, by=c("YJS_name"="Strain"))
pheno <- read.csv("Data/Peter_2018/pheno.csv")
stats <- stats[stats$Standardized_name %in% pheno$ID, ] # keep only diploids
length(unique(stats$Standardized_name)) # only have replicate info for 124 diploids
write.csv(stats, "Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h_only_diploids.csv", quote=F, row.names=F)
summary(stats[4:10]) # interested in the sd, these should be small
ggplot() + geom_point(aes(x=sd, y=mean, group=Condition, alpha=0.1), stats)
ggsave("Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h_only_diploids_stdev.pdf")

