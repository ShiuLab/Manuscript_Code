library(ggplot2)
library(reshape2)
library(data.table)

##########################
## Grid Search Analysis ##
##########################

# Prep data
gs <- fread('20190104_14_gridsearch_merged.csv')
gs <- data.frame(gs)
algs <- unique(gs$alg)
species <- unique(gs$species)
gs$mse_valid <- abs(gs$mse_valid)
gs[gs==''] <- NA
gs[, c('Epochs', 'V1', 'Train_Loss', 'Valid_PCC', 'kernel')] <- NULL

# Calculate Variance in MSE across grid search for each species/trait
gs$paramspace <- paste(gs$ActFun, gs$Arch, gs$C, gs$L2, gs$LearnRate, gs$degree, gs$dropout, gs$gamma, gs$learning_rate, gs$max_depth, gs$max_features, sep=':')
gs_mean <- aggregate(gs$mse_valid, gs[,c('species','alg','paramspace')], FUN=mean)
gs_var <- aggregate(gs$mse_valid, gs[,c('species','alg')], FUN=var)
reo_alg <- function(x) { factor(x, levels = c('SVM','SVMpoly','SVMrbf','RF','GB','ann'))}
ggplot(gs_var, aes(x=reo_alg(alg), y=log10(x))) + geom_violin() + theme_bw()+ geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = species))
ggsave('00_Figures/190305_GridS_varianceMSE.pdf', width=3, height = 3, device='pdf', useDingbats=FALSE)  # 1 column
gs_var_median <- aggregate(gs_var$x, list(gs_var$alg), FUN=median)

# ANOVA analysis
gs_res <- data.frame(sp = character(), alg = character(), p = character(), stringsAsFactors=FALSE)
for(al in algs){
  gs_al <- subset(gs, alg == al)
  gs_al <- gs_al[, colSums(is.na(gs_al)) != nrow(gs_al)]
  for(sp in species){
    gs_al_sp <- subset(gs_al, species == sp)
    gs_al_sp <- gs_al_sp[, !(names(gs_al_sp) %in% c('alg', 'species'))]
    if(dim(gs_al_sp)[1] >=9){
      aov <- summary(aov(mse_valid ~ ., gs_al_sp))
      p_len <- length(aov[[1]][["Pr(>F)"]]) - 1
      p_str <- paste(aov[[1]][["Pr(>F)"]][1:p_len], collapse=',')
      gs_res <- rbind(gs_res, list(sp=sp, alg=al, p = p_str), stringsAsFactors=F)
    }
  }
}
write.csv(gs_res, 'GridSearch_stats.csv', quote=FALSE)


#######################
## Feature Selection ##
#######################

#### Test RF, BayesA, and Elastic net as feature selection algorithms for height in maize

# Format results
lr <- read.table('FS_results_lin.txt', sep=',', header=T)
lr$sp <- as.character(lapply(strsplit(as.character(lr$species), split="_"),head, n=1))
lr <- subset(lr, model == 'rrBLUP')
lr <- lr[c('sp','model','feat_method','feat_num','accuracy_PCC')]
names(lr) <- c('sp','Alg','FSm','FSn','PCC')

ml <- read.table('FS_results_ml.txt', sep=',', header=T)
ml$sp <- as.character(lapply(strsplit(as.character(ml$species), split="_"),head, n=1))
ml$Tag <- gsub("fs_HT_", "", ml$Tag) 
ml$Tag <- gsub("_", "", ml$Tag) 
ml$Tag <- gsub("[0-9]", "", ml$Tag) 
ml <- ml[c('sp','Alg','Tag','FeatureNum','PCC_ho')]
names(ml) <- c('sp','Alg','FSm','FSn','PCC')

ann <- read.table('FS_results_ann_v2.txt', sep=',', header=T)
ann$sp <- as.character(lapply(strsplit(as.character(ann$species), split="_"),head, n=1))
ann <- ann[c('sp','NumFeat','PCC_test')]
names(ann) <- c('sp','FSn','PCC')
ann$Alg <- 'ANN'
ann$FSm <- 'na'

# Stats
maize <- rbind(lr, ml, ann)
maize <- subset(maize, sp=='maize')
summary(aov(PCC ~ Alg * FSm *FSn, maize))

# Background expectation for FS overlap 
d <-fread('~/03_GenomicSelection/maize_WDP_Hirsch/01_Data/geno_withDups.csv', sep=',')
m <- names(d)
m <- m[2:length(m)]
for(i in seq(1, 1000)){
  tmp <- sample(m, 8000, replace=F)
  sn <- paste('rand_',i,'.txt', sep='')
  write.table(tmp, file = sn, quote=F, row.names=F, col.names=F)
}



#### Using RF for feature selection, test rrBLUP, SVR, RF, GTB, and ANN

# Format results
lr <- read.table('FS_results_lin.txt', sep=',', header=T)
lr$sp <- as.character(lapply(strsplit(as.character(lr$species), split="_"),head, n=1))
lr <- subset(lr, lr$feat_method %in% c('RF','none','pass'))
lr <- lr[c('sp','model','cv_Num','feat_num','accuracy_PCC')]
names(lr) <- c('sp','Alg','cv','FSn','PCC')
lr$cv <- gsub("[^0-9]", "", lr$cv) 
lr <- subset(lr, Alg == 'rrBLUP')

ml <- read.table('FS_results_ml_noALL.txt', sep='\t', header=T)
ml$sp <- as.character(lapply(strsplit(as.character(ml$species), split="_"),head, n=1))
ml$cv <- as.character(lapply(strsplit(as.character(ml$Tag), split="_"),tail, n=1))
ml <- ml[!grepl('bayesa', ml$Tag),]
ml <- ml[!grepl('EN_', ml$Tag),]
ml <- ml[c('sp','cv','Alg','Tag','FeatureNum','PCC_ho')]
ml <- subset(ml, Alg %in% c('SVM', 'GB', 'RF')) # Removing poly and rbf SVM from experiment!!!
ml$Tag <- gsub('HT_', '', ml$Tag)
ml$Tag <- gsub('fs_HT_', '', ml$Tag)
ml$FS <- as.character(lapply(strsplit(as.character(ml$Tag), split="_"),head, n=1))
ml$FS <- NULL
names(ml) <- c('sp','cv','Alg','Tag','FSn','PCC')
ml$Tag <- NULL
ml <- rbind(ml, ml2)

ann <- read.table('FS_results_ann_v2.txt', sep = ',', header=T)
ann$sp <- as.character(lapply(strsplit(as.character(ann$species), split="_"),head, n=1))
ann$cv <- as.character(lapply(strsplit(as.character(ann$Holdout), split="/"),tail, n=1))
ann <- ann[c('sp','cv','model','FeatSel','NumFeat','PCC_test')]
ann$FeatSel <- as.character(lapply(strsplit(as.character(ann$FeatSel), split="/"),tail, n=1))
ann$FeatSel <- gsub("[a-z]*_HT_", '', ann$FeatSel)
ann$FS <- gsub("_[0-9]*", '', ann$FeatSel)
ann <- ann[c('sp', 'cv','model', 'NumFeat', 'PCC_test')]
names(ann) <- c('sp','cv', 'Alg', 'FSn', 'PCC')
ann <- rbind(ann, ann2)
ann$cv <- gsub("[^0-9]", "", ann$cv) 


#Stats
fs <- rbind(lr, ml, ann)
stats_results <- data.frame(alg = character(), sp= character(), FSn = character, num_all= character(),pval=character(), stringsAsFactors = FALSE)
for(al in unique(fs$Alg)){
  for(s in unique(fs$sp)){
    if(s %in% c('soy','spruce')){
      n <- 4000
    }
    else{
      n <- 8000
    }
    test <- subset(fs, Alg==al & sp==s & FSn %in% c('all',n))
    num_all <- dim(subset(test, FSn=='all'))[1]
    num_fs <- dim(subset(test, FSn !='all'))[1]
    if((num_all > 4) & (num_fs > 4)){
      resl <- wilcox.test(PCC~FSn, test, alternative = 'greater')
      stats_results <- rbind(stats_results, 
                             list(alg = al, species=s, FSn=n, pval = resl$p.value, num_all=num_all), stringsAsFactors=F)
    }}}

stats_results$q <- (p.adjust(stats_results$pval, "BH"))

# Stdev between replicates for each algorithm
stdev <- aggregate(PCC ~ sp+Alg+FSn, fs, FUN=sd)
summary(aov(PCC~Alg*FSn, subset(stdev, FSn !='all')))
wilcox.test(PCC~Alg, subset(stdev, FSn !='all'))

##########################
## Algorithm Comparison ##
##########################
setwd('/Volumes/azodichr/03_GenomicSelection/06_otherTraits/')
bglr <- read.csv('bglr/BGLR_RESULTS.csv', sep=',')
bglr <- bglr[,c('model', 'tag', 'y', 'PCC')]
rrb <- read.csv('rrblup/rrBLUP_RESULTS.csv')
rrb <- rrb[,c('model', 'tag', 'y', 'PCC')]
ml <- read.csv('ML/RESULTS_reg.txt', sep='\t')
ml <- ml[,c('Alg', 'Tag', 'Y', 'PCC_ho')]
ml$Tag <- as.character(lapply(strsplit(as.character(ml$Tag), split="_"),head, n=1))
names(ml) <- c('model', 'tag', 'y', 'PCC')

# ANN Feature Selection Results
annFS <- read.csv('mlp_use/RESULTS.txt', sep='\t')
annFS <- annFS[,c('DFs', 'Trait', 'PCC_test')]
annFS$DFs <- as.character(lapply(strsplit(as.character(annFS$DFs), split="03_GenomicSelection/"),tail, n=1))
annFS$DFs <- as.character(lapply(strsplit(as.character(annFS$DFs), split="_"),head, n=1))
names(annFS) <- c('tag', 'y', 'PCC')
annFS$model <- 'ANN_FS'

# Ensemble Results
en <- read.csv('ensem/RESULTS_ens.txt', sep='\t')
en$model <- 'EN5'
en$model[en$NumEnsemble > 5] = 'EN11'
en <- en[,c('Species','Y','PCC','model')]
names(en) <- c('tag', 'y', 'PCC', 'model')

# Seeded starting weight ANN Results
lrann <- read.csv('mlp_ens_2/fs_8000/02_UseInManuscript/RESULTS.txt', sep='\t')
lrann$sp <- as.character(lapply(strsplit(as.character(lrann$Tag), split="_"),head, n=1))
lrann$seed <- as.character(lapply(strsplit(as.character(lrann$Tag), split="_"),tail, n=1))
lrann$model <- paste('ANN_', lrann$seed, '_8k', sep='')
lrann <- aggregate(PCC_test ~ sp+Trait+model+Holdout, lrann, FUN=mean)
lrann <- lrann[,c('model','sp', 'Trait', 'PCC_test')]
names(lrann) <- c('model','tag', 'y', 'PCC')
lrann$tag <- gsub('sorgh', 'sorg', lrann$tag)

# Calulate mean PCC and % change from the best for each model
d <- rbind(bglr, rrb, ml, annFS, en, lrann)
count <- aggregate(PCC~tag+model+y, d, FUN=length)

keep = c('rrBLUP', 'BRR', 'BayesA', 'BayesB', 'BL', 'SVM','SVMpoly','SVMrbf', 'RF', 'GB', 'ANN_FS', 'EN5', 'EN11', 'ANN_rrB_8k','ANN_BayesB_8k', 'ANN_BL_8k','ANN_RF_8k', 'ANN_rrB_50p','ANN_BayesB_50p', 'ANN_BL_50p','ANN_RF_50p')

d_mean <- aggregate(PCC~tag+model+y, d, FUN=mean)
d_mean <- subset(d_mean, model %in% keep)
d_best <- aggregate(PCC~tag+y, d_mean, function(x) x[which.max(abs(x))])
d_mean <- merge(d_mean, d_best, by=c('tag', 'y'))
names(d_mean) <- c('species', 'y', 'alg', 'r', 'r_best')

d_mean$species <- gsub('sorgh', 'sorg', d_mean$species)
d_mean$pchange <- (d_mean$r_best-d_mean$r)/d_mean$r_best
d_mean$pchange_cut <- d_mean$pchange
d_mean$pchange_cut[d_mean$pchange_cut>0.2] <- 0.2

reo_meth <- function(x) { factor(x, levels = keep)}

# Plot heatmap of results
ggplot(d_mean, aes(x=reo_meth(alg), y=sp_y)) + 
  geom_tile(aes(fill=pchange_cut), colour='white') +
  scale_y_discrete(limits = unique(rev(d_mean$sp_y))) + 
  geom_text(aes(label = round(r, 2))) +
  labs(x="method", y="species") + theme_minimal(10) +
  scale_fill_gradient2(low='firebrick3', mid = "darkgoldenrod1", high = "gold",midpoint = 0.10, limits = c(0,0.20))


# Calcuate median % of optimal (outputs % decrease from best, just subtract from 100%)
tmp <- subset(d_mean, y == 'HT')[,c('alg','pchange')]
tmp$pchange <- as.numeric(tmp$pchange)
ht_stats <- aggregate(pchange ~ alg, tmp, FUN='median')
tmp <- subset(d_mean, y %in% c('FT', 'AN', 'R8'))[,c('alg','pchange')]
tmp$pchange <- as.numeric(tmp$pchange)
ft_stats <- aggregate(pchange ~ alg, tmp, FUN='median')
tmp <- subset(d_mean, y == 'YLD')[,c('alg','pchange')]
tmp$pchange <- as.numeric(tmp$pchange)
yld_stats <- aggregate(pchange ~ alg, tmp, FUN='median')
tmp <- subset(d_mean, y %in% c('MO', 'DBH', 'DE', 'ST'))[,c('alg','pchange')]
tmp$pchange <- as.numeric(tmp$pchange)
othr_stats <- aggregate(pchange ~ alg, tmp, FUN='median')


## Cluster algorithms by performance
library(clValid)
d$ID <- paste(d$tag, d$y, sep=":")
keep = c('rrBLUP', 'BRR', 'BayesA', 'BayesB', 'BL', 'SVM','SVMpoly','SVMrbf', 'RF', 'GB', 'ANN_FS','ANN_BayesB_8k', 'EN5', 'EN11')
d4clust <- subset(d, d$model %in% keep)
d4clustm <- dcast(d4clust[,c('model', 'ID', 'PCC')], model ~ ID,  fun.aggregate=mean)
row.names(d4clustm) <- d4clustm$model
d4clustm$model <- NULL

distance <- dist(d4clustm, method='euclidean')
hc1 <- hclust(distance, method='ward.D')
plot(hc1, cex=0.6, hang=-1)




###################################
## Seeded ANN Results for Height ##
###################################
setwd('/Volumes/azodichr/03_GenomicSelection/06_otherTraits/mlp_ens/')

### Plot for height
bestx <- rbind(bglr, rrb, ml, annFS)
bestx$tag <- gsub('sorgh', 'sorg', bestx$tag)
lrann_ht <- rbind(lrann, lrann2, annFS)
lrann_ht$tag <- gsub('sorgh', 'sorg', lrann_ht$tag)

best_mean <- aggregate(PCC~tag+model+y, bestx, FUN=mean)
best <- aggregate(PCC~tag+y, best_mean, function(x) x[which.max(abs(x))])
best$set <- 'best'

lrann_best <- merge(lrann_ht, best, by=c('tag', 'y'))
names(lrann_best) <- c('species', 'y','alg', 'r', 'r_best', 'x2')
lrann_best_ht <- subset(lrann_best, y=='HT')
reo_ann <- function(x) { factor(x, levels = c('ANN_FS', 'ANN_rrB_8k','ANN_BayesB_8k', 'ANN_BL_8k', 'ANN_RF_8k'))}

lrann_best_8k <- subset(lrann_best_ht, alg %in% c('ANN_FS', 'ANN_rrB_8k','ANN_BayesB_8k', 'ANN_BL_8k', 'ANN_RF_8k'))
ggplot(lrann_best_8k, aes(x=reo_ann(alg), y=r)) +
  geom_violin() + theme_bw() + facet_grid(.~species) +
  geom_boxplot(outlier.colour = NA, width=.2) +
  geom_hline(data=lrann_best_ht, aes(yintercept = r_best)) +
  coord_cartesian(ylim=c(0,0.7))
ggsave("../00_Figures/190301_HT_ANN_LR.pdf", width = 4, height = 2)

lrann_best_8k_spMean <- aggregate(lrann_best_8k$r, by=list(lrann_best_8k$species, lrann_best_8k$alg), FUN='mean')
colnames(lrann_best_8k_spMean) <- c('species','alg','r')
summary(with(lrann_best_8k_spMean, aov(r ~ alg + Error(species/alg))))

# Stats
summary(aov(r~alg*species, lrann_best_8k))
summary(aov(r~alg, subset(lrann_best_8k, species=='maize')))
summary(aov(r~alg, subset(lrann_best_8k, species=='rice')))
summary(aov(r~alg, subset(lrann_best_8k, species=='sorg')))
summary(aov(r~alg, subset(lrann_best_8k, species=='soy')))
summary(aov(r~alg, subset(lrann_best_8k, species=='spruce')))
summary(aov(r~alg, subset(lrann_best_8k, species=='swgrs')))

wilcox.test(r~alg, subset(lrann_best_8k, species=='rice' & alg %in% c('ANN_FS', 'ANN_rrB_8k')& r != 0.5106800))
wilcox.test(r~alg, subset(lrann_best_8k, species=='rice' & alg %in% c('ANN_FS', 'ANN_BayesB_8k')& r != 0.5106800))
wilcox.test(r~alg, subset(lrann_best_8k, species=='rice' & alg %in% c('ANN_FS', 'ANN_BL_8k') & r != 0.5106800))
wilcox.test(r~alg, subset(lrann_best_8k, species=='rice' & alg %in% c('ANN_FS', 'ANN_RF_8k')& r != 0.5106800))


# Stdev between replicates for each algorithm
stdev <- aggregate(r ~ species+alg, lrann_best_8k, FUN=sd)
ggplot(stdev, aes(x=reo_ann(alg), y=r)) + geom_violin() + theme_bw()+
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = species))

summary(aov(r~alg, stdev))
wilcox.test(r~alg, subset(stdev, alg %in% c('ANN_FS', 'ANN_rrB_8k'), paired=FALSE))
wilcox.test(r~alg, subset(stdev, alg %in% c('ANN_FS', 'ANN_BayesB_8k'), paired=TRUE))
wilcox.test(r~alg, subset(stdev, alg %in% c('ANN_FS', 'ANN_BL_8k'), paired=TRUE))
wilcox.test(r~alg, subset(stdev, alg %in% c('ANN_FS', 'ANN_RF_8k'))) #, paired=TRUE))



### Compare seeding 25% and 50% of nodes
r25 <- read.csv('mlp_ens_2/fs_8000/02_UseInManuscript/RESULTS.txt',sep='\t', header=T, row.names = NULL)
r25$seed <- '25p'
r50 <- read.csv('mlp_ens_2/fs_8000/RESULTS.txt',sep='\t', header=T, row.names = NULL)
r50$seed <- '50p'
r <- rbind(r25, r50)
r$pc_TrVal <- (r$MSE_Train - r$MSE_Valid) / r$MSE_Train

summary(aov(MSE_test ~ seed, r))
summary(aov(MSE_Valid ~ seed, r))
aggregate(MSE_test ~ seed, r, FUN=mean)
aggregate(MSE_Valid ~ seed, r, FUN=mean)
aggregate(PCC_test ~ seed, r, FUN=mean)

##########################################################################
###########  % winner comparison for all species/traits ##################
##########################################################################
p <- read.csv('PercentOutperformResults.csv', header=T)

keep <- c('rrBLUP', 'BRR', 'BayesA', 'BayesB', 'BL', 'SVM','SVMpoly','SVMrbf', 'RF', 'GB', 'ANN','ANN_BayesB','ANN_rrB','ANN_BL','ANN_RF')
reo_meth <- function(x) { factor(x, levels = keep)}
reo_meth2 <- function(x) { factor(x, levels = rev(keep))}

ggplot(p, aes(x=reo_meth(m2), y=reo_meth2(m1))) + 
  geom_tile(aes(fill=m1_per_wins), colour='white') +
  facet_wrap(~ID) +
  labs(y="winner", x="looser") + 
  scale_fill_gradient2(low='turquoise4', mid = "white", high = "firebrick3",midpoint = 0.5, limits = c(0,1))+
  theme_minimal(10) 
ggsave("../00_Figures/190302_percent_winner.pdf", width = 10, height = 10)

