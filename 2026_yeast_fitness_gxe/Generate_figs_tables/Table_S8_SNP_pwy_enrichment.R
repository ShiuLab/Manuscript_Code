################################################################################
# Pathway enrichment (overrepresentation) analysis of top random forest 
# features in each condition       
################################################################################

rm(list=ls())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GSEABase))

########### Reshape pathway data (doesn't have pathway descriptions) ###########
setwd("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/")
pwy <- read.csv("Co-function/Data/MetaCyc/All-genes-pathways-S288c.txt", sep="\t")
pwy2 <- pwy %>% separate_rows(Pathways.of.gene, sep="[/]+") %>% as.data.frame()
write.csv(pwy2, "Co-function/Data/MetaCyc/All-genes-pathways-S288c_pivoted.txt", quote=F, row.names=F)

################################################################################
#                         Pathway Enrichment Analysis                          #
#                    contingency table for each GO term:                       #
#                                  | Gene in top | Gene not top                #
#                   -----------------------------------------                  #
#                    In pathway    |             |                             #
#                    Not in pathway|             |                             #
################################################################################
# Prep gene-pathway map for background set
setwd("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")
all_genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv", sep="\t", header=F)
colnames(all_genes) <- c("snp", "chr", "pos", "gene")
all_genes <- all_genes[all_genes$gene!="intergenic",] # drop intergenic snps
all_genes <- all_genes[!grepl(",", all_genes$gene),] # drop snps that mapped to multiple genes

pwys <-  read.csv("../Co-function/Data/MetaCyc/All-genes-pathways-S288c_pivoted.txt")
all_genes <- left_join(all_genes, pwys, by=c("gene"="Accession.1"), relationship="many-to-many")
all_genes <- all_genes[,c("gene", "Pathways.of.gene")]
all_genes <- all_genes[!duplicated(all_genes),]
colnames(all_genes) <- c("gene", "pathway")
all_genes$pathway <- trimws(all_genes$pathway, which="both") # strip white spaces
all_genes$pathway <- na_if(all_genes$pathway, "") # replace "" with NA
write.csv(all_genes, "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED_pathway_map.csv", quote=F, row.names=F)

enrichment <- function(k, n, C, G){ 
    # determine direction of enrichment
    # if >= 1: + (overrepresented)
    # if < 1: - (underrepresented)
    # k: number of genes in target_list with GO
    # n: total number of genes in target_list
    # C: total number of genes (in target_list + background) with GO
    # G: total number of genes (in target_list + background)
    return((k/C)/(n/G))
}

ora <- function(with_pwy, top, bg, path){
    # Overrepresentation Analysis
    # with_pwy: dataframe of all genes and GO annotations
    # top: dataframe of genes of interest
    # bg: dataframe of genes in background set
    # path: file path and name to save as

    # create contingency table
    cols <- c("PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY", "Gene_top_no_PWY",
              "Gene_not_top_no_PWY", "direction", "p.val", "odds ratio", "qvalues")
    contingency <- data.frame(matrix(nrow=1, ncol=9))
    colnames(contingency) <- cols

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (pwy in unique(all_genes$pathway)){
        a <- length(unique(top[which(top$pathway==pwy),]$gene)) # Genes in top features and have `pwy`
        b <- length(unique(bg[which(bg$pathway==pwy),]$gene)) # Genes not in top features and have `pwy`
        c <- length(unique(top$gene)) - a # Genes in top features and do not have `pwy`
        d <- length(unique(bg$gene)) - b # Genes not in top features and do not have `pwy`
        tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        if (a+b!=0){
            if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # direction of enrichment
            contingency <- rbind(contingency, list(pwy, a, b, c, d, direction, res$p.value, res$estimate, 'NA'))
        }
    }
    contingency <- contingency[!is.na(contingency$PWY),] # drop first row with NAs
    contingency <- contingency[contingency$PWY!="",] # drop rows with empty PWY
    
    # Calculate q-values
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # add biological process, cellular component, and molecular function info
    if (nrow(sub)!=0){
        print("   Calculating q values...")
        qvals <- p.adjust(sub$p.val, method="BH")
        sub$qvalues <- qvals

        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, path, sep="\t", quote=F, row.names=F)
    }
}

pwy_enrichment <- function(f){
    ### ORA of top gene features
    # read in top gene feature file
    top <- read.csv(f)
    
    # drop intergenic snps
    top <- top[top$gene!="intergenic",]
    top <- top[complete.cases(top$gene),]
    
    # add pathway information
    top <- left_join(top, all_genes, by=c("gene"="gene"), relationship="many-to-many")

    # make background set
    bg <- all_genes[!(all_genes$gene %in% top$gene),] # remove top from background set
    print(paste("   Top Genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    print(paste("   Total number of genes is correct: ", length(unique(top$gene))+length(unique(bg$gene))==length(unique(all_genes$gene))))
    
    ## Overrepresentation Analysis
    path <- gsub("Genes_", "ORA_PWY_Genes_", f)
    path <- gsub("csv", ".tsv", path)
    path <- gsub("GO_Enrichment", "PWY_Enrichment", path)
    ora(all_genes, top, bg, path)
}


# Read in top features' (FS) average SHAP values files
dir <- "/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/GO_Enrichment/SNPs_fs" # path to optimized gini files
files <- list.files(path=dir, pattern="^Genes_", full.names=TRUE, recursive=FALSE)

mclapply(X=files, FUN=pwy_enrichment, mc.cores=35) #35 match go to orfs

