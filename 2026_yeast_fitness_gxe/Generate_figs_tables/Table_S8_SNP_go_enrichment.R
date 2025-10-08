################################################################################
# Go enrichment analysis of top random forest features in each condition       #
################################################################################

rm(list=ls())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
# suppressPackageStartupMessages(library(GSEABase))

setwd("/mnt/research/glbrc_group/shiulab/kenia/Shiu_Lab/Project")

#### 2022-09-20T19:04 S288C GO term annotation file from GO database
go <- fread("Data/yeast_GO/sgd.gaf", sep="\t", skip=35, header=F, fill=T)
go <- as.data.frame(go)
colnames(go) <- c("DB", "DB.ID", "Protein", "Qualifier", "GO.ID", "DF.Ref",
    "Evidence", "With.or.From", "Aspect", "Gene.Description", "Gene.Synonym",
    "Type", "Taxon", "Date", "Assigned.By", "Annotation.Extension",
    "Gene.Product") # GAF Format Column Names

# Biological process go term gene counts
tmp <- aggregate(go$Gene.Synonym, by=list(go$GO.ID), FUN=length)
tmp[tmp$Group.1=='GO:0008150',]
#         Group.1    x
# 1899 GO:0008150 1075

# Filter GO terms by evidence code
exp <- c("IDA","IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP") # EXP & HTP evidence codes
go <- subset(go, Evidence %in% exp)

# Extract gene name from Gene.Synonym
unique(go$Type)
Gene <- str_extract(go$Gene.Synonym, "[A-Z0-9-]+|")
go <- cbind(go, Gene)
# write.table(go, "Data/yeast_GO/sgd_GO_BP_no_info.txt", sep="\t", quote=F, row.names=F)
go <- fread("Data/yeast_GO/sgd_GO_BP_no_info.txt", sep="\t")

############################# Map GO terms to SNPs #############################
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_CORRECTED.tsv",
                  sep="\t", header=F)
colnames(genes) <- c("snp", "chr", "pos", "gene")
genes <- genes[genes$gene!="intergenic",] # drop intergenic snps
genes <- genes[!grepl(",", genes$gene),] # drop snps that mapped to multiple genes
out <- left_join(genes, go, by=c("gene"="Gene"))
# write.table(out, 
#     "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go_all.tsv", sep="\t")

################################################################################
#                           GO Enrichment Analysis                             #
#                    contingency table for each GO term:                       #
#                                  | ORF in top | ORF not top                  #
#                   -----------------------------------------                  #
#                    In GO term    |            |                              #
#                    Not in GO term|            |                              #
################################################################################
all_genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go_all.tsv", sep="\t")

get_genes <- function(f, baseline=F) {
    # add the genes and GO annotations to each SHAP file
    df <- read.delim(f, sep="\t") # read in shap values
    out <- right_join(all_genes, df, by=c("snp"="actual_feature")) # add gene & go information
    name <- str_extract(f, "[A-Z0-9]+_[a-z_]+_[0-9]+_imp") # extract file name
    name <- paste("Genes_", name, sep="")
    path <- as.character("Scripts/Data_Vis/Section_3/GO_Enrichment/SNPs_fs")
    out <- out[,c(1,4,9,33)]
    write.csv(out, paste(file.path(path, name), ".csv", sep=""), quote=F, row.names=F)
    return(out)
}

enrichment <- function(k, n, C, G){ 
    # determine direction of enrichment
    # if >= 1: + (overrepresented)
    # if < 1: - (underrepresented)
    # k: number of genes in target list with GO
    # n: total number of genes in target list
    # C: total number of genes (in target list + background) with GO
    # G: total number of genes (in target list + background)
    return((k/C)/(n/G))
}

ora <- function(all_genes, top, bg, path){
    # Overrepresentation Analysis
    # all_genes: dataframe of all genes and GO annotations
    # top: dataframe of genes of interest
    # bg: dataframe of genes in background set
    # path: file path and name to save as

    # create contingency table
    cols <- c("GO", "Gene_top_has_GO", "Gene_not_top_has_GO", "Gene_top_no_GO",
              "Gene_not_top_no_GO", "direction", "p.val", "odds ratio", "qvalues")
    contingency <- data.frame(matrix(nrow=1, ncol=9))
    colnames(contingency) <- cols
    
    # fill in contingency table for each gene
    print("   Running ORA...")
    for (go in unique(all_genes$GO)){
        if (!is.na(go)){
            a <- length(unique(top[which(top$GO==go),]$gene)) # Genes in top features and have `go`
            b <- length(unique(bg[which(bg$GO==go),]$gene)) # Genes not in top features and have `go`
            c <- length(unique(top$gene)) - a # Genes in top features and do not have `go`
            d <- length(unique(bg$gene)) - b # Genes not in top features and do not have `go`
            tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
            res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
            if (a+b!=0){
                if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # direction of enrichment
                contingency <- rbind(contingency, list(go, a, b, c, d, direction, res$p.value, res$estimate, 'NA'))
            }
        }
    }
    contingency <- contingency[!is.na(contingency$GO),] # drop first row with NAs

    # Calculate q-values
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # add biological process, cellular component, and molecular function info
    if (nrow(sub)!=0){
        # print("   Grabbing BP, CC, and MF info...")
        # sub$BP <- "" # biological process
        # sub$CC <- "" # cellular component
        # sub$MF <- "" # molecular function
        # for(i in 1:nrow(sub)){
        #     tryCatch({
        #         if(!is.null(getGOTerm(sub[i,1])$BP[1])) sub[i,10] <- getGOTerm(sub[i,1])$BP[1]
        #         if(!is.null(getGOTerm(sub[i,1])$CC[1])) sub[i,11] <- getGOTerm(sub[i,1])$CC[1]
        #         if(!is.null(getGOTerm(sub[i,1])$MF[1])) sub[i,12] <- getGOTerm(sub[i,1])$MF[1]
        #     }, error = function(e){print(paste("no GO for ", sub[i,1])); NaN},
        #         finally = {})
        # }
        print("   Calculating q values...")
        qvals <- p.adjust(sub$p.val, method="BH")
        sub$qvalues <- qvals

        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, paste(path, ".tsv", sep=""), sep="\t", quote=F, row.names=F)
    }
}

go_enrichment <- function(f){
    ### ORA and GSEA of top gene features
    # add gene information to top snp file
    top <- get_genes(f)
    colnames(top) <- c("snp", "gene", "GO", "mean_imp")
    
    # drop intergenic snps from top
    top <- top[which(top$gene!="intergenic"),]
    
    # prepare background set
    bg <- all_genes[!(all_genes$gene %in% top$gene),] # remove top from background set
    bg <- bg[c(4,9)] # keep only gene, GO.ID
    colnames(bg) <- c("gene", "GO")
    bg <- bg[!duplicated(bg),] # remove duplicates
    print(paste("   Top Genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    print(paste("   Total number of genes is correct: ", length(unique(top$gene))+length(unique(bg$gene))==length(unique(all_genes$gene))))
    
    ## Overrepresentation Analysis
    path <- paste("Scripts/Data_Vis/Section_3/GO_Enrichment/SNPs_fs/ORA_Genes_",
                  str_extract(f, "[A-Z0-9]+_[a-z_]+_[0-9]+_imp"), sep="")
    ora(all_genes, top, bg, path)
}

# Read in top features' (after FS) RF importance score files
dir <- "/mnt/research/glbrc_group/shiulab/kenia/yeast_project/SNP_yeast_RF_results/fs" # path to optimized gini files
files <- list.files(path=dir, pattern="_imp_with_actual_feature_names_09102025$", full.names=TRUE)

mclapply(X=files, FUN=go_enrichment, mc.cores=35) #35 match go to orfs
