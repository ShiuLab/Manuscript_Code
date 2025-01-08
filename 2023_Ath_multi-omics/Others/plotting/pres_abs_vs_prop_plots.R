# To use this script:
# Parallelization doesn't work via HPCOnDemand, therefore this script has to be
# run from the command line. To do this, you first need to run 
# module load GCC/10.2.0  OpenMPI/4.0.5
# followed by 
# module load R/4.0.3
# If you just do module load R, it will point to a version of R that doesn't
# have ggpubr installed & will throw an error.
# After doing this, run the script via 
# Rscript pres_abs_vs_prop_plots.R

library(data.table)
library(parallel)
library(ggplot2)
library(ggpubr)


path = '/mnt/home/lotrecks/Shiu_lab/multi-omics/outputs/rrblup'
outpath = '/mnt/home/lotrecks/Shiu_lab/multi-omics/experimental_plotting'
setwd(path)

files = c(
      'Coef_multi_omics_binned_methylation_data_CG_pres_abs_mean_feature_matrix_FT10_mean.csv',
      'Coef_multi_omics_binned_methylation_data_CG_prop_median_feature_matrix_FT10_mean.csv',
      'Coef_multi_omics_binned_methylation_data_CHG_pres_abs_mean_feature_matrix_FT10_mean.csv',
      'Coef_multi_omics_binned_methylation_data_CHG_prop_median_feature_matrix_FT10_mean.csv',
      'Coef_multi_omics_binned_methylation_data_CHH_pres_abs_mean_feature_matrix_FT10_mean.csv',
      'Coef_multi_omics_binned_methylation_data_CHH_prop_median_feature_matrix_FT10_mean.csv')

dtypes = c('pres_abs_mean', 'prop_median')
methtypes = c('CG', 'CHG', 'CHH')

plot_corrs <- function(files, dtypes, methtypes){
   
   for (methtype in methtypes){
      
      # Get files for this set of graphs
      cat(sprintf('\n\nGetting files for graph for %s...', methtype))
      methtype_files = files[grepl(methtype, files, fixed = TRUE)]
      
      # Get the dtype and pheno name for each file
      cat(sprintf('\nGetting dataset names...'))
      df_names = get_data_names(methtype_files, dtypes)
      
      # Read in the files
      # NOTE: This does not work on OnDemand, must be run from the command line
      cat(sprintf('\nReading in the files...'))
      clust <- makeCluster(detectCores())
      df_list <- parLapply(clust, methtype_files, fread)
      stopCluster(clust)
      
      # Assign names to files 
      cat(sprintf('\nNaming files...'))
      names(df_list) = df_names
      
      # Average the rows and make into one df
      cat(sprintf('\nAveraging rows...'))
      df_means_list = lapply(df_list, colMeans)
      plot_data = data.frame(df_means_list[[1]], df_means_list[[2]])
      
      # Make plots
      cat(sprintf('\nMaking plots...'))
      # plot_geom_point(plot_data, "df_means_list..1..", "df_means_list..2..", df_names, "pearson", outpath, methtype)
      # plot_geom_point(plot_data, "df_means_list..1..", "df_means_list..2..", df_names, "spearman", outpath, methtype)
      # plot_geom_bins(plot_data, "df_means_list..1..", "df_means_list..2..", df_names, "pearson", outpath, methtype)
      # plot_geom_bins(plot_data, "df_means_list..1..", "df_means_list..2..", df_names, "spearman", outpath, methtype)
      
      # Hacky but works 
      save_path = sprintf('%s/%s_pearson_coefficient_correlation_point.pdf', outpath, methtype)
      pdf(save_path)
      p = ggplot(plot_data, aes(df_means_list..1.., df_means_list..2..)) +
         geom_point(size=1,alpha = 0.01) +
         ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
         labs(y=df_names[2], x=df_names[1]) +
         stat_cor(method="pearson")
      plot(p)
      dev.off()
      
      save_path = sprintf('%s/%s_spearman_coefficient_correlation_point.pdf', outpath, methtype)
      pdf(save_path)
      p = ggplot(plot_data, aes(df_means_list..1.., df_means_list..2..)) +
         geom_point(size=1, alpha=0.01) +
         ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
         labs(y=df_names[2], x=df_names[1]) +
         stat_cor(method="spearman")
      plot(p)
      dev.off()
      
      save_path = sprintf('%s/%s_pearson_coefficient_correlation_geom.pdf', outpath, methtype)
      pdf(save_path)
      p = ggplot(plot_data, aes(df_means_list..1.., df_means_list..2..)) +
         geom_bin2d(bins=100) +
         ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
         labs(y=df_names[2], x=df_names[1]) +
         stat_cor(method="pearson")
      plot(p)
      dev.off()
      
      save_path = sprintf('%s/%s_spearman_coefficient_correlation_geom.pdf', outpath, methtype)
      pdf(save_path)
      p = ggplot(plot_data, aes(df_means_list..1.., df_means_list..2..)) +
         geom_bin2d(bins=100) +
         ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
         labs(y=df_names[2], x=df_names[1]) +
         stat_cor(method="spearman")
      plot(p)
      dev.off()
   }
}

get_data_names <- function(methtype_files, dtypes){

   dtype_matches = sapply(methtype_files, mygreplsapply, patterns=dtypes)
   col_names = colnames(dtype_matches)
   names = unname(sapply(col_names, get_true_row_name, matches=dtype_matches))
   
   return(names)   
}

mygreplsapply <- function(patterns, str_to_match){
   return(sapply(patterns, grepl, x=str_to_match))
}

get_true_row_name <- function(matches, col){
   true_name = row.names(matches)[which(matches[,col]==TRUE)]
   
   return(true_name)
}

# For some reason, compartmentalizing this code into function doesn't work - 
# I believe it's because when I do it this way, it doesn't let me pass the column
# names to aes as objects, and makes me quote them - I think something bad is 
# happening there, but just going to do by copy-paste for now
# plot_geom_point <- function(data, aes1, aes2, df_names, corr, outpath, methtype){
#    save_path = sprintf('%s/%s_%s_coefficient_correlation_point.pdf', outpath, methtype, corr)
#    pdf(save_path)
#    p = ggplot(data, aes(aes1, aes2)) +
#       geom_point(size=1,alpha = 0.01) +
#       ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
#       labs(y=df_names[2], x=df_names[1]) +
#       stat_cor(method=corr)
#    plot(p)
#    dev.off()
# }
# 
# plot_geom_bins <- function(data, aes1, aes2, df_names, corr, outpath, methtype){
#    save_path = sprintf('%s/%s_%s_coefficient_correlation_geom.pdf', outpath, methtype, corr)
#    pdf(save_path)
#    p = ggplot(data, aes(aes1, aes2)) +
#       geom_bin2d(bins=1000) +
#       ggtitle(sprintf('Correlation of rrBLUP coefficients for %s methylation', methtype)) +
#       labs(y=df_names[2], x=df_names[1]) +
#       stat_cor(method=corr)
#    plot(p)
#    dev.off()
# }

# Call function
plot_corrs(files, dtypes, methtypes)