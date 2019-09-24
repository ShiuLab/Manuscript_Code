# Transcriptome-based prediction of complex traits in maize

## Data Processing/Exploration

Code to process genotype, transcriptome, and phenotype data as described in (Azodi et al. TBD)[link coming!]: scripts data_preprocessing.R

Code to generate figures demonstrating the relationship between G, T, and P: scripts/plot_correlations.R


## Modeling 

For each algorithm, run model using kinship (K), genotype (G), expression correlation (eCOR), transcript (T),
T+K, and T+G as input features. Showing command line examples just using transcript (T) data

**Bayesian LASSO** 

Packages needed:

- data.table 
- BGLR

```
for i in $(seq 1 100); do Rscript scripts/predict_BGLR.R data/transcripto.csv data/pheno.csv FT data/CVFs.csv $i BL results/ ; done
```

**rrBLUP**

Packages needed:

- data.table
- rrBLUP

```
for i in $(seq 1 100); do Rscript scripts/predict_rrBLUP.R data/transcripto.csv data/pheno.csv FT data/CVFs.csv $i /results/ ; done
```

**Random Forest**

For the most recent version of the ML pipeline (and detailed instructions on how to use the ML pipeline for regression/classification see the [ML_Pipeline Repository](https://github.com/ShiuLab/ML-Pipeline).

```
for i in $(seq 1 100); do python scripts/ML_regression.py -df data/transcripto.csv -df2 data/pheno.csv -y_name FT -sep ',' -alg RF -gs T -cv_set data/CVFs.csv -save trans_RF_FT; done
```

**Ensemble**

Need output from rrBLUP, BL, and RF to run

```
python scripts/ensemble_pred.py ensem_transcripto_FT_yhat.csv data/pheno.csv FT transcripto rrBLUP_transcripto_FT_yhat.csv BL_transcripto_FT_yhat.csv transcriptoFT_RF_scores.txt
python scripts/ensemble_imp.py ensem_transcripto_FT_imp rrBLUP_transcripto_FT_coef.csv BL_transcripto_FT_coef.csv transcriptoFT_RF_imp
```

**Plots**

*Code:* scripts/plots.R



## eQTL Analysis

Find significant cis and trans eQTL using [MatrixeQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/).

*Code:* scripts/eQTL_analysis.R



## Comparison of importance scores

The importance scores are output from rrBLUP, BL, RF, En scripts above.

**Between T:G pairs**

```
python scripts/find_T_G_pairs.py -g BL_geno_FT_coef.csv -t BL_transcripto_FT_coef.csv -key data/v3_v4_xref.txt -save BL_FT_SNP2kb -LD 1000
Rscript ../../scripts/plot_importance_topSNPs_ggp.R BL_FT_SNP2kb.csv BL_FT_coef_comp
```

**Between T:eQTL pairs**

```
python scripts/find_T_eQTL_pairs.py -cis eqtl_cis.csv -trans eqtl_trans.csv -g ensem_geno_FT_imp -t ensem_transcripto_FT_imp -key data/v3_v4_xref.txt -type imp -save ensem_eqtl_all_pairs```
Rscript scripts/plot_impCOR_eQTL.R ensem_FT_SNP2kb.csv ensem_eqtl_all_pairs_eQTL.csv ensem_eqtl_all_pairs_cor
```

**Between algorithms**

```
Rscript scripts/plot_impCOR_betweenModels.R FT rrBLUP_FT_SNP2kb.csv BL_FT_SNP2kb.csv RF_FT_SNP2kb.csv
```



## Flowering Time Benchmark Analysis

**Pull percentile, rank, and z-score for each algorithm type for SNPs and transcripts**

```
Rscript scripts/importance_statistics.R rrBLUP_transcripto_FT_coef.csv coef
Rscript scripts/importance_statistics.R rrBLUP_geno_FT_coef.csv coef
Rscript scripts/importance_statistics.R ensem_transcripto_FT_imp imp
Rscript scripts/importance_statistics.R ensem_geno_FT_imp imp
```

**Run t-tests and linear models on allele type and transcript levels VS. flowering time**

```Rscript scripts/summary_stats_GvFT_TvFT.R```


**Generate figures**

*Code:* scripts/plot_benchmarkAnalysis.R




## Side Analysis

### Feature selection for FT

*Define hold out and run Feature Selection using Random Forest*

```
for i in $(seq 1 10); do python ML-Pipeline/holdout.py -df data/pheno.csv -type r -p 0.1 -sep ',' -save holdout_"$i".txt; done
```

*Run feature selection using Random Forest*

```
for i in $(seq 1 10); do python ML-Pipeline/Feature_Selection.py -f RF -type r -n 50,100,250,500,750,1000,2500,5000 -df data/transcripto.csv -df2 data/pheno.csv -y_name FT -sep ',' -ho holdout_"$i".txt -save t_ho"$i"; done
```

*Then subset dataframe to include only markers/transcripts selected by RF and run rrBLUP, BL, and RF as above*




