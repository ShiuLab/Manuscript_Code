# Benchmarking algorithms for genomic prediction of complex traits

**Manuscript Goal:** Using data from 18 traits across six plant species with different marker densities and training population sizes, compare the performance of six linear and five non-linear algorithms, including artificial neural networks. 


## Data Format Required

**Phenotype Data:**

Matrix with rows representing lines and columns representing the phenotypes. The first column should be the unique line identifiers and should match the genotype data set. Phenotypes should be normalized to improve performance of certain algorithms such as SVR and ANNs.

**Genotype Data:**

Matrix with rows representing lines and columns representing allele calles for each genetic marker. The first column should be the unique line identifiers and should match the phenotype data set. Allele calls should be coded as follows: homozygous major/reference allele = 1, heterozygous = 0, homozygous minor/non-reference allele = -1.


## Define test set

Models are generated using a training, validation, testing scheme. To facilitate comparing new algorithms using this benchmark dataset, the testing sets for each species are available [here]() (Note, link will go live at the time of publication)

Example code to randomly generate a test set using the [ML_Pipeline](https://github.com/ShiuLab/ML-Pipeline/):
```
python ML_Pipeline/test_set.py -df soy_geno.csv -df2 soy_pheno.csv -p 0.2 -type r -y_name HT -sep ',' -save soy_testset_1.txt
```


## rrBLUP

The rrBLUP algorithm was implemented in R using [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf) [cite](https://dl.sciencesocieties.org/publications/tpg/abstracts/4/3/250). Packages needed: rrBLUP, data.table

Rscript predict_rrBLUP.R [Geno_file] [Pheno_file] [Features2Use] [Y_name] [test_set] [save_ID] [save_dir]

Example code to train and test rrBLUP models using [predict_rrBLUP.R](https://github.com/ShiuLab/GenomicSelection):
```
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv all YLD soy_testset_1.txt soy rrblup_results/
```

## BGLR

The Bayesian LASSO (BL), Bayesian ridge regression (BRR), Bayes A (BayesA), and Bayes B (BayesB) algorithms were implemented in R using [BGLR](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf) [cite](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4196607/). Packages needed: BGLR, data.table

Rscript predict_rrBLUP.R [Geno_file] [Pheno_file] [Features2Use] [Y_name] [BGLR_model] [test_set] [save_ID] [save_dir]

Example code to train and test BGLR models using [predict_BGLR.R](https://github.com/ShiuLab/GenomicSelection):
```
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv all YLD BL soy_testset_1.txt soy bl_results/
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv all YLD BRR soy_testset_1.txt soy brr_results/
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv all YLD BayesA soy_testset_1.txt soy ba_results/
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv all YLD BayesB soy_testset_1.txt soy bb_results/
```


## Machine Learning Pipeline

The support vector regression (SVR), Random Forest (RF), and Gradient Tree Boosting (GTB) algorithms were implemented in python using [Scikit-learn](https://scikit-learn.org/stable/) [cite](http://www.jmlr.org/papers/volume12/pedregosa11a/pedregosa11a.pdf) in the Shiu Lab's ML Pipeline. For a detailed overview of how to use this pipeline for genomic prediction or any kind of machine learning problem, please check out the [tutorial](https://github.com/ShiuLab/ML-Pipeline/tree/master/Tutorial).

For a summary of all options for the ML pipeline, run the command without any arguments. 

Example code
```
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVM -gs T -cv 5 -n 10 -tag soy_HT_SVRlin -save soy_HT_SVRlin 
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVMpoly -gs T -cv 5 -n 10 -tag soy_HT_SVRpoly -save soy_HT_SVRpoly
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVMrbf -gs T -cv 5 -n 10 -tag soy_HT_SVRrbf -save soy_HT_SVRrbf
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg RF -gs T -cv 5 -n 10 -tag soy_HT_RF -save soy_HT_RF
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg GTB -gs T -cv 5 -n 10 -tag soy_HT_GTB -save soy_HT_GTB
```




