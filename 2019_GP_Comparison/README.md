# Benchmarking algorithms for genomic prediction of complex traits

**Manuscript Goal:** Using data from 18 traits across six plant species with different marker densities and training population sizes, compare the performance of six linear and five non-linear algorithms, including artificial neural networks. 

**Code Versions:** Code in this repository is the version used for the GP Comparision manuscript. For updated versions of code, see other Shiu lab repositories designated for each type of algorithm (links below).


## Data Format Required

**Phenotype Data:**

Matrix with rows representing lines and columns representing the phenotypes. The first column should be the unique line identifiers and should match the genotype data set. Phenotypes should be normalized to improve performance of certain algorithms such as SVR and ANNs.

**Genotype Data:**

Matrix with rows representing lines and columns representing allele calles for each genetic marker. The first column should be the unique line identifiers and should match the phenotype data set. Allele calls should be coded as follows: homozygous major/reference allele = 1, heterozygous = 0, homozygous minor/non-reference allele = -1.

**Data used in original benchmark available at: LINK COMING SOON**

## Define test set

Models are generated using a training, validation, testing scheme. To facilitate comparing new algorithms using this benchmark dataset, the testing sets for each species are available [here]() (Note, link will go live at the time of publication)

Example code to randomly generate a test set using the [ML_Pipeline](https://github.com/ShiuLab/ML-Pipeline/):
```
python ML_Pipeline/test_set.py -df soy_geno.csv -df2 soy_pheno.csv -p 0.2 -type r -y_name HT -sep ',' -save soy_testset_1.txt
```


## rrBLUP

The rrBLUP algorithm was implemented in R using [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/rrBLUP.pdf) [cite](https://dl.sciencesocieties.org/publications/tpg/abstracts/4/3/250). Packages needed: rrBLUP, data.table.

Rscript predict_rrBLUP.R [Geno_file] [Pheno_file] [Features2Use] [Y_name] [test_set] [save_ID] [save_dir]

For more details about how to train and test rrBLUP models using [predict_rrBLUP.R](https://github.com/ShiuLab/GenomicSelection):
```
Rscript predict_rrBLUP.R soy_geno.csv soy_pheno.csv all HT soy_testset_1.txt soy rrblup_results/
```

## BGLR

The Bayesian LASSO (BL), Bayesian ridge regression (BRR), Bayes A (BayesA), and Bayes B (BayesB) algorithms were implemented in R using [BGLR](https://cran.r-project.org/web/packages/BGLR/BGLR.pdf) [cite](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4196607/). Packages needed: BGLR, data.table

Rscript predict_rrBLUP.R [Geno_file] [Pheno_file] [Features2Use] [Y_name] [BGLR_model] [test_set] [save_ID] [save_dir]

For more details about how to train and test BGLR models using [predict_BGLR.R](https://github.com/ShiuLab/GenomicSelection):
```
Rscript predict_BGLR.R soy_geno.csv soy_pheno.csv all HT BL soy_testset_1.txt soy bl_results/
Rscript predict_BGLR.R soy_geno.csv soy_pheno.csv all HT BRR soy_testset_1.txt soy brr_results/
Rscript predict_BGLR.R soy_geno.csv soy_pheno.csv all HT BayesA soy_testset_1.txt soy ba_results/
Rscript predict_BGLR.R soy_geno.csv soy_pheno.csv all HT BayesB soy_testset_1.txt soy bb_results/
```


## Machine Learning Pipeline

The support vector regression (SVR), Random Forest (RF), and Gradient Tree Boosting (GTB) algorithms were implemented in python using [Scikit-learn](https://scikit-learn.org/stable/) [cite](http://www.jmlr.org/papers/volume12/pedregosa11a/pedregosa11a.pdf) in the Shiu Lab's ML Pipeline. For a detailed overview of how to use this pipeline for genomic prediction or any kind of machine learning problem, please check out the [tutorial](https://github.com/ShiuLab/ML-Pipeline/tree/master/Tutorial). See the [Shiu Lab ML Pipeline](https://github.com/ShiuLab/ML-Pipeline) for more options for data-preprocessing and downstream analysis of ML regults.

For a summary of all options for the ML pipeline, run the command without any arguments. 

Example code
```
python ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVM -gs T -cv 5 -n 10 -tag soy_HT_SVRlin -save soy_HT_SVRlin 
python ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVMpoly -gs T -cv 5 -n 10 -tag soy_HT_SVRpoly -save soy_HT_SVRpoly
python ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVMrbf -gs T -cv 5 -n 10 -tag soy_HT_SVRrbf -save soy_HT_SVRrbf
python ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg RF -gs T -cv 5 -n 10 -tag soy_HT_RF -save soy_HT_RF
python ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg GTB -gs T -cv 5 -n 10 -tag soy_HT_GTB -save soy_HT_GTB
```

## Artificial Neural Network 

The ANN algorithms were implemented in python using [Tensorflow](https://www.tensorflow.org/) in the Shiu Lab's [ANN Pipeline](https://github.com/ShiuLab/ANN_Pipeline).


### Fully connected (Multi-layer perceptron) ANNs

python ANN_mlp.py -f full -x [Geno_file] -y [Pheno_file] -y_name [Y_Col_Name] -ho [test_set] -save [OUTPUT_NAME]

Example code:
```
# Run 10 replicates of the grid search using train/validation
for i in $(seq 1 10); do python ANN_Pipeline/ANN_mlp.py -f gs -gs_reps 1 -x soy_geno.csv -y soy_pheno.csv -sep , -y_name HT -ho soy_testset_1.txt -save ANN_soy_HT -epoch_thresh 0.001 -burnin 10; done

python ANN_Pipeline/ANN_mlp.py -f run -x soy_geno.csv -y soy_pheno.csv -sep , -n 10 -y_name HT -ho soy_testset_1.txt -save ANN_soy_HT -params ANN_soy_HT_GridSearch -tag ANN_soy_HT -epoch_thresh 0.001 -burnin 10 -s_yhat t -s_losses t 
```


### Fully connected (Multi-layer perceptron) ANNs with seeded starting weights

Available weighing methods: random, xavier (default), RF importance scores (RF), or coefficients from rrBLUP (rrB), BayesA, BayesB, Bayesian LASSO (BL), or Bayesian ridge regression (BRR)

python ANN_Pipeline/ANN_mlp.py -f full -x [Geno_file] -y [Pheno_file] -y_name [Y_Col_Name] -ho [test_set] -save [OUTPUT_NAME] weights [weighing_method]

Example code:
```
python ANN_Pipeline/ANN_mlp.py -f full -x soy_geno.csv -y soy_pheno.csv -sep , -n 10 -y_name HT -test soy_testset_1.txt -save ANN_soy_HT_seededRF -tag ANN_soy_HT_seedRF -epoch_thresh 0.001 -burnin 10 -s_yhat t -s_losses t -weights RF 
```


### Convolutional Neural Network 

The ANN algorithm was implemented in python using [Tensorflow v. 2.0](https://www.tensorflow.org/beta).

[MSU HPCC](https://wiki.hpcc.msu.edu/display/ITH/Installing+TensorFlow+2.0+Alpha+with+virtualenv) users can set up their own environment and run TF using the following environment:
```
module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load CUDA/10.0.130 cuDNN/7.5.0.56-CUDA-10.0.130
module load Python/3.6.4
source /mnt/home/azodichr/tf2env/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2
```

python CNN_regression.py -x [Geno_file] -y [Pheno_file] -y_name [Y_Col_Name] -test [test_set] -save [OUTPUT_NAME]

Example code for grid search:
```
# Run 100 replicates of the grid search using train/validation
for i in $(seq 1 10); do python CNN_regression.py -x rice_geno.csv -y rice_pheno.csv -sep , -x_sort alpha -y_name HT -test rice_testset_1.txt -save CNN_rice_HT -run f -gs_reps 10 -params gs -onehot t; done
```

We selected the best set of parameters on rice_HT_rep1 and use those for all other species/traits:
```
python CNN_regression.py -x soy_geno.csv -y soy_pheno.csv -sep , -x_sort alpha -y_name HT -test soy_testset_1.txt -save CNN_soy_HT -run t  -onehot t -cnn_type simple -filters 8 -kernel_l 16 -learn_rate 0.01 -pool_size 16 -activation relu -optimizer adam -run t -n 20 -x_sort alpha -onehot t 
```
or 
```
python CNN_regression.py -x soy_geno.csv -y soy_pheno.csv -sep , -x_sort alpha -y_name HT -test soy_testset_1.txt -save CNN_soy_HT -run t  -onehot t -params CNN_soy_HT1_GridSearch.txt -run t -n 20 -x_sort alpha -onehot t 
```



## Feature Selection

To select the best markers to include using the Feature_Selection.py script. For new versions of the Feature Selection code see the Shiu Lab's [ML_Pipeline](https://github.com/ShiuLab/ML-Pipeline/).

Example Code:
```
python Feature_Selection.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep , -f RF -n 500 -type r -test soy_testset_1.txt -save soy_test1_RF_500.txt
```

**How to use**

To run rrBLUP or a BGLR algorithm with only selected features, provide this output file inplace of 'all'. For example:
```
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv soy_test1_RF_500.txt HT soy_testset_1.txt soy rrblup_results/
Rscript GenomicSelection/predict_rrBLUP.R soy_geno.csv soy_pheno.csv soy_test1_RF_500.txt HT BL soy_testset_1.txt soy bl_results/
```

To run the ML_Pipeline or ANN_Pipeline using only the selected features, include the argument -feat. For example:
```
python ML-Pipeline/ML_regression.py -df soy_geno.csv -df2 soy_pheno.csv -y_name HT -sep ',' -test soy_testset_1.txt -alg SVM -gs T -cv 5 -n 10 -tag soy_HT_SVRlin -save soy_HT_SVRlin_featRF500 -feat soy_test1_RF_500.txt

python ANN_Pipeline/ANN_mlp.py -f full -x soy_geno.csv -y soy_pheno.csv -sep , -n 10 -y_name HT -test soy_testset_1.txt -save ANN_soy_HT_seededRF_featRF500 -tag ANN_soy_HT_seedRF -epoch_thresh 0.001 -burnin 10 -s_yhat t -s_losses t -weights RF -feat soy_test1_RF_500.txt 

```
