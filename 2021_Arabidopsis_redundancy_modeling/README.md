# Predictive models of genetic redundancy in *Arabidopsis thaliana*


## Files
- Feature matrices for the training datasets primarily discussed in this paper are [available on Zenodo](https://doi.org/10.5281/zenodo.3987384) as Dataset_1.txt and Dataset_2.txt
- Feature matrices for predicting genetic redundancy of selected additional gene pairs are [also available on Zenodo](https://doi.org/10.5281/zenodo.3987384) as Dataset_3.txt and Dataset_4.txt
- [ML classification script](scripts/ML_classification.py), go [here](https://github.com/ShiuLab/ML-Pipeline) for the most recent version
- [feature selection script](scripts/Feature_Selection.py), go [here](https://github.com/ShiuLab/ML-Pipeline) for the most recent version
- [script for selection of test set instances](scripts/test_set.py), go [here](https://github.com/ShiuLab/ML-Pipeline) for the most recent version  

## Overview

### Get test sets

#### Input
- feature matrix with all training instances (e.g. RD9_full_feature_matrix.txt, available [here](https://doi.org/10.5281/zenodo.3987384) as Dataset_2.txt)

#### Example
```
python test_set.py -df RD9_full_feature_matrix.txt -type c -p 0.1
```

#### Output
- file containing names of instances selected as test set (e.g. RD9_test_set.txt, available [here](../data/))


### Feature selection

#### Input
- feature matrix with all training instances (e.g. RD9_full_feature_matrix.txt)
- test set file (e.g. RD9_test_set.txt)

#### Example
```
python Feature_Selection.py -df RD9_full_feature_matrix.txt -test RD9_test_set.txt -alg RF -n 1000 -scores True 
```

#### Output
- selected feature file (e.g. selected_features_RD9_RF_BT_200.txt, available [here](../data/))


### Model building

#### Input
- feature matrix with all training instances (e.g. RD9_full_feature_matrix.txt)
- test set file (e.g. RD9_test_set.txt)
- selected feature file (e.g. selected_features_RD9_RF_BT_200.txt)

#### Example
```
python ML_classification.py -df RD9_full_feature_matrix.txt -alg SVM -cl_train positive,negative -test RD9_test_set.txt  -gs True -cm True -plots True -feat selected_features_RD9_RF_BT_200.txt -x_norm force_false
```
