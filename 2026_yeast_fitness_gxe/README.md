# Predictive models of the genetic bases underlying budding yeast fitness in multiple environments
Manuscript code and datasets

---

### Table of Contents
- [Datasets](#datasets)
- [Data-processing](#data-processing)
- [Models](#models)
- [Generate figures and tables](#generate-figures-and-tables)

---

## Datasets
- Raw datasets from Peter et al. 2018 can be downloaded from http://1002genomes.u-strasbg.fr/files/
- Supplementary files generated for this project were uploaded to Zenodo at https://doi.org/10.5281/zenodo.17245961

---

## Data-processing
This section describes how to generate the processed feature tables (Files S1, S2, S4, S5, and S6) used in the modeling steps.
| Data Type   | Script                      | Description                                                                                        |
| ----------- | --------------------------- | -------------------------------------------------------------------------------------------------- |
| Fitness     | `pre-process_pheno_data.py` | Filter and clean the raw phenotype dataset.                                                        |
| SNPs        | `pre-process_snp_data.sh`   | Pre-process the single nucleotide polymorphism (`SNP`) genotypes and build the SNP feature matrix. |
| PCs         | `get_PCs.py`                | Perform a principal component analysis on the SNP feature matrix (`PC`).                           |
| PAVs & CNVs | `pre-process_orf_data.py`   | Pre-process the presence/absence (`PAV`) and copy number variation (`CNV`) genotype data.          |

---

## Models 
### Conda environment requirements
- For training RF models, see the package versions in `ml-pipeline_environment.yml`.
- For training XGBoost models, see the `shap_environment.yml` file. This conda environment was also utilized for estimating SHAP values of RF models and conducting the analysis described in [Generate figures and tables](#generate-figures-and-tables).

### Scripts by algorithm
| Script                                | Description                                                                     |
| ------------------------------------- | ------------------------------------------------------------------------------- |
| `BL_Model.R`                          | Implements the Bayesian LASSO (BL) model.                                       |
| `BayesC_Model.R`                      | Implements the BayesC model.                                                    |
| `XGB_Model.py` and `xgboost_model.py` | Implement the Extreme Gradient Boosting (XGBoost) model.                        |
| `rrBLUP_Model.R`                      | Implements the Ridge Regression Best Linear Unbiased Prediction (rrBLUP) model. |

The Random Forest regression models were trained using the pipeline available at https://github.com/ShiuLab/ML-Pipeline

### SLURM batch scripts for training models
These scripts follow a general naming convention {alg}\_{feature table type}\_{genetic variant type}.sb, where:
- `{alg}` – algorithm used (`RF`, `XGB`, `rrBLUP`, `BayesC`, `BL`) to train models.
- `{feature table type}` – type of feature table:
  - `complete`: models trained using all features in the input dataset.
  - `optimized`: models trained using optimized or feature-selected subsets.
- `{genetic variant type}` – genetic variant dataset type (`PC` [first 5 principal components of the SNP data], `SNP`, `PAV`, `CNV`).

### Test set isolates
The test set isolates which were used to evaluate the models are provided in the file `Test.txt`.

### Input feature table example format for training models
Phenotype matrix (S1_File.xlsx):
<img width="410" height="136" alt="image" src="https://github.com/user-attachments/assets/4d1a9ef5-8d14-465c-8891-0b262cf52f03" />

SNP matrix (S2_File.csv):
<img width="482" height="141" alt="image" src="https://github.com/user-attachments/assets/2fa9aa8a-df4d-4e8d-89fc-6593df95f2bc" />

Example .txt file containing a list of features (one feature per line) for subsetting the target features for model training:
<img width="140" height="185" alt="image" src="https://github.com/user-attachments/assets/d78f9ac0-a227-4525-9449-001f232cbff5" />

Example CVF file (`CVFs_20rep.csv`; cross-validation fold assignments file for running rrBLUP, BL, and BayesC models):

This file was generated using the `07_make_CVs.py` script available at https://github.com/ShiuLab/Manuscript_Code/blob/master/2022_GP_in_Switchgrass/07_make_CVs.py
<img width="673" height="110" alt="image" src="https://github.com/user-attachments/assets/d365862a-efaf-4bd4-8ee9-859f343a0dab" />


### Model interpretation
| Script                                | Description                                                                                  |
| ------------------------------------- | ---------------------------------------------------------------------------------------------|
| `SHAP_training_saving_interaction.py` | Calculate SHAP values and SHAP interaction values using a trained RF model.                  |
| `Sum_SHAP_interaction.py`             | Sum the SHAP interaction values across instances.                                            |
| `Calculate_RF_SHAP_values.sb`         | Calculate the SHAP values for all RF models built on the complete or optimized feature sets. |

---

## Generate figures and tables
This directory contains scripts used to generate all main figures, supplementary figures, supplementary tables, and the remaining supplementary files.

| Script                                                        | Description                                                                                                 |
| ------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `Figure_1A_heatmaps.R`                                        | Generates figure 1A.                                                                                        |
| `Figure_1B-E_data_corr.py`                                    | Generates figures 1B, 1C, 1D, and 1E.                                                                       |
| `Figure_2_model_performance.py`                               | Generates figure 2.                                                                                         |
| `Figure_2b–g_S3_fitness_distributions.py`                     | Generates figures 2B, 2C, 2D, 2E, 2F, 2G and the supplementary figure S3.                                   |
| `Figure_3a_S4a–e_compare_imp_measures.py`                     | Generates figure 3A and supplementary figures S4A–E.                                                        |
| `Figure_3b_S4f_shared_genes_across_envs.py`                   | Generates figure 3B and supplementary figure S4F.                                                           |
| `Figure_4a–c_TableS12_SFile11_benchmark_gene_models.py`       | Generates figures 4A, 4B, 4C, supplementary tables S12, and supplementary file S11.                         |
| `Figure_4a–c_benchmark_gene_models.sb`                        | SLURM script for submitting models shown in figure 4A, 4B, and 4C.                                          |
| `Figure_4d–e_S5_TableS13_SFile13–14_lab_strain_similarity.py` | Generates figures 4D and 4e and supplementary tables S13 and S14.                                           |
| `Figure_5_S6–10_TableS14_SFile15–16_SHAP_clustering.py`       | Generates figures 5, S6–10, supplementary table S14, and supplementary files S15 and S16.                   |
| `Figure_6_S11–12_TableS16–19_SHAP_interactions.py`            | Generates figures 6, S11, and S12 and supplementary tables S16–19.                                          |
| `Figure_6_S11–12_TableS16–19_SHAP_interactions.sb`            | SLURM batch script for figure 6 models.                                                                     |
| `Figure_S1_data_pair_corr.py`                                 | Computes data pair correlations for Supplementary Figure S1.                                                |
| `Table_S3_performance_comparisons.py`                         | Generates supplementary table S3.                                                                           |
| `Table_S4_factors_regression.py`                              | Generates supplementary table S4.                                                                           |
| `Table_S5_S9_SFile7_gini_vs_shap_rho.py`                      | Generates supplementary tables S5 and S9 and supplementary file S7.                                         |
| `Table_S6_rank_rho_btwn_variants.py`                          | Generates supplementary table S6.                                                                           |
| `Table_S7_rank_per_corr_btwn_envs.py`                         | Generates supplementary table S7.                                                                           |
| `Table_S8_ORF_go_enrichment.R`                                | GO enrichment results of the PAV and CNV RF models saved in the go enrichment sheet of table S8.            |
| `Table_S8_ORF_pwy_enrichment.R`                               | Pathway wanrichment results of the PAV and CNV RF models saved in the pathway_enrichment sheet of table S8. |
| `Table_S8_SNP_go_enrichment.R`                                | GO enrichment results of the SNP RF models saved in the go enrichment sheet of table S8.                    |
| `Table_S8_SNP_pwy_enrichment.R`                               | Pathway wanrichment results of the SNP RF models saved in the "pathway_enrichment" sheet of table S8.       |
| `Table_S8_combine_go_pwy_enrichment.py`                       | Combines GO and pathway enrichment outputs into supplementary table S8.                                     |
| `Table_S11_lit_gene_enrichment.py`                            | Generates supplementary table S11.                                                                          |
| `Calculate_heritability_sommer.R`                             | Estimates heritability using the `sommer` package.                                                          |
| `S12_S14_S17_Table_additional_analyses.py`                    | Comparative analysis of SHAP values and SHAP interaction values. Results added to S12, S14, and S17 Tables. |
| `Table_S1_genetic_variation.py`                               | Estimates of genetic variation using SNPs, PAVs, and CNVs. Results added to Table S1.                       |

