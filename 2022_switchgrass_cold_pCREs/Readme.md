# Temporal Regulation of Cold Transcriptional Response in Switchgrass
**Authors: Thilanka Ranaweera, Brianna N.I. Brown, Peipei Wang, Shin-Han Shiu**

## Abstract
In switchgrass, low-land ecotypes have significantly higher biomass but lower cold tolerance compared to high-land ecotypes. Understanding the molecular mechanisms underlying cold response, including the ones at transcriptional level, is important for improving tolerance of high-yield switchgrass under chilling and freezing environmental conditions. Although the cold transcriptional responses in plants are well studied, it remains unclear how temporal cold responses are regulated by regulatory switches, known as cis-regulatory elements (CREs). Here, by investigating the transcriptional responses of switchgrass genes over a time course after cold treatment, we found that the number of cold-responsive genes and enriched Gene Ontology terms increased, suggesting an amplified response/cascading effect in gene expression. Machine learning models built using enriched k-mers were able to distinguish the cold-responsive genes from non-responsive genes, suggesting that these k-mers (or putative cis-regulatory elements, pCREs) may be regulators of cold-response in switchgrass. Our models allowed identification of known cold-responsive CREs such as the C-repeat/dehydration response element, W-box, CGCG-box, and GCC box. Performance comparison between ML models built using k-mers and experimentally validated CREs indicated that known CREs alone cannot explain cold stress response in switchgrass. In addition, the bulk of pCREs identified among different timepoint models were similar in sequence, indicating their regulatory role for genes that show in complex expression patterns. Our study suggested a presence of well-known cold responsive pathways in switchgrass, i.e., the CBF-dependent pathways. Moreover we showcase how well the computational approaches can be utilized to uncover novel condition-specific regulatory switches and their regulatory patterns.

*Key words*: Regulation of cold stress, machine learning, switchgrass, CBF-dependent pathways

----
## Overview
#### 1. Data collection and preprocessing
#### 2. Identification of cold-responsive putative cis-regulatory elements
#### 3. Selection of minimal pCRE sets as features and determining relationships between pCREs
#### 4. Identification of transcription factors (TFs) with binding sites similar to pCREs

## Software requirements
#### Machine Learning Pipeline
- biopython v1.73
- matplotlib v1.5.1
- numpy v1.16.2
- pandas v0.24.2
- python v3.4.4
- scikit-learn v0.20.3
- scipy v1.2.1

#### Motif Mapping
- [motility](https://github.com/ctb/motility)
- TAMO v1.0
- R v3.5.3
- agilp v3.8.0

#### Other
- add add add
---
### 1. Data collection and preprocessing
- The sequenced transcriptomes related with switchgrass stress responses were collected from published studies
  - Switchgrass cold stress data: [Meng et al, (2021)](https://www.pnas.org/doi/10.1073/pnas.2026330118)
  - Switchgrass dehydration stress data: [Zhang et al, (2018)](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-018-1088-x)
  - Switchgrass salt stress data [Zhang et al, (2021)](https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-018-1088-x)
  - Switchgrass drought stress data [Zuo et al, (2018)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204426)

- The RNA-seq data were processed using [RNA-seq_data_processing](https://github.com/ShiuLab/RNA-seq_data_processing.git)

### 2. Identification of cold-responsive putative cis-regulatory elements
##### &emsp; a. GO enrichment analysis
&emsp;&emsp;&emsp;&emsp; - **GO_and_Pathway_enrichment_analysis.r** was used to find BP GO terms enriched in set of genes.
 &emsp;&emsp;&emsp;&emsp; - **Heatmap_for_GO_enrichment.R** was used to make the final heatmap of enriced GO terms.
 ##### &emsp; b. Identifying pCREs responsible for cold stress
 &emsp;&emsp;&emsp;&emsp; - pCREs responsible for cold stress regulation were identified using [pCRE_identification](https://github.com/peipeiwang6/pCRE_identification.git) pipeline.
### 3. Selection of minimal pCRE sets as features and determining relationships between pCREs
##### &emsp; a. Feature selection
  * To automate feature selction **do_feature_selection_mod.py** was used
  This script uses the feature importance file to pick top n features and run models recursively
```
python do_feature_selection_mod.py train_enriched_kmers_matrix.txt_RF_cv5_score_roc_auc_y_imp
```
  * To run model with specific number of features **do_feature_selection_with_specific_number_fets.py** was used.
  within the code you need to specify how many features you want to run the model with.
  * To output a file with number of features and F1 performance **get_CV_score.py** was used. To run this script it should be in the  same directory the feature selection as done.

##### &emsp; b.Determining the similarity between pCREs  
&emsp;&emsp;&emsp;&emsp; I. Making list of kmers with only forward sequences (Because we used forward and reverse complementary sequnces)
```
sed "s/\..*//" top_4_kmer.txt > top_4_kmers_no_RC.txt
```
&emsp;&emsp;&emsp;&emsp; II. Generating TAMO file for kmer list
Using Shiu lab [Motif discovery pipeline's](https://github.com/ShiuLab/MotifDiscovery.git) **generate_PWM.py**

&emsp;&emsp;&emsp;&emsp; III. Making distance matrix based on PCC distances
Used **3.calculate_distance_matrix.py** in [CRE-Pipeline](https://github.com/ShiuLab/CRE-Pipeline.git)

&emsp;&emsp;&emsp;&emsp; IV. Creating a tree and getting cultures of similar pCREs
Used **getting_kmer_similarities.R** to obtain a distances tree and clusters based on similarity threshould of 0.39. Also similar script was used to calculate the propotion of general and time specific pCREs in a pCRE clusters.

### 4.Identification of transcription factors (TFs) with binding sites similar to pCREs

##### &emsp; a. Assessing the sequence similarity between pCREs and known transcription factor binding sites (TFBSs)
* To identify pCREs similar to known TFBS scripts in [Motif Discovery Pipeline]((https://github.com/ShiuLab/MotifDiscovery.git) and the instructions provided in [2019_CRC_HeatDrought](https://github.com/ShiuLab/Manuscript_Code/tree/master/2019_CRC_HeatDrought) repository were used.
*  pCREs similar to known TFBS identification was completed in batches for pCRE lists using **batch_merge_results.py** and **batch_get_sig_TF.py** scripts.

##### &emsp; b. Assessing how well the binding sites of TFs known to regulate cold response might predict cold response

&emsp;&emsp;&emsp;&emsp; I. Converting meme formatted dapseq pwms and CISBP pwms for TFs into pwms that can be used in C.Azodi mapping pipeline was done using **convert_dap_pwm.py** and **convert_CISBP_pwms.py** respectively.
&emsp;&emsp;&emsp;&emsp; II. Since mapping with large number of sequnces (> 100) is time consuming we binned all upregulated genes in to ~900 bins. It was done using **binning_genes.R** script.

&emsp;&emsp;&emsp;&emsp; III. Binned genes were mapped in batches using **batch_mapper.py**

&emsp;&emsp;&emsp;&emsp; IV. **check_results.py** and **clening_results.py** scripts were used to check mapped results and clean unmapped or error files generated using mapping.

&emsp;&emsp;&emsp;&emsp; V. Feature table was generated using mapping results using **batch_feat_tbl.py**.

&emsp;&emsp;&emsp;&emsp; VI. **get_all_feature_tables.py** was used to pick fetaure tables made in the batch analysis.

&emsp;&emsp;&emsp;&emsp; VII. Final feature table was created for all the upregulated genes based on mapping of known TFBSs using **dict_tf.py**

&emsp;&emsp;&emsp;&emsp; VIII. **get_new_feature_tbl.py** script was used to take instances and lables from feature tables of initial timepoint models (because we wanted to keep same training and testing instances as initial timepoint models) and merge with new features from mapping.

&emsp;&emsp;&emsp;&emsp; IX. **make_models.py** was used to run
