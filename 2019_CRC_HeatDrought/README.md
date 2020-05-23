# The cis-regulatory codes of response to combined heat and drought stress in Arabidopsis thaliana

Authors: Christina B. Azodi, John P. Lloyd, Shin-Han Shiu.

## Abstract
Plants respond to their environment by dynamically modulating gene expression. A powerful approach for understanding how these responses are regulated is to integrate information about cis-regulatory elements (CREs) into models called cis-regulatory codes. Transcriptional response to combined stress is typically not the sum of the responses to the individual stresses. However, cis-regulatory codes underlying combined stress response have not been established. Here we modeled transcriptional response to single and combined heat and drought stress in Arabidopsis thaliana. We grouped genes by their pattern of response (independent, antagonistic, synergistic) and trained machine learning models to predict their response using putative CREs (pCREs) as features (median F-measure = 0.64). We then developed a deep learning approach to integrate additional omics information (sequence conservation, chromatin accessibility, histone modification) into our models, improving performance by 6.2%. While pCREs important for predicting independent and antagonistic responses tended to resemble binding motifs of transcription factors associated with heat and/or drought stress, important synergistic pCREs resembled binding motifs of transcription factors not known to be associated with stress. These findings demonstrate how in silico approaches can improve our understanding of the complex codes regulating response to combined stress and help us identify prime targets for future characterization. 


## Overview

- Step 1: Define single and combined heat and drought stress response groups
	- Downloading and processing microarray data
	- Differential expression analysis
	- Defining response groups
	- GO analysis
- Step 2: Identify known TFBMs and putative CREs for each response group
	- pCRE finding
	- Enriched TFBMs
	- Test for enrichment
- Step 3: Machine Learning with Random Forest
- Step 4: Convolutional Neural Network (multi-omic model)
- Step 5: Get best matching known TFBMs for pCREs 

## Software requirements

Machine Learning Pipeline (aka Random Forest):

	- biopython v1.73
	- matplotlib v1.5.1
	- numpy v1.16.2
	- pandas v0.24.2
	- python v3.4.4
	- scikit-learn v0.20.3
	- scipy v1.2.1

Convolutional Neural Networks (CNNs):

	- Tensorflow v2.0 
	- Keras 2.2.4
	- Python v3.7.3 
	- h5py 2.9.0
	- scipy v1.3.0
	- pandas v0.24.2
	- numpy v1.16.4
	- scikit-learn v0.21.2

Other:

	- R v3.5.3
	- agilp v3.8.0
	- limma v3.38.3
	- [motility](https://github.com/ctb/motility)
	- TAMO v1.0
	- agricolae v1.3.1


Experiments in this study were conducted on [MSU's High Performance Computing Cluster](https://icer.msu.edu/hpcc/hardware)

## Data

Example datasets are available in the /data/ directory:

**example_data/NNU_pcres_df.txt**: Rows are genes, first column is the class (i.e. NNU=1 or NNN=0), and the remaining columns are the pCRE features (1=present in promoter of gene, 0=not present)

**example_data/NNU_pcres_MultiOmic.txt**: Same general format as pCRE only data, except each additional omic data has it's own column named like: pcreA_omic1, pcreA_omic2, etc. Note that missing omic data is not allowed (i.e. must have the same number of columns for each pCRE).

Full datasets are available on [Zenodo](https://zenodo.org/record/3840714#.XsiT7xMzY3E).

## Step 1: Define single and combined heat and drought stress response groups

### Downloading and processing microarray data

Heat+Drought microarray expression data (Prasch and Sonnewald, 2013, Plant Phys) was downloaded (GSE46757) as log2-transformed, normalized to the 75th percentile and corrected to the median of all samples (Agilent-021169 Arabidopsis 4 Oligo Microarray (V4)). (see: https://github.com/azodichr/Combined_Stress_Response/blob/master/combine_Agilent_key.py)

Converted Agilent array probe IDs (https://earray.chem.agilent.com/earray/, Design ID 021169) to TAIR10 gene names using IDswop from agilp (v3.8.0), if multiple probes overlapped the same gene, IDswop takes the mean or tosses the gene out if the values > 20% different.

### Differential Expression analysis

Used limma v3.38.3 to calculate differential expression using 3 contrasts:
	- heat vs. control
	- drought vs. control
	- heat & drought vs. control

### Defining response groups

Identified response groups using adj.p.value <=0.05 and logFC >= or <= 1 for each pattern. For example (NNU):

```
nnu <- subset(sigDGEs, Drought > -1 & Drought < 1 & Heat > -1 & Heat < 1 & Combined >= 1 & Dual_adj.P.Val <= 0.05)
```

Then, filter to only include genes with non-overlapping promoters (based on TAIR10 gff file).


### GO Analysis

Instructions: https://github.com/ShiuLab/GO-term-enrichment



## Step 2: Identify known TFBMs and putative cis-regulatory elements

For more details see: https://github.com/azodichr/MotifDiscovery

### pCRE finding

Use TAIR10 to get gene sequences for genes in each response group. Then use "pCRE_Finding_FET.py" to iteratively find enriched k-mers.

```
python ~/GitHub/Utilities/FastaManager.py -f NNU_genes.txt -fasta TAIR10_Non_Overlapping_Promoters.fa
python ~/GitHub/Utilities/FastaManager.py -f NNN_genes.txt -fasta TAIR10_Non_Overlapping_Promoters.fa
python ~/GitHub/MotifDiscovery/pCRE_Finding_FET.py -k 6mers.txt -neg NNN_genes.txt.fa -pos NNU_genes.txt.fa -save NNU_pcres_FET01
```

### Enriched TFBMs

**DAP-Seq data**

Data frame of the presence/absence of DAP-Seq peaks in promoters (1 kb) of A. thaliana genes from TAIR10.
- peaks_parsed_ampRemoved.csv.mapped.regulator.matrix.altsplicRemoved.txt .

**Cis-BP TFBMs**

CISbp_to_include.txt: Removed all that were already included in DAP-Seq dataset. Then removed those that only had inferred evidence (leaving only those with direct evidence). Finally, if there were still duplicates of a TF I kept Weirauch over JASPAR or Transfac unless the JASPAR data was from ChIP-chip method. Then I kept JASPAR (2014) over Transfac (2006). If there were 2 motifs from Weirauch I kept the one with more motif info (longer, higher certainty). Then map motifs to TAIR10 and make a presence/absence dataframe like for DAP-Seq

```
# Get separate fasta files for each non-overlapping promoter:
python ~GitHub/ShiuLab/Utilities/FastaManager.py -f indiv -fasta TAIR10_Non_Overlapping_Promoters.fa -odir 01_ProcessingCISbp/

# Make mapping jobs and run:
python scripts/make_mapping_runcc.py -fasta 01_ProcessingCISbp/ -pwm /mnt/home/azodichr/CIS_BP/pwms_all_motifs/00_pwm_format/ -include CISbp_to_include.txt
python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c runcc_pwm_mapping -w 10 -m 4 -n 1000
python scripts/convert_to_df.py -mapping all_hits.txt
CISBP_Hits2NonOverlappingPromoters_AT.matrix.txt
```

**Test for enrichment** 

Using Feature Selection code from ML Pipeline (https://github.com/ShiuLab/ML-Pipeline/blob/master/Feature_Selection.py)

```
python ~/GitHub/ML-Pipeline/Feature_Selection.py -df knownDF_NNU_df.txt -alg FET -p 0.01 -list t
```

## Step 3: Machine learning (pCREs only)

For an extensive tutorial on the ML pipeline and the most recent version of the ML Pipeline, [see here](https://github.com/ShiuLab/ML-Pipeline/).

**Example of ML**
```
python ~/GitHub/ML-Pipeline/ML_classification.py -df example_data/NNU_pcres_df.txt -gs T -alg RF -plots F -n 100 -p 5 -threshold_test auroc -tag pCREs

# Also run with KnownTFBM and pCRE+KnownTFBM features
# Compare the classifications from each model:
python ~/GitHub/ML-Pipeline/scripts_PostAnalysis/compare_classifiers.py -ids pCREs,KnownTFBM -save NNU_pCRE_vs_Known -scores NNU_pcres_df.txt_RF_scores.txt,NNU_knownTFBM_df.txt_RF_scores.txt,NNU_combined_df.txt_RF_scores.txt
```


## Step 4: Convolutional Neural Network (multi-omic model)

CNN models were run on MSU's [HPC](https://icer.msu.edu/hpcc/hardware) in a virtual environment installed as follows:

```
module load GCC/6.4.0-2.28  OpenMPI/2.1.2
module load CUDA/10.0.130 cuDNN/7.5.0.56-CUDA-10.0.130
module load Python/3.6.4
virtualenv -p python3 tf2env
source tf2env/bin/activate
pip install --upgrade --force-reinstall tensorflow-gpu==2.0.0-alpha0 
```

Example of implemeting CNN code.
```
source tf2env/bin/activate

# Convert dataframe (see example) into format needed for CNN
python scripts/input_converter.py -input example_data/NNU_pcres_MultiOmic.txt -out NNU -method omic_stack -y_name Class

# Train and test CNN
python scripts/CNN_TF2_omic.py -x NNU_x.npy -y NNU_y.npy -save NNU -n 100 -run t -params NNU_GridSearch.txt -imp_m t -imp_k t
```

*Note: for our largest response group, the parameter search required 40 GPU hours and training required <2 GPU hours.*

## Step 5: Get best matching known TFBMs for pCREs 

Find scripts here: https://github.com/ShiuLab/CRE-Pipeline

*Software required: TAMO v1.0, python2*
```
# Generate tamo file
python ~/GitHub/CRE-Pipeline/generate_PWM.py pcres_all.txt.tm

## Get distance matrices
# Convert CIS-BP motifs into tamo file: Athaliana_TFBM_v1.01.tm.index.tm
python ~/GitHub/CRE-Pipeline/pcc_merge_CC.py create_cc_2 -t pcres_all.txt.tm -t2 Athaliana_TFBM_v1.01.tm.index.tm
bash runcc

# Convert DAP-Seq motifs into tamo file: DAP_motifs.txt.tm
python ~/GitHub/CRE-Pipeline/pcc_merge_CC.py create_cc_2 -t pcres_all.txt.tm -t2 DAP_motifs.txt.tm
bash runcc

# Merge results
python ~/GitHub/CRE-Pipeline/pcc_merge_CC.py combine_distance_matrix_2 -t pcres_all.txt.tm -t2 Athaliana_TFBM_v1.01.tm.index.tm
python ~/GitHub/CRE-Pipeline/pcc_merge_CC.py combine_distance_matrix_2 -t pcres_all.txt.tm -t2 DAP_motifs.txt.tm

# Get best matches
python scripts/motif_pcc_best_match.py -t pcres_all.txt.tm -pcc_cis pcres_all.txt.tm-Athaliana_TFBM_v1.01.tm.index.tm.dm -pcc_dap pcres_all.txt.tm-DAP_motifs.txt.tm.dm -save pcre

python scripts/TFBM_pcc_significance.py -top pcre_TopHits.txt
```
