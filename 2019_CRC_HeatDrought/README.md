# The cis-regulatory codes of response to combined heat and drought stress in Arabidopsis thaliana

Authors: Christina B. Azodi, John P. Lloyd, Shin-Han Shiu.

## Abstract
Plants respond to their environment by dynamically modulating gene expression. A powerful approach for understanding how these responses are regulated is to integrate information about cis-regulatory elements (CREs) into models called cis-regulatory codes. Transcriptional response to combined stress is typically not the sum of the responses to the individual stresses. However, cis-regulatory codes underlying combined stress response have not been established. Here we modeled transcriptional response to single and combined heat and drought stress in Arabidopsis thaliana. We grouped genes by their pattern of response (independent, antagonistic, synergistic) and trained machine learning models to predict their response using putative CREs (pCREs) as features (median F-measure = 0.64). We then developed a deep learning approach to integrate additional omics information (sequence conservation, chromatin accessibility, histone modification) into our models, improving performance by 6.2%. While pCREs important for predicting independent and antagonistic responses tended to resemble binding motifs of transcription factors associated with heat and/or drought stress, important synergistic pCREs resembled binding motifs of transcription factors not known to be associated with stress. These findings demonstrate how in silico approaches can improve our understanding of the complex codes regulating response to combined stress and help us identify prime targets for future characterization. 


## Overview

- Step 1: Define single and combined heat and drought stress response groups
	- Process microarray data
	- Differential expression analysis
	- Defining groups
	- GO analysis
- Step 2: Identify known TFBMs and putative CREs for each response group
- Step 3:

## Step 1: Define single and combined heat and drought stress response groups

### Downloading and processing microarray data

Heat+Drought microarray expression data (Prasch and Sonnewald, 2013, Plant Phys) was downloaded (GSE46757) as log2-transformed, normalized to the 75th percentile and corrected to the median of all samples (Agilent-021169 Arabidopsis 4 Oligo Microarray (V4)). (see: https://github.com/azodichr/Combined_Stress_Response/blob/master/combine_Agilent_key.py)

Converted Agilent array probe IDs (https://earray.chem.agilent.com/earray/, Design ID 021169) to TAIR10 gene names using IDswop, if multiple probes overlapped the same gene, IDswop takes the mean or tosses the gene out if the values > 20% different.

### Differential Expression analysis

Used edgeR to calculate differential expression using 3 contrasts:
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
python ~/GitHub/MotifDiscovery/pCRE_Finding_FET.py -k 6mers.txt -neg NNN_genes.txt.fa -pos NNU_genes.txt.fa -save NNU_pcres_FET01
```

### Enriched TFBMs

**DAP-Seq data**

Master data frame of all gene’s promoters and which contain a DAP-peak:
- peaks_parsed_ampRemoved.csv.mapped.regulator.matrix.altsplicRemoved.txt .

**Cis-BP TFBMs**

A. CISbp_to_include.txt: Removed all that were already included in DAP-Seq dataset. Then removed those that only had inferred evidence (leaving only those with direct evidence). Finally, if there were still duplicates of a TF I kept Weirauch over JASPAR or Transfac unless the JASPAR data was from ChIP-chip method. Then I kept JASPAR (2014) over Transfac (2006). If there were 2 motifs from Weirauch I kept the one with more motif info (longer, higher certainty).

B. Split cluster fasta files of all non-overlapping promoters into individual files:
```python ~shius/codes/FastaManager.py -f indiv -fasta TAIR10_Non_Overlapping_Promoters.fa -odir 01_ProcessingCISbp/```

C. Make and run mapping jobs and split into 18 files with 1000 jobs in each file:
```
python ~/GitHub/Combined_Stress_Response/make_mapping_runcc.py -fasta 01_ProcessingCISbp/ -pwm /mnt/home/azodichr/CIS_BP/pwms_all_motifs/00_pwm_format/ -include CISbp_to_include.txt
python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c runcc_pwm_mapping -w 90 -m 4 -n 1000 -wd /mnt/home/azodichr/01_CombinedStress/
```


$ head -1 AT1G01020.fa_hits.txt > all_hits.txt
$ tail -n +2 -q AT*_hits.txt >> all_hits.txt

Move all individual hit files & runcc files to: /mnt/home/azodichr/01_CombinedStress/01_ProcessingCISbp/00_At_nonOverlapPromoters_Hits/

D. Convert hits results into ML like data frame
$ python ~/GitHub/Combined_Stress_Response/convert_to_df.py -mapping all_hits.txt

**Test for enrichement:** Using Feature Selection code from ML Pipeline (https://github.com/ShiuLab/ML-Pipeline/blob/master/Feature_Selection.py)

```
python ML-Pipeline/Feature_Selection.py -df knownDF_NNU_df.txt -alg FET -p 0.01 -list t
```



### Compare classifiers:
```
python ~/GitHub/ML-Pipeline/compare_classifiers.py -ids pCREs,KnownFET -save NNU_pCRE_vs_Known -scores NNU_Ras_FET01_df_p0.01.txt_RF_pCREs_scores.txt,knownDF_NNU_df.txt_RF_KnownFET_scores.txt
```
