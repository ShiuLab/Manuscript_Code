# Spatial Expression Evolution of Post Whole Genome Duplicates in Cultivated Tomato


**This is the code block used to analyze the dataset in different section of the manuscript**

## ABSTRACT
***Large-scale duplication events, such as whole genome duplications (WGD), are key drivers of plant genome evolution, enabling functional diversification and innovation. In cultivated tomato, the most recent whole genome triplication (WGT) has led to the retention of relatively young paralogous genes, offering a unique opportunity to study the evolutionary dynamics of gene expression following such events. In this study we aim to explore expression evolution of larger sets of tissue types through comparing extant expression profiles as well as reconstructed ancestral states of post WGT paralogous genes. Our analysis of expression patterns across 23 tissue types and ancestral state reconstruction in post-WGT paralogs reveals a high frequency of expression divergence, with ~6% of gene branches showing inferred gain events compared to ~20% with losses. The genes that gained expression were significantly enriched with the functions that lead to traits that favors the domestication of tomato.  Notably, asymmetrical expression divergence was prevalent among retained post-WGT genes, with expression retention more common when ancestral genes were either non-expressed or highly expressed. We also observed strong organ-specific expression evolution, particularly tissue types originating from anatomical regions or following a sequential developmental process suggesting a tight spatial specific expression regulation. Ancestral state analysis of tissue expression dependence shows that organ-specific expression is influenced by expression loss in other tissue types. Syntenic gene block analysis further highlighted independent evolutionary trajectories among genes retained after the WGT event.***

***Key words: Ancestral State Reconstruction, Expression Divergence, Organ-Specific Expression, Post-WGD Paralogs, Solanaceae-Specific Genome Triplication, Spatial Dependence of Expression, Syntenic Gene Blocks***

## MATERIALS AND METHODS

### 1. Expression data collection, processing, and determining expression levels and differential expression

&emsp;&emsp;**Analysis in this section was mainly carried out using RNA-seq data processing pipeline (https://github.com/ShiuLab/RNA-seq_data_processing.git)**

### 2. Determining absolute expression threshold

&emsp;&emsp;Mixed modeling of TPM from all the tissue types was carried out by following code
```
01_Mixed_models.R
```
### 3. Determining the paralogous genes evolved from Solanacea specific whole genome triplication event

&emsp;&emsp;i. First, to identify syntenic blocks retained after duplication events, we employed the MCScan2 algorithm.

&emsp;&emsp;ii. Next to distinguish syntenic blocks arising from different duplication events, the median Ks distribution of the syntenic blocks was modeled using a mixture modeling approach . 

```
02_Get_KS_plots.R
```

###  4. Identifying the extant expression divergence of the post WGT paralogous gene pairs

&emsp;&emsp; Extant expression divergence and simulation for random expectation was does using below code.

```
03_Looking_at_extant_expression_divergence.R
```
&emsp;&emsp; Next we looked at the extreme extant expression divergence scenarios 

```
04_Observing_the_extreme_divergence_sen.R
```
&emsp;&emsp; All the GO term enrichment analysis during this study was performed using below code.

```
05_GO_enrichment_on_different_gene_sets.R
```

### 5. Determining  the spatial dependence of extant expression of paralogous gene pairs

&emsp;&emsp; To assess the spatial dependence of paralogous gene expression, pairwise tissue comparisons were performed to determine the expression association between paralogous gene pairs across different tissue types using below code.

```
06_Tissue_partition_enrichment.R
```

### 6. Phylogenetic Reconstruction

&emsp;&emsp; i. To reconstruct the phylogenies first out group for post WGT paralogoues genes were determined from a BLAST search. The γ duplicates showing the highest sequence similarity to the T syntenic paralog–tandem duplicate genes were selected as potential outgroups.

```
07_Iteratively_pick_OG_based_of_BLAST.ipynb
```
&emsp;&emsp; ii. CDS FASTA files were then used to construct phylogenies 
```
08_Get_fasta_from_guide_convert_to_AA.py
```
&emsp;&emsp; iii. Phylogenies were constructed using RAxML. The itarative jobs were submitted to MSU HPCC using below code

```
09_Submit_jobs_to_get_alignmnts_and_trees.py
```
&emsp;&emsp; iv. Incomplete tree were iteratively searched using below code.

```
10_Check_for_incomplete_trees_runs.py
```

&emsp;&emsp; v. Since BLAST was used for finding OGs sometimes RAxML reconstrction of the tree might not place OG outside the ingroup. If the OG placement was wrong, those trees were chosen using below code

```
11_Check_ML_tree_topology.py
```
&emsp;&emsp; vi. For the trees with wrong OG placement below code was used to repick a OG

```
12_Re_pick_OGs_after_ML_run.py
```
### 7. The ancestral state reconstruction for expression levels

***The ancestral expression levels were reconstructed using <Github_Repo>***

### 8.Analyzing the ancestral expression divergence among extant paralogs

&emsp;&emsp; To test how ancestral gene expression was diverged among post WGT paralogs, we compared the reconstructed ancestral expression states to the extant expression states using below code. Both expression divergence on single branches and paralogous gene pairs are assessed using the same code.

```
13_Analyzing_ancestral_spatial_exp_divergence.R
```

### 9.Determining the ancestral spatial expression association of syntenic paralogs

&emsp;&emsp; To evaluate the spatial dependence in expression evolution of ancestral genes that arose from  WGT events below code was used.

```
14_Analysing_ancestral spatial expression association.py
```

### 10. Calculating block level expression partitioning

&emsp;&emsp; The Block level partition score based of foldchange of absolute expression was calculated using below code.

```
15_Calculate_block_level_partition_based_of_FC.py
```

### 11. Analyzing the gene block level expression correlation

&emsp;&emsp; To assess whether genes within syntenic gene blocks are co-expressed, we calculated the expression correlation among genes retained in these blocks using below code.

```
16_Get_block_genes_correlations.py
```