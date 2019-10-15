Data sets included:

**geno.csv:** Allele for each genetic markers for each maize line (column 1) in the format [-1,0,1]  corresponding to [aa, Aa, AA] with the more common allele (AA) designated as 1. The genetic marker positions were converted from maize B73 reference genome A Golden Path v2 (AGPv2) to AGPv4.37. The AGPv2 genetic markers that did not map to AGPv4.37 and genetic markers with a minor allele frequency less than 5% were removed, resulting in 332,178 genetic markers.

**transcripto.csv:** RNA-Seq derived transcriptomic data from whole-seedling tissue (i.e. root and shoot) at the V1 stage from each maize line (column 1)(Hirsch et al., 2014). Only transcripts that map to AGPv4.37 and had variation in transcript levels in >5% of maize lines were included (31,238 transcripts). Values are loge + 1 transformed raw transcripts per million count data. 

**pheno.csv:** Height (HT), flowering-time (FT), and grain yield (YLD) values across all maize lines (column 1: ID) from (Hansey et al., 2011)

**CVFs.csv:*** Key for cross-validation fold assignments (using 5-fold cross validation) for each maize line (column 1: ID) across the 100 replicates used in the study. 

**v3_v4_xref.txt:** Key for the mapping of AGPv2/3 to AGPv4.37. Includes the v3 gene model name (v3_gene_model, canonical_transcript), v3 location (v3_start, v3_end, v3_chr), other names (genbank_accession, entrez_gene), v4 name (v4_gene_model), and v4 location (v4_start, v4_end, v4_chr).
