# RNAseq-ERBB2-breast-cancer
DGE analysis of ERBB2-positive breast cancer: A mechanistic perspective on cancer-promoting phenotypes RNASeq data from cBioPortal

--------------------------
1. Formatting data
--------------------------
The script takes two data files:
  A. data_Rnaseq  : Contains the number of transcript reads across the genome of each patient with breast invasive carcinoma.
  B. data_cna     : Contains the copy number alteration (CNA) number for each gene of each patient with breast invasive carcinoma.

The first part of the script restructures and merges the data from these two tables to create:
  A. metadata       : The ERBB2 copy number alterations for each patient, normalized to a binary 0-1 indexing.
  B. assay_filtered : The RNA Sequencing read counts for each gene in the genome of each patient with breast invasive carcinoma.

--------------------------
2. DESeq2 analysis
--------------------------
DESeq2 performs a few key operations:
  A. Filtering: filters genes to those with greater than 10 counts in at least 3 groups.
  B. Normalization: DESeq2 estimates size factors for read count normalization, fits a negative binomial generalized linear model to the data, and computes statistics to assess differential expression between ERBB2 amplification (1) or non-amplification (0) corresponding with the metadata.
  C. Variance-Stabilized Transformation: transforms the normalized values for downstream clustering and dimensionality reduction (via PCA) with stabilized variance.

--------------------------
3. Visualization
--------------------------
Visualizations produced in this script include...
  A. a volcano plot, graphing p-value significance over log-Fold change for each gene, communicating patterns and clusters of significantly differentially expressed genes in both ERBB2+/ERBB- groups.
  B. a PCA plot, enabling investigation of clusters and trends between ERBB2+/ERBB2- differential gene expression in breast cancer. Included in this visual analysis is a deep dive into the genes contributing highest to PC1 variance.

--------------------------
4. KEGG analysis
--------------------------
KEGG enriched pathway analysis was performed on...
  A. All significantly differentially expressed genes (regardless of up/down regulation)
  B. Upregulated gene pathway sets
  C. Downregulated gene pathway sets

Dotplot visualizations were produced to assess the scope/size of how a pathway was affected by ERBB2+ amplification.
