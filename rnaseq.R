

# Script analysis tailored by Vivek Kanpa, adapted from template
# provided by Dr. Luis Iglesias Martinez.


# SET PROJECT DIRECTORY PATH (path to this project folder) HERE
path = "Documents/UCD/Bio_Principles/Assignment_2"
setwd(path)


# Inspect ERBB2-mutant patients
her2 = read.table(file = 'brca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep = '\t', header = TRUE)
her2_patients = her2$Patient.ID

# Inspect Clinical File.
data_patient  = read.delim("data_clinical_patient.txt")

# Inspect RNASeq File.
data_Rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# Inspect CNA File, subset to ERBB2 CNA data, and reformat to clean ERBB2_cna
data_cna = read.delim("data_cna.txt")
ERBB2_cna = data_cna[data_cna$Hugo_Symbol=="ERBB2",]
rownames(ERBB2_cna) <- ERBB2_cna$Hugo_Symbol
ERBB2_cna = subset(ERBB2_cna, select = -c(Hugo_Symbol, Entrez_Gene_Id) )

# What are the copy numbers associated with ERBB2 for each patient?
erbb2_indx = which(data_cna[,1] == 'ERBB2')
hist(as.numeric(data_cna[erbb2_indx,-c(1,2)]))

# MAKING SURE THE DIMENSIONS ARE EQUAL BETWEEN ASSAY AND METADATA
assay = as.matrix(data_Rnaseq[,-c(1,2)])
common_columns <- intersect(colnames(assay), colnames(ERBB2_cna))
assay_filtered <- assay[, common_columns]
ERBB2_cna_filtered <- ERBB2_cna[, common_columns]

# transpose ERBB2_cna_filtered, preserve row and column names
library(data.table)
ERBB2_cna_t <- transpose(ERBB2_cna_filtered)
rownames(ERBB2_cna_t) <- colnames(ERBB2_cna_filtered)
colnames(ERBB2_cna_t) <- rownames(ERBB2_cna_filtered)

# store the erbb2 copy number values as a list and build metadata.
metadata <- ERBB2_cna_t
erbb2 <- metadata$ERBB2

# Replace copy number alteration values less than or equal to 0 with 0, 
# and values greater than 1 with 1. Refill the metadata with these binary-ized values.
erbb2[erbb2 < 0] <- 0
erbb2[erbb2 > 1] <- 1
metadata$ERBB2 <- erbb2
row.names(metadata) <- 1:nrow(metadata)

cna_num = which(colnames(ERBB2_cna_t) =="ERBB2")

metadata[is.na(metadata)] =0
colnames(metadata) = "Amplified"

# Install BiocManager if not yet installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install DeSeq2 if not yet installed
BiocManager::install("DESeq2")

library(DESeq2)

# Build DESeq Object
assay_filtered[is.na(assay_filtered)] = 0  # Impute with zeros the NA
assay_filtered[assay_filtered<0] = 0

# build a dds object with metadata, assay data, and a design feature 
dds <- DESeqDataSetFromMatrix(countData = round(assay_filtered),
                              colData = metadata,
                              design = ~ Amplified)

# Filter data
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Normalize data
dds <- DESeq(dds)

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
pc = prcomp(assay(rld))

# Plot 
plot(pc$rotation[,1], pc$rotation[,2], col = 1+(metadata$Amplified), pch = 19,
     main = "PCA plot of amplified gene copy numbers for ERBB2+ BC patients")
legend("topright", legend = unique(metadata$Amplified), col = 1+unique(metadata$Amplified), pch = 1)

# Variance explained by each principal component
var_explained <- pc$sdev^2
total_var <- sum(var_explained)

# Calculate percent variance explained by each PC
percent_var_explained <- var_explained / total_var * 100

# Display percent variance explained by each PC
percent_var_explained

# Get Results
res <- results(dds)

# Summary
summary(res)
rownames(res) = data_Rnaseq[keep,1]

# Significantly Differentially Expressed
signif = which(res$padj<0.05)
deg = res[signif,]

# Install EnhancedVolcano if not installed yet
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

# Build a volcano plot of p value across lFC for genes differentially expressed
# in amplified/not amplified ERBB2 samples
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

# Separate the significantly differentially expressed genes
dup = deg[deg[,2]>0.,]

ddown = deg[deg[,2]<0.,]

# For Pathway Enrichment we need Entrez IDs
entrez_ids = data_Rnaseq[keep,2]

entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]

# Pathway Enrichment
BiocManager::install("clusterProfiler")

library(clusterProfiler)

# Do a KEGG pathway over-representation analysis (p < 0.05)
all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

# plot the KEGG pathways by rank of over-representation
dotplot(all_paths, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)


# plot the KEGG pathways of OVEREXPRESSED GENES (p < 0.1; none signif. @ p < 0.05)
up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.1)
head(up_paths)
dotplot(up_paths, 
        showCategory = 10, 
        title = "Overexpressed Genes Pathways",
        font.size = 8)


# plot the KEGG pathways of UNDEREXPRESSED GENES (p < 0.1; none signif. @ p < 0.05)
down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.1)
head(down_paths)
dotplot(down_paths, 
        showCategory = 10, 
        title = "Underexpressed Genes Pathways",
        font.size = 8)
