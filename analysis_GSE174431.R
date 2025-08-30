#if(!requireNamespace("BiocManager", quietly=TRUE))
      #install.packages("BiocManager")
      #BiocManager::install(c("GEOquery", "DESeq2", "pheatmap"))
if(!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")


if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(pheatmap)


      
# Fetch GEO dataset (expression + metadata)
gse <- getGEO("GSE174431", GSEMatrix = TRUE)

# Extract expression matrix + phenotype data
exprSet <- exprs(gse[[1]])   # expression data
pheno   <- pData(gse[[1]])   # metadata (conditions, etc.)

dim(exprSet)
head(pheno[,1:10])

getGEOSuppFiles("GSE174431")
# Point to your RAW folder
path <- "/Users/haniehroodashty/GSE174431/GSE174431_RAW"
# List all counts files (with or without .gz)
files <- list.files(path, pattern = "counts", full.names = TRUE)
files

library(dplyr)

read_count_file <- function(f) {
  df <- read.delim(gzfile(f), header = TRUE)
  df <- df %>% 
    select(feature_id, counts) %>% 
    group_by(feature_id) %>%        # group by feature_id
    summarise(counts = sum(counts)) # sum counts if duplicates exist
  colnames(df)[2] <- basename(f)
  return(df)
}
# List all files ending with counts.txt.gz
files <- list.files(path, pattern = "counts.txt.gz", full.names = TRUE)
# Read each one
count_list <- lapply(files, read_count_file)
# Merge by feature_id
counts <- Reduce(function(x, y) full_join(x, y, by = "feature_id"), count_list)

# Convert tibble to data.frame
counts <- as.data.frame(counts)


# Row names = feature_id
rownames(counts) <- counts$feature_id
counts <- counts[ , -1]   # drop feature_id column

#Clean Column Names (Optional)
colnames(counts) <- gsub("_rna.*", "", colnames(counts))  # keep only GSM IDs
head(counts[ , 1:5])

#If you want gene-level counts, you’ll need to summarize exons → genes (usually by collapsing exon IDs with the same gene name).
counts$gene <- sub("_exon.*", "", rownames(counts))  # strip exon info
counts_gene <- counts %>% 
  group_by(gene) %>% 
  summarise(across(everything(), sum))
#Now you have a gene-level count matrix, ready for DESeq2.

#Get phenotype info (conditions, labels, etc.):
# See available columns
colnames(pheno)
# Make sure rownames of pheno = sample IDs
rownames(pheno) <- pheno$geo_accession

# Keep only samples present in counts
pheno <- pheno[colnames(counts), ]

colnames(pheno)   # to see possible columns like "title", "source_name_ch1", etc.
head(pheno$title) # many studies store condition here

#Create a clean group column
#This will give you the counts of samples in each group
pheno$group <- ifelse(grepl("Lin\\+", pheno$title), "LinPositive", "LinNegative")
table(pheno$group)

#counts_gene and pheno should match perfectly
#Remove the NA row:
pheno <- pheno[!is.na(rownames(pheno)), ]


# Convert tibble to data.frame first
counts_gene <- as.data.frame(counts_gene)

# Move gene IDs into rownames
rownames(counts_gene) <- counts_gene$gene
counts_gene$gene <- NULL   # drop the "gene" column

# Make sure it's numeric
counts_gene <- as.matrix(counts_gene)
#Now you have a matrix with genes as rownames and GSM IDs as column names.

pheno <- pheno[colnames(counts_gene), ]
pheno$group <- factor(pheno$group)

dds <- DESeqDataSetFromMatrix(
  countData = counts_gene,
  colData   = pheno,
  design    = ~ group
)
dds <- DESeq(dds)
res <- results(dds)
# Order genes by significance
resOrdered <- res[order(res$padj), ]
# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
# Pick top 20 DEGs
topGenes <- head(rownames(resOrdered), 20)
mat <- assay(vsd)[topGenes, ]
# Heatmap
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=TRUE, show_colnames=TRUE,
         main="Top 20 Differentially Expressed Genes")
#Volcano Plot
res$threshold <- as.factor(abs(res$log2FoldChange) > 1 & res$padj < 0.05)

ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot: Metastatic vs Healthy PBMCs",
       x="Log2 Fold Change",
       y="-Log10 Adjusted p-value")

#PCA Plot
pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal() +
  labs(title="PCA Plot: GSE174431 Breast Cancer PBMCs")
