# Differential Expression Analysis of Breast Cancer PBMCs (GSE174431)

## 📌 Project Overview
This project demonstrates a **bioinformatics pipeline** for analyzing RNA-seq data using **R and DESeq2**.  
The dataset comes from **GSE174431** (GEO), which contains RNA-seq profiles of **lineage-positive (Lin+)** and **lineage-negative (Lin–)** populations isolated from PBMCs of metastatic breast cancer patients.  

The objective was to identify **differentially expressed genes (DEGs)** and visualize sample clustering and expression patterns.

---

## 🔬 Methods
1. Downloaded raw exon count files from **GEO** (`GSE174431_RAW.tar`).  
2. Processed and merged count files into a matrix in **R**.  
3. Annotated samples using metadata (`LinPositive vs LinNegative`).  
4. Performed **DESeq2 analysis**:
   - Normalization
   - PCA for QC
   - Differential expression testing
   - Multiple testing correction (FDR < 0.05)  
5. Generated visualizations:
   - PCA plot
   - Volcano plot
   - Heatmap of top 20 DEGs  

---

## 📊 Results

### PCA Plot
Samples cluster clearly by lineage group (Lin+ vs Lin–), confirming biological separation.  
![PCA Plot](Rplot05.png)

### Volcano Plot
Differentially expressed genes between **Lin+** and **Lin–**.  
![Volcano Plot](Rplot04.png)

### Heatmap of Top 20 DEGs
The top 20 genes show distinct expression patterns between groups.  
![Heatmap](Rplot03.png)

### Top 20 DE Genes (Table)
| Gene | log2FC | Adj.P.Val |
|------|--------|-----------|
| ExampleGene1 | +2.45 | 3.2e-06 |
| ExampleGene2 | -1.78 | 7.1e-05 |
| ... | ... | ... |

---

## ✅ Conclusion
- **Lin– cells** (lineage-negative) display expression consistent with **stem-like / progenitor states**.  
- **Lin+ cells** (lineage-positive) show **differentiated gene expression profiles**.  
- Results align with the hypothesis that Lin– cells may act as **cancer stem-like populations** in breast cancer.  

---

## 🛠 Tools & Packages
- **R / Bioconductor**  
  - `DESeq2`, `GEOquery`, `pheatmap`, `ggplot2`  
- **Data Source:** GEO (GSE174431)  

---

## 📂 Files in Repository
- `results.csv` → full DESeq2 results  
- `Rplot01.png` → Heatmap of top 20 DEGs  
- `Rplot02.png` → Volcano plot  
- `Rplot03.png` → PCA plot  
- `GSE174431_DESeq2_Report_with_Top20.pdf` → Final report (publication-style)

---

## ✨ Portfolio Value
This project demonstrates:
- RNA-seq preprocessing
- DESeq2 workflow
- Data visualization (PCA, Volcano, Heatmap)
- Clear scientific reporting  

It serves as a **sample analysis report** for freelance bioinformatics or academic consulting.

---