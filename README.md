# Removing batch effects in single-cell RNA-seq data.

Project during the summer of 2018, Dr. Zemin Zhang's lab, Peking University

**Goal:**
Compare the effectiveness of different approaches to remove batch effects in single-cell RNA-seq datasets.

**Quantification:**
1. Normalized Mutual Information: How well each method separates the cell types and mixes the batch
2. Spearman Correlation: To what extent each method distorts the original data

**Methods used:**
 - in R: Mutual Nearest Neighbors (MNN), Seurat, Limma, ComBat, Quantile Normalization
 - in Python: Batch balanced K nearest neighbors (BBKNN) 

**Workflow:**
Real dataset (expression matrix) --> simulated batches --> combined batches
--> removed batch effects --> compared results 

**An example output from BBKNN:**
The first row shows the clustering results before removing batch, while the second row shows results after batch removal.
![Before vs. After batch removal using BBKNN](https://github.com/miko-798/batch_effects_removal/blob/master/bbknn.png)


