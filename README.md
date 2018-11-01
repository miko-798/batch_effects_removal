# batch_effects_removal
This repo contains my experiments with removing batch effects in single-cell RNA-seq data, using multiple different approaches.

Goal:
Compare the effectiveness of different approaches to remove batch effects in single-cell RNA-seq datasets.

Quantification: 
1. Normalized Mutual Information: How well each method separates the cell types and mix the batch
2. Spearman Correlation: To what extent each method distort the original data

Methods used:
 - in R: MNN, Seurat, Limma, ComBat, Quantile Normalization
 - in Python: BBKNN 

Workflow:
Real dataset (expression matrix) --> simulated batches --> combined batches
--> removed batch effects --> compared results 


