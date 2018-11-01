library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)


batch1_seu <- CreateSeuratObject(Batch1)
batch2_seu <- CreateSeuratObject(Batch2)

# violin plot to see the distribution
VlnPlot(
  object = batch1_seu, 
  features.plot = c("nGene", "nUMI"), 
  nCol = 2
)
GenePlot(
  object = batch1_seu, 
  gene1 = "nUMI", 
  gene2 = "nGene"
)

VlnPlot(
  object = batch2_seu, 
  features.plot = c("nGene", "nUMI"), 
  nCol = 2
)
GenePlot(
  object = batch2_seu, 
  gene1 = "nUMI", 
  gene2 = "nGene"
)

# filter cells
batch1_seu <- FilterCells(
  object = batch1_seu, 
  subset.names = c("nUMI"), 
  high.thresholds = c(2e7)
)
batch2_seu <- FilterCells(
  object = batch2_seu, 
  subset.names = c("nUMI"), 
  high.thresholds = c(2e7)
)

# Normalize the data
batch1_seu <- NormalizeData(
  object = batch1_seu, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)
batch2_seu <- NormalizeData(
  object = batch2_seu, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

# highly variable genes
batch1_seu <- FindVariableGenes(
  object = batch1_seu,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)
batch2_seu <- FindVariableGenes(
  object = batch2_seu,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

# scale data for dimensionality reduction
batch1_seu <- ScaleData(
  object = batch1_seu, 
  vars.to.regress = c("nUMI")
)
batch2_seu <- ScaleData(
  object = batch2_seu, 
  vars.to.regress = c("nUMI")
)

batch1_seu <- RunPCA(
  object = batch1_seu, 
  pc.genes = batch1_seu@var.genes, 
  do.print = TRUE, 
  pcs.print = 1:5, 
  genes.print = 5
)

PrintPCA(object = batch1_seu, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = batch1_seu, pcs.use = 1:2)
PCAPlot(object = batch1_seu, dim.1 = 1, dim.2 = 2)
PCHeatmap(
  object = batch1_seu, 
  pc.use = 1:6, 
  cells.use = 500, 
  do.balanced = TRUE, 
  label.columns = FALSE,
  use.full = FALSE
)

# significant PCs
batch1_seu <- JackStraw(
  object = batch1_seu, 
  num.replicate = 100, 
  do.print = FALSE
)

JackStrawPlot(object = batch1_seu, PCs = 1:9)

# clustering cells
batch1_seu <- FindClusters(
  object = batch1_seu, 
  reduction.type = "pca", 
  dims.use = 1:8, 
  resolution = 1.0, 
  print.output = 0, 
  save.SNN = TRUE
)
PrintFindClustersParams(object = batch1_seu)



# visualize clustering results
batch1_seu <- RunTSNE(
  object = batch1_seu,
  dims.use = 1:8,
  do.fast = TRUE
)
TSNEPlot(object = batch1_seu)


