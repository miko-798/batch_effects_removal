gc()
rm(list = ls())
sce <- readRDS("./analysis_results/10X/clustering/sce.qc.CD45P.analyzed.rds")
# Creat Seurat
library(Seurat)
cellnames.bytissue <- split(colnames(sce), sce$tissue)
normal.sr <- CreateSeuratObject(counts(sce)[, cellnames.bytissue$Normal], project = "Normal" )
tumor.sr <- CreateSeuratObject(counts(sce)[, cellnames.bytissue$Tumor], project = "Tumor" )
blood.sr <-  CreateSeuratObject(counts(sce)[, cellnames.bytissue$Blood], project = "Blood" )
normal.sr <- AddMetaData(normal.sr, as.data.frame(colData(sce)[cellnames.bytissue$Normal, c("donor", "tissue", "Detection")]))
tumor.sr <- AddMetaData(tumor.sr, as.data.frame(colData(sce)[cellnames.bytissue$Tumor, c("donor", "tissue", "Detection")]))
blood.sr <- AddMetaData(blood.sr, as.data.frame(colData(sce)[cellnames.bytissue$Blood, c("donor", "tissue", "Detection")]))

rm(sce); gc();


mito.genes <- grep(pattern = "^MT-", x = rownames(x = normal.sr@data), value = TRUE)
percent.mito <- Matrix::colSums(normal.sr@raw.data[mito.genes, ])/Matrix::colSums(normal.sr@raw.data)
normal.sr <- AddMetaData(object = normal.sr, metadata = percent.mito, col.name = "percent.mito")


mito.genes <- grep(pattern = "^MT-", x = rownames(x = tumor.sr@data), value = TRUE)
percent.mito <- Matrix::colSums(tumor.sr@raw.data[mito.genes, ])/Matrix::colSums(tumor.sr@raw.data)
tumor.sr <- AddMetaData(object = tumor.sr, metadata = percent.mito, col.name = "percent.mito")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = blood.sr@data), value = TRUE)
percent.mito <- Matrix::colSums(blood.sr@raw.data[mito.genes, ])/Matrix::colSums(blood.sr@raw.data)
blood.sr <- AddMetaData(object = blood.sr, metadata = percent.mito, col.name = "percent.mito")


# Normalization
normal.sr@meta.data$stim <- 'Normal'
tumor.sr@meta.data$stim <- 'Tumor'
blood.sr@meta.data$stim <- 'Blood'

normal.sr <- FilterCells(normal.sr, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
normal.sr <- NormalizeData(normal.sr)
normal.sr <- ScaleData(normal.sr, vars.to.regress = c("donor", "nUMI", "percent.mito"))

tumor.sr <- FilterCells(tumor.sr, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
tumor.sr <- NormalizeData(tumor.sr)
tumor.sr <- ScaleData(tumor.sr, vars.to.regress = c("donor", "nUMI", "percent.mito"))

blood.sr <- FilterCells(blood.sr, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
blood.sr <- NormalizeData(blood.sr)
blood.sr <- ScaleData(blood.sr, vars.to.regress = c("donor", "nUMI", "percent.mito"))


# Gene selection for input to CCA
normal.sr <- FindVariableGenes(normal.sr, do.plot = F)
tumor.sr <- FindVariableGenes(tumor.sr, do.plot = F)
blood.sr <- FindVariableGenes(blood.sr, do.plot = F)


g.1 <- rownames(normal.sr@hvg.info[1:1500,])
g.2 <- rownames(tumor.sr@hvg.info[1:1500,])
g.3 <- rownames(blood.sr@hvg.info[1:1500,])

genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(normal.sr@scale.data))
genes.use <- intersect(genes.use, rownames(tumor.sr@scale.data))

# Perform a canonical correlation analysis (CCA)
immune.combined <- RunCCA(normal.sr, tumor.sr, genes.use = genes.use, num.cc = 30)
MetageneBicorPlot(immune.combined, grouping.var = "tissue", dims.eval = 1:25, 
                  display.progress = FALSE)
DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 19:27, do.balanced = TRUE)
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", 
                                 grouping.var = "tissue", 
                                 dims.align = 1:25)

immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", 
                           dims.use = 1:25, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = .6, dims.use = 1:25)


# Visualization

immune.combined <- RunUMAP(immune.combined, reduction.use = "cca.aligned", 
                           dims.use = 1:25, min_dist = 1)
markers <- FindAllMarkers(immune.combined, only.pos = T, print.bar = T, min.pct = 0., logfc.threshold = .5)
library(dplyr)
topmarkers <- markers %>% 
  group_by(cluster) %>% 
  top_n(5, avg_logFC)

saveRDS(immune.combined, "./analysis_results/normal_tumor.CCA.rds")
saveRDS(immune.combined, "./analysis_results/normal_tumor.markers.rds")


# --- Normal and Blood

g.1 <- rownames(normal.sr@hvg.info[1:1500,])
g.3 <- rownames(blood.sr@hvg.info[1:1500,])

genes.use <- unique(c(g.1, g.3))
genes.use <- intersect(genes.use, rownames(normal.sr@scale.data))
genes.use <- intersect(genes.use, rownames(blood.sr@scale.data))

# Perform a canonical correlation analysis (CCA)
immune.combined <- RunCCA(normal.sr, blood.sr, genes.use = genes.use, num.cc = 30)
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", 
                                 grouping.var = "tissue", 
                                 dims.align = 1:25)

immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", 
                           dims.use = 1:25, 
                           do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned", 
                                resolution = .6, dims.use = 1:25)


# Visualization
library(reticulate)
use_python("/home/heyao/tools/basic/anaconda3/bin/python")
umap <- import("umap")
immune.combined <- RunUMAP(immune.combined, reduction.use = "cca.aligned", 
                           dims.use = 1:25, min_dist = 1)
markers <- FindAllMarkers(immune.combined, only.pos = T, print.bar = T, min.pct = 0., logfc.threshold = .5)
library(dplyr)
topmarkers <- markers %>% 
  group_by(cluster) %>% 
  top_n(5, avg_logFC)
DimPlot(immune.combined, reduction.use = "umap", group.by = "tissue")
saveRDS(immune.combined, "./analysis_results/normal_blood.CCA.rds")
saveRDS(immune.combined, "./analysis_results/normal_blood.markers.rds")
# table(immune.combined@meta.data$tissue, immune.combined@meta.data$res.0.6)
FeaturePlot(object = immune.combined, 
            features.plot = c("CD3D", "MKI67", "SELL", "GPR183", "CD8A", "GNLY", "CD79A", "FCGR3A", "MS4A1", "CST3"), 
            min.cutoff = "q9", reduction.use = "umap",
            cols.use = c("lightgrey", "blue"), pt.size = 0.5)
# FeaturePlot(object = immune.combined, 
#             features.plot = c("CD3D", "GNLY", "CD79A", "FCGR3A"), 
#             min.cutoff = "q9", reduction.use = "tsne",
#             cols.use = c("lightgrey", "blue"), pt.size = 0.5)
# immune.combined <- Seurat::SetAllIdent(immune.combined, id = "res0.6")

immune.combined <- readRDS("./analysis_results/normal_tumor.CCA.rds")
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "tissue")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

table(immune.combined@meta.data$res.0.6, immune.combined@meta.data$tissue)
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", 
                     "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", 
                     "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", 
                     "HBB", "TSPAN13", "IL3RA")
sdp <- SplitDotPlotGG(immune.combined, 
                      genes.plot = rev(markers.to.plot), 
                      cols.use = c("blue", "red"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, 
                      grouping.var = "stim")



DoHeatmap(object = immune.combined, genes.use = topmarkers$gene, 
          slim.col.label = TRUE, remove.key = TRUE)
saveRDS(immune.combined, )
