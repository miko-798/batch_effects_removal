# continued from project_batch.R

library(Seurat)
load("~/savedFile/second_dataset.RData")

s3_ncells <- 295 # only cell type 1
s4_ncells <- 250 # only cell type 3
total_ncells <- s3_ncells + s4_ncells

# transpose for Seurat purposes
sample1 <- t(sample3)
sample2 <- t(sample4)

# Set up seurat objects
s1_seu <- CreateSeuratObject(raw.data = sample1, project = "Batch3")
s1_seu@meta.data$stim <- "SAMPLE1"
#s1_seu <- NormalizeData(s1_seu)
#s1_seu <- ScaleData(s1_seu, display.progress = F)

s2_seu <- CreateSeuratObject(raw.data = sample2, project = "Batch4")
s2_seu@meta.data$stim <- "SAMPLE2"
#s2_seu <- NormalizeData(s2_seu)
#s2_seu <- ScaleData(s2_seu, display.progress = F)

# Gene selection for input to CCA
"""
s1 <- FindVariableGenes(
  object = s1_seu,
  do.plot = F
)
s2 <- FindVariableGenes(
  object = s2_seu,
  do.plot = F
)

g.1 <- head(rownames(s1_seu@hvg.info), 1000)
g.2 <- head(rownames(s2_seu@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2)) # 1585
genes.use <- intersect(genes.use, rownames(s1_seu@data))
genes.use <- intersect(genes.use, rownames(s2_seu@data))  # 1585

"""
# Perform CCA (identify common sources of variation) # use all genes
batch.combined <- RunCCA(s1_seu, s2_seu, genes.use = rownames(s1_seu@data),
                         num.cc = 15, 
                         scale.data = FALSE)

# Assign stim values for grouping
#batch.combined@meta.data$stim <- "SAMPLE1"
#batch.combined@meta.data[(s1_ncells+1):total_ncells,ncol(batch.combined@meta.data)] <- "SAMPLE2"

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = batch.combined, reduction.use = "cca", group.by = "stim",
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = batch.combined, features.plot = "CC1", group.by = "stim",
              do.return = TRUE)
plot_grid(p1, p2)
PrintDim(object = batch.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

p3 <- MetageneBicorPlot(batch.combined, grouping.var = "stim", dims.eval = 1:15, 
                        display.progress = FALSE)
DimHeatmap(object = batch.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)

# Align the CCA subspaces
genes.use = rownames(s1_seu@data)
batch.combined <- AlignSubspace(batch.combined, reduction.type = "cca", grouping.var = "stim", 
                            dims.align = 1:15)
# stored here
#batch.combined@dr$cca.aligned

#visualize
p1 <- VlnPlot(object = batch.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = batch.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

# integrated analysis
# t-SNE and Clustering
batch.combined <- RunTSNE(batch.combined, reduction.use = "cca.aligned", dims.use = 1:15, 
                           do.fast = T)
batch.combined <- FindClusters(batch.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:15)

# Visualization
p1 <- TSNEPlot(batch.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(batch.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

# NMI
"""
seurat_matrix <- as.matrix(batch.combined@data)
seurat_dist_matrix <- as.matrix(dist(t(seurat_matrix)))

set.seed(0)
tsne_seu <- Rtsne(seurat_dist_matrix, is_distance=TRUE)
plot(tsne_seu$Y, main='seurat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
df_tsne_seu <- as.data.frame(tsne_seu$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
seu_kmeans=kmeans(df_tsne_seu, 3)

NMI(c(seu_kmeans$cluster), cell_type)  # 0.5901593
NMI(c(seu_kmeans$cluster), c(batch.id))  # 0.6559789

set.seed(0)
tsne_embed <- Rtsne(cell_embed, perplexity=15)
plot(tsne_embed$Y, main='seurat_embed', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

## GET K means clusters ##
df_tsne_embed <- as.data.frame(tsne_embed$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
embed_kmeans=kmeans(df_tsne_embed, 3)

NMI(c(embed_kmeans$cluster), cell_type)  # 0.3219745
NMI(c(embed_kmeans$cluster), c(batch.id)) # 0.1328245
"""

###################################

# REAL data after Seurat run: cell_gene

cell_embed <- batch.combined@dr$cca.aligned@cell.embeddings
gene_load <- t(batch.combined@dr$cca@gene.loadings)
cell_gene <- cell_embed %*% gene_load

#save(file="~/savedFile/cell_gene.RData", cell_gene)

#load("~/savedFile/cell_gene.RData")

set.seed(0)
tsne_cg <- Rtsne(cell_gene, perplexity=15)

# plot colored by batch
plot(tsne_cg$Y, main='Seurat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# plot colored by cell type
plot(tsne_cg$Y, main='Seurat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
df_tsne_cg <- as.data.frame(tsne_cg$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
cg_kmeans=kmeans(df_tsne_cg, 2)

# Draw tsne showing the kmeans clustering results
plot(df_tsne_cg, main='Seurat clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=cg_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(cg_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(cg_kmeans$cluster), c(batch.id)))

cat("6. Seurat", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)


###############################


# correlation for each cell
rho_stat <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=cell_gene[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_stat <- c(rho_stat, curr_stat)  # append
}
#rho_stat

dotchart(rho_stat,cex=.7,
         main="Spearman correlation rho: Seurat",
         xlab = "rho",
         color = "black")

#save(file="~/savedFile/batch_seurat.RData", full_sample, seurat_matrix, batch.id)


