library(Rtsne)
require("scran")
library("R.utils")
require(aricode) 

dataset_Kolod <- loadToEnv("~/input/dataset_Kolod.RData")
matrix1 <- dataset_Kolod$in_X

## START OF CLUSTERING ##
{
matrix1_seu <- CreateSeuratObject(matrix1)

# highly variable genes
matrix1_seu <- FindVariableGenes(
  object = matrix1_seu,
  mean.function = ExpMean, 
  dispersion.function = LogVMR, 
  x.low.cutoff = 0.0125, 
  x.high.cutoff = 3, 
  y.cutoff = 0.5
)

length(x = matrix1_seu@var.genes) # 211

matrix1_seu@scale.data <- matrix1_seu@data
matrix1_seu <- RunPCA(object = matrix1_seu, 
                      pc.genes = matrix1_seu@var.genes, 
                      do.print = TRUE,
                      pcs.print = 1:5, 
                      genes.print = 5)

VizPCA(object = matrix1_seu, pcs.use = 1:5)

matrix1_seu <- JackStraw(
  object = matrix1_seu, 
  num.replicate = 100,
  num.pc = 10,
  prop.freq = 0.8
)

JackStrawPlot(object = matrix1_seu, PCs = 1:5)

## PC: 1-3
# clustering of cells
matrix1_seu <- FindClusters(
  object = matrix1_seu, 
  reduction.type = "pca", 
  dims.use = 1:3, 
  resolution = 0.2, 
  print.output = 0, 
  save.SNN = TRUE
)
PrintFindClustersParams(object = matrix1_seu)

# run tSNE to visualize clustering results
matrix1_seu <- RunTSNE(
  object = matrix1_seu,
  dims.use = 1:3,
  do.fast = TRUE
)
TSNEPlot(object = matrix1_seu)

saveRDS(matrix1_seu, file = "~/savedRDS/matrix1_seu.rds")

}
## END OF CLUSTERING ##

#####################################

## START OF BATCH SIMULATION ##
{
# get the indices for cell types 1-3
#which(dataset_Kolod$true_labs == 3)

# two batches, cell by gene
s1_ncells <- 454
s2_ncells <- 409
ngenes <- 5000
sample1 <- matrix1[1:454,]    # cell type 1 and 2
sample2 <- matrix1[296:704,]  # cell type 2 and 3

sample1 <- data.matrix(sample1)
sample2 <- data.matrix(sample2)

# create a palette to show diff. colors for diff. cell types
palette <- c("pink", "blue", "gold")
# clust1: the colors for each cell by its cell type, for sample1
clust1 <- palette[dataset_Kolod$true_labs[1:454]]
# clust2: the colors for each cell by its cell type, for sample2
clust2 <- palette[dataset_Kolod$true_labs[296:704]]
# integers indicating the cell type (1,2 or 3) in the entire sample
cell_type <- c(dataset_Kolod$true_labs[1:454], dataset_Kolod$true_labs[296:704])


##  CREATE BATCH EFFECTS  ##

# Add batch effects and random noise.
sample1 <- sample1 + matrix(rep(rnorm(ngenes), each=s1_ncells), ncol=ngenes) # gene-specific batch effect (genes are columns)
sample1 <- sample1 + rnorm(ngenes*s1_ncells)#, mean = 0, sd = 0.1) # noise
sample1 <- sample1 + 7 # make the values positive

sample2 <- sample2 + matrix(rep(rnorm(ngenes), each=s2_ncells), ncol=ngenes) # gene-specific batch effect (genes are columns)
sample2 <- sample2 + rnorm(ngenes*s2_ncells)#, mean = 0, sd = 0.1) # noise
sample2 <- sample2 + 7

# No log transform to avoid negative numbers

# plot simulated batch 1 and 2
plot(sample1, main="Batch1", pch=16, cex=1.5, col=clust1)
plot(sample2, main="Batch2", pch=17, cex=1.5, col=clust2)

# Assign cell and gene names
rownames(sample1) <- paste0("Cell", seq_len(s1_ncells), "-1")
rownames(sample2) <- paste0("Cell", seq_len(s2_ncells), "-2")

colnames(sample1) <- paste0("Gene", seq_len(ngenes))
colnames(sample2) <- paste0("Gene", seq_len(ngenes))
}
save(file="~/savedFile/batch_sim.RData", sample1, sample2, clust1, clust2, cell_type)

## END OF BATCH SIMULATION: sample1, sample2, clust1, clust2 ##


### FOR PYTHON METHOD: BBKNN ###
{
# export batches
write.csv(sample1, "sample1.csv")
system("gzip sample1.csv")
write.csv(sample2, "sample2.csv")
system("gzip sample2.csv")

# export cell type and batch info
write.table(cell_type, "cell.type.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
batch.id <- rep(1:2, c(nrow(sample1), nrow(sample2)))
write.table(batch.id, "batch.id.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
}
##############################


## START OF BATCH COMBINATION #
{
load("~/savedFile/batch_sim.RData")

# combine the two batches, and get info
full_sample <- rbind(sample1, sample2)
clust.cols <- c(clust1, clust2)   # cell type colors
batch.id <- rep(1:2, c(nrow(sample1), nrow(sample2)))

# draw TSNE before correction
set.seed(0)
tsne_unc <- Rtsne(full_sample)

# color by batch
plot(tsne_unc$Y, main='Uncorrected', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_unc$Y, main="Uncorrected", xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
df_tsne_unc <- as.data.frame(tsne_unc$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
unc_kmeans=kmeans(df_tsne_unc, 3)

# Draw tsne showing the kmeans clustering results
plot(df_tsne_unc, main='Original sample clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=unc_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(unc_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(unc_kmeans$cluster), c(batch.id)))

cat("None: Uncorrected", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)

}
## END OF BATCH COMBINATION ##

##################################

## START OF BATCH REMOVAL ##
{
#create an output file
cat("Method", "NMI_by_cell_type", "NMI_by_batch", "\n",
    file="results.txt", sep="\t", append=TRUE)

# Method 1: MNN correction

# k = 20, sigma = 0.5
mnn_corrected <- mnnCorrect(t(sample1),t(sample2), sigma = 0.5) # need the transpose
mnn_matrix <- cbind(mnn_corrected$corrected[[1]], mnn_corrected$corrected[[2]])
set.seed(0)
tsne_mnn <- Rtsne(t(mnn_matrix), verbose=TRUE, max_iter = 1000)

# color by batch
plot(tsne_mnn$Y, main='MNN', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_mnn$Y, main='MNN', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_mnn <- as.data.frame(tsne_mnn$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
mnn_kmeans=kmeans(df_tsne_mnn, 3)

# Centroid Plot against 1st 2 discriminant functions
#library(fpc)
#plotcluster(df_tsne_mnn, mnn_kmeans$cluster, main="mnn clusters")

plot(df_tsne_mnn, main='MNN clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=mnn_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(mnn_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(mnn_kmeans$cluster), c(batch.id)))

cat("1. MNN", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)


# Method 2: limma

library(limma)
lm_matrix <- removeBatchEffect(t(full_sample), factor(batch.id)) # 5000*863
set.seed(0)
tsne_lm <- Rtsne(t(lm_matrix))

# color by batch
plot(tsne_lm$Y, main='limma', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_lm$Y, main='limma', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_lm <- as.data.frame(tsne_lm$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
lm_kmeans=kmeans(df_tsne_lm, 3)

# Centroid Plot against 1st 2 discriminant functions
#library(fpc)
#plotcluster(df_tsne_lm, lm_kmeans$cluster, main="lm clusters")

# Draw tsne showing the kmeans clustering results
plot(df_tsne_lm, main='limma clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=lm_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(lm_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(lm_kmeans$cluster), c(batch.id)))

cat("2. limma", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)


# Method 3: Quantile normalization

require(preprocessCore)
qn_sample1 <- normalize.quantiles(t(sample1),copy=TRUE)
qn_sample2 <- normalize.quantiles(t(sample2),copy=TRUE)
qn_batch <- cbind(qn_sample1, qn_sample2) 
#qn_dist_matrix <- as.matrix(dist(t(qn_batch)))

set.seed(0)
tsne_qn<-Rtsne(t(qn_batch))

# color by batch
plot(tsne_qn$Y, main='Quantile normalization', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_qn$Y, main='Quantile normalization', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_qn <- as.data.frame(tsne_qn$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
qn_kmeans=kmeans(df_tsne_qn, 3)

# Centroid Plot against 1st 2 discriminant functions
#plotcluster(df_tsne_qn, qn_kmeans$cluster, main="qn clusters")

# Draw tsne showing the kmeans clustering results
plot(df_tsne_qn, main='Quantile normalization clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=qn_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(qn_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(qn_kmeans$cluster), c(batch.id)))

cat("3. Quantile normalization", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)

# Method 4: ComBat 

library(sva)
combat_matrix <- ComBat(t(full_sample), factor(batch.id), mod=NULL, prior.plots = FALSE)
set.seed(0)
tsne_combat<-Rtsne(t(combat_matrix))

# color by batch
plot(tsne_combat$Y, main='ComBat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_combat$Y, main='ComBat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_combat <- as.data.frame(tsne_combat$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
combat_kmeans=kmeans(df_tsne_combat, 3)

# Centroid Plot against 1st 2 discriminant functions
#plotcluster(df_tsne_combat, combat_kmeans$cluster, main="combat clusters")

# Draw tsne showing the kmeans clustering results
plot(df_tsne_combat, main='Combat clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=combat_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(combat_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(combat_kmeans$cluster), c(batch.id)))

cat("4. ComBat", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)

# Method 4.1: QN + ComBat 

qn_combat_matrix <- ComBat(qn_batch, factor(batch.id), mod=NULL, prior.plots = FALSE)
set.seed(0)
tsne_qn_combat<-Rtsne(t(qn_combat_matrix))

# color by batch
plot(tsne_qn_combat$Y, main='QN + ComBat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id)

# color by cell type
plot(tsne_qn_combat$Y, main='QN + ComBat', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols)

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_qn_combat <- as.data.frame(tsne_qn_combat$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
qn_combat_kmeans=kmeans(df_tsne_qn_combat, 3)

# Centroid Plot against 1st 2 discriminant functions
#plotcluster(df_tsne_qn_combat, qn_combat_kmeans$cluster, main="qn+combat clusters")

# Draw tsne showing the kmeans clustering results
plot(df_tsne_qn_combat, main='QN+Combat clusters', xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=qn_combat_kmeans$cluster)

# output NMI values
out1 <- capture.output(NMI(c(qn_combat_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(qn_combat_kmeans$cluster), c(batch.id)))

cat("5. QN + ComBat", out1, out2, "\n", file="results.txt", 
    sep="\t", append=TRUE)

# Method 5: Seurat: See project_seurat.R

}
## END OF BATCH REMOVAL ##

#################################
## Calculate Spearman correlation of the expression values ##
# before vs. after batch removal

{
total_ncells = 863

# 1: mnn
rho_mnn <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=t(mnn_matrix)[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_mnn <- c(rho_mnn, curr_stat)  # append
}
#rho_mnn

# visualize all rho
dotchart(rho_mnn,cex=.7,
         main="Spearman correlation rho: MNN",
         xlab = "rho",
         color = "black")

# 2: limma
rho_limma <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=t(lm_matrix)[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_limma <- c(rho_limma, curr_stat)  # append
}
#rho_limma

# visualize all rho
dotchart(rho_limma,cex=.7,
         main="Spearman correlation rho: limma",
         xlab = "rho",
         color = "black")


# 3: combat
rho_combat <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=t(combat_matrix)[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_combat <- c(rho_combat, curr_stat)  # append
}
#rho_combat

# visualize all rho
dotchart(rho_combat,cex=.7,
         main="Spearman correlation rho: ComBat",
         xlab = "rho",
         color = "black")

# 4: QN
rho_qn <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=t(qn_batch)[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_qn <- c(rho_qn, curr_stat)  # append
}
#rho_qn

# visualize all rho
dotchart(rho_qn,cex=.7,
         main="Spearman correlation rho: Quantile normalization",
         xlab = "rho",
         color = "black")

# 5: QN + combat
rho_qn_combat <- c()
for (i in 1:total_ncells) {  # total_ncells
  curr_stat <- cor.test(x=full_sample[i,], y=t(qn_combat_matrix)[i,], method = 'spearman')
  curr_stat <- curr_stat$estimate     # get rho
  rho_qn_combat <- c(rho_qn_combat, curr_stat)  # append
}
#rho_qn_combat

# visualize all rho
dotchart(rho_qn_combat,cex=.7,
         main="Spearman correlation rho: QN + ComBat",
         xlab = "rho",
         color = "black")
}


