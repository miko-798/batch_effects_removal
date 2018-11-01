
load("~/savedFile/batch_sim.RData")

# create an output file

#cat("k", "sigma", "NMI_by_cell_type", "NMI_by_batch", "\n",
#    file="nmi_of_mnn_runs.txt", sep="\t", append=TRUE)

# function call to run MNN
mnn <- function(k, sigma){

mnn_corrected <- mnnCorrect(t(sample1),t(sample2), k=k, sigma=sigma)
mnn_matrix <- cbind(mnn_corrected$corrected[[1]], mnn_corrected$corrected[[2]])
set.seed(0)
tsne_mnn <- Rtsne(t(mnn_matrix), verbose=TRUE, max_iter = 1000)

# color by batch
plot(tsne_mnn$Y, xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=batch.id,
     main=paste0("MNN, k=", k, ", sigma=", sigma))

# color by cell type
plot(tsne_mnn$Y, xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=clust.cols,
     main=paste0("MNN, k=", k, ", sigma=", sigma))

## GET K means clusters ##
# getting the two dimension matrix
df_tsne_mnn <- as.data.frame(tsne_mnn$Y)
# Creating k-means clustering model, and assigning the result to the data used to create the tsne
mnn_kmeans=kmeans(df_tsne_mnn, 3)

# Centroid Plot against 1st 2 discriminant functions
#library(fpc)
#plotcluster(df_tsne_mnn, mnn_kmeans$cluster, main=paste0("MNN clusters, k=", k, ", sigma=", sigma))

# Draw tsne showing the kmeans clustering results
plot(df_tsne_mnn, main=paste0("MNN clusters, k=", k, ", sigma=", sigma),
     xlab = 't-SNE 1', ylab = 't-SNE 2', 
     pch=c(16, 17)[batch.id], col=cg_kmeans$cluster)


# NMI(c(mnn_kmeans$cluster), cell_type)   
# NMI(c(mnn_kmeans$cluster), c(batch.id))

out1 <- capture.output(NMI(c(mnn_kmeans$cluster), cell_type))
out2 <- capture.output(NMI(c(mnn_kmeans$cluster), c(batch.id)))

#cat(k, sigma, out1, out2, "\n", file="nmi_of_mnn_runs.txt", 
#    sep="\t", append=TRUE)

}

# Run MNN with different parameters

#k_list <- c(20)
k_list <- c(5, 10, 40, 50, 70, 80, 100, 120)
sigma_list <- c(0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 1.0)

for (i in 1:length(k_list)) {
  for (j in 1:length(sigma_list)){
    mnn(k_list[i], sigma_list[j])
  }
}

# good results

# 80	0.1	[1] 0.691821	[1] 0.2861519	
### 40	1	  [1] 0.5900187	[1] 0.2655852	
# 20	0.5	[1] 0.5742215	[1] 0.2600867
# 10	0.3	[1] 0.6405602	[1] 0.2852963	
# 10	0.1	[1] 0.5631225	[1] 0.2561514

mnn(40, 1)
# mnn(20, 0.5) # main parameter used
