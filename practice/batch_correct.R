require("scran")
require("Rtsne")
require("Seurat")

source("/home/miko/simulateBatches.R")
#source("/Users/Miko/Desktop/MNN/runMN.R")
#source("/Users/Miko/Desktop/MNN/fastMNN.R")

Batch1 <- B1
Batch2 <- B2i
# combine the two batches for plotting
Batch <- cbind(Batch1, Batch2)

dim(Batch1)
dim(Batch2)
dim(Batch)

summary(factor(Batch1))
summary(factor(Batch2))

batch.id <- rep(1:2, c(ncol(Batch1), ncol(Batch2)))
clust.cols <- c(clust1, clust2)

# TSNE before correction
tsne_batch <- Rtsne(t(Batch), dims = 2, perplexity=30, verbose=TRUE, max_iter = 50)
plot(tsne_batch$Y, main="Uncorrected", xlab = 't-SNE 1', ylab = 't-SNE 2', col=clust.cols)


# MNN correction
corrected <- mnnCorrect(Batch1,Batch2, k=100) # mnnCorrect will take several minutes to run
#head(corrected$pairs[[2]])

# fastMNN:
#fast_mnn <- fastMNN(Batch1,Batch2)

#total_pairs <- nrow(corrected$pairs[[2]])
#n_unique_batch2 <- length(unique((corrected$pairs[[2]][,1])))
#n_unique_batch1 <- length(unique((corrected$pairs[[2]][,2])))

mnn_correct <- cbind(corrected$corrected[[1]], corrected$corrected[[2]])
tsne_mnn_correct <- Rtsne(t(mnn_correct),  dims = 2, perplexity=30, verbose=TRUE, max_iter = 50)
plot(tsne_mnn_correct$Y, main='MNN', xlab = 't-SNE 1', ylab = 't-SNE 2', col=clust.cols)

# Seurat


# limma
library(limma)
Xlm <- removeBatchEffect(Batch, factor(batch.id))
#all.dists2.lm <- as.matrix(t(Xlm))
all.dists2.lm <- as.matrix(dist(t(Xlm)))

set.seed(0)
tsne.lm <- Rtsne(all.dists2.lm, is_distance=TRUE)#, perplexity = 0.9)
plot(tsne.lm$Y, main='limma', xlab = 't-SNE 1', ylab = 't-SNE 2', col=clust.cols)


#Combat 
library(sva)
cleandat <- ComBat(Batch, factor(batch.id), mod=NULL, prior.plots = FALSE)
all.dists.combat <- as.matrix(dist(t(cleandat)))

set.seed(0)
tsne.combat<-Rtsne(all.dists.combat, is_distance=TRUE)
plot(tsne.combat$Y, main='ComBat', xlab = 't-SNE 1', ylab = 't-SNE 2', col=clust.cols)




