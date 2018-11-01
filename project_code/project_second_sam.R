library(Rtsne)
require("scran")
library("R.utils")
require(aricode)

dataset_Kolod <- loadToEnv("~/input/dataset_Kolod.RData")
matrix1 <- dataset_Kolod$in_X

# Batch simulation
{
s3_ncells <- 295 # only cell type 1
s4_ncells <- 250 # only cell type 3
ngenes <- 5000
sample3 <- matrix1[1:295,]    # cell type 1 
sample4 <- matrix1[455:704,]  # cell type 3

sample3 <- data.matrix(sample3)
sample4 <- data.matrix(sample4)

# create a palette to show diff. colors for diff. cell types
palette_1 <- c("pink", "gold")
color1 <- palette[dataset_Kolod$true_labs[1:295]]
color2 <- palette[dataset_Kolod$true_labs[455:704]]
# integers indicating the cell type in the entire sample
cell_type_1 <- c(dataset_Kolod$true_labs[1:295], dataset_Kolod$true_labs[455:704])


##  CREATE BATCH EFFECTS  ##

# Add batch effects and random noise.
sample3 <- sample3 + matrix(rep(rnorm(ngenes), each=s3_ncells), ncol=ngenes) # gene-specific batch effect (genes are columns)
sample3 <- sample3 + rnorm(ngenes*s3_ncells)#, mean = 0, sd = 0.1) # noise
sample3 <- sample3 + 7 # make the values positive

sample4 <- sample4 + matrix(rep(rnorm(ngenes), each=s4_ncells), ncol=ngenes) # gene-specific batch effect (genes are columns)
sample4 <- sample4 + rnorm(ngenes*s4_ncells)#, mean = 0, sd = 0.1) # noise
sample4 <- sample4 + 7

# No log transform to avoid negative numbers

# plot simulated batch 3 and 4
plot(sample3, main="Batch3", pch=16, cex=1.5, col=color1)
plot(sample4, main="Batch4", pch=17, cex=1.5, col=color2)

# Assign cell and gene names
rownames(sample3) <- paste0("Cell", seq_len(s3_ncells), "-3")
rownames(sample4) <- paste0("Cell", seq_len(s4_ncells), "-4")

colnames(sample3) <- paste0("Gene", seq_len(ngenes))
colnames(sample4) <- paste0("Gene", seq_len(ngenes))
}
save(file="~/savedFile/second_dataset.RData", sample3, sample4, color1, color2, cell_type_1)

### FOR PYTHON METHOD: BBKNN ###
{
# export batches
write.csv(sample3, "sample3.csv")
system("gzip sample3.csv")
write.csv(sample4, "sample4.csv")
system("gzip sample4.csv")

# export cell type and batch info
write.table(cell_type_1, "cell.type_1.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
batch.id.1 <- rep(1:2, c(nrow(sample3), nrow(sample4)))
write.table(batch.id.1, "batch.id.1.txt", sep="\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
}

## START OF BATCH COMBINATION #
{
  load("~/savedFile/second_dataset.RData")
  
  #create an output file
  cat("Method", "NMI_by_cell_type", "NMI_by_batch", "\n",
      file="results.txt", sep="\t", append=TRUE)
  
  # combine the two batches, and get info
  full_sample <- rbind(sample3, sample4)
  clust.cols <- c(color1, color2)   # cell type colors   # pink or gold
  batch.id <- rep(1:2, c(nrow(sample3), nrow(sample4)))  # 1 or 2
  cell_type <- cell_type_1                               # 1 or 3
  
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
  unc_kmeans=kmeans(df_tsne_unc, 2)
  
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
  # Method 1: MNN correction
  
  # k = 20, sigma = 0.5
  mnn_corrected <- mnnCorrect(t(sample3),t(sample4), sigma = 0.5) # need the transpose
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
  mnn_kmeans=kmeans(df_tsne_mnn, 2)
  
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
  lm_kmeans=kmeans(df_tsne_lm, 2)
  
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
  qn_sample3 <- normalize.quantiles(t(sample3),copy=TRUE)
  qn_sample4 <- normalize.quantiles(t(sample4),copy=TRUE)
  qn_batch <- cbind(qn_sample3, qn_sample4) 
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
  qn_kmeans=kmeans(df_tsne_qn, 2)
  
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
  combat_kmeans=kmeans(df_tsne_combat, 2)
  
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
  qn_combat_kmeans=kmeans(df_tsne_qn_combat, 2)
  
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
  total_ncells = 545
  
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




