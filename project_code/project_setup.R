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

