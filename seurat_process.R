.libPaths( c( "/share/hormozdiarilab/Experiments/earlyPredictionCNV/R/x86_64-pc-linux-gnu-library/4.0/" , .libPaths() ) )

library(dplyr)
library(Seurat)  # 3.2.1, but changed to 2.3.4
library(patchwork)
library(Matrix)

load("/share/hormozdiarilab/Data/scRNASeq_AllenBrain/codex/raw_counts_mat.rdata")

input_data <- CreateSeuratObject(raw_counts_mat, project = "Polioudakis", min.cells = 3, min.genes = 200)

input_data <- NormalizeData(
  input_data, normalization.method = "LogNormalize", scale.factor = 10000)

input_data <- FindVariableGenes(input_data)

# all.genes <- rownames(input_data)
input_data <- ScaleDataR(input_data)

input_data <- RunPCA(input_data, pc.genes = input_data@var.genes, pcs.compute = 40)
# input_data <- FindNeighbors(input_data, dims = 1:40)

input_data = FindClusters(input_data, dims.use = 1:40, save.SNN=TRUE)
# print(input_data)
# write.table(as.matrix(input_data@snn), sep=",", file="RNA_snn.csv", quote=FALSE)

# doing t-SNE just for confirmation that the distance matrix will be the same
input_data = RunTSNE(input_data, dims.use = 1:40)
# TSNEPlot(input_data, do.label = TRUE, pt.size = 1, label.size = 4, 
#          vector.friendly = TRUE, png.file = "plotting.png")

write.table(as.matrix(input_data@ident), sep=",", file="cell_identity.csv", quote=FALSE)
# print(WhichCells(input_data, cells.use))

# write.table(as.matrix(input_data$RNA_snn), sep=",", file="RNA_snn.csv", quote=FALSE)
# input_data <- FindClusters(input_data, resolution = 0.5, save.SNN = TRUE)
# print(input_data)