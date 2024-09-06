#if("Seurat" %in% rownames(installed.packages()) == FALSE) {remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)}
#install.packages("SeuratObject")
#install.packages("Seurat")
library(Seurat)

clusterCells = function(counts, resolution = 0.8, dims = 1:15){
  
  cellMix <- CreateSeuratObject(counts = counts, 
                                project = "cellMix", min.cells = 1)#, min.features = 200)
  #cellMix <- subset(cellMix, subset = nFeature_RNA > nFeature_RNA_threshold & nFeature_RNA < nFeature_RNA_threshold_max)
  #nFeature_RNA_threshold = 1000, nFeature_RNA_threshold_max = 9000,
  cellMix <- NormalizeData(cellMix, normalization.method = "LogNormalize", scale.factor = 10000)
  cellMix <- FindVariableFeatures(cellMix, selection.method = "vst", nfeatures = 2500)
  all.genes <- rownames(cellMix)
  cellMix <- ScaleData(cellMix, features = all.genes)
  cellMix <- RunPCA(cellMix, features = VariableFeatures(object = cellMix))
  dim = dim(cellMix@reductions$pca)[2]
  if(dim < 15){ dims = 1:dim}
  cellMix <- FindNeighbors(cellMix, dims = dims)
  cellMix <- FindClusters(cellMix, resolution = resolution)
  cellMix <- RunUMAP(cellMix, dims = dims)
  
  return(cellMix)
}