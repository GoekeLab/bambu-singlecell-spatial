library(Seurat)

clusterCells = function(counts, resolution = 0.8, dim = 15){
  
  cellMix <- CreateSeuratObject(counts = counts, 
                                project = "cellMix", min.cells = 1)#, min.features = 200)
  #cellMix <- subset(cellMix, subset = nFeature_RNA > nFeature_RNA_threshold & nFeature_RNA < nFeature_RNA_threshold_max)
  #nFeature_RNA_threshold = 1000, nFeature_RNA_threshold_max = 9000,
  cellMix <- NormalizeData(cellMix, normalization.method = "LogNormalize", scale.factor = 10000)
  cellMix <- FindVariableFeatures(cellMix, selection.method = "vst", nfeatures = 2500)
  all.genes <- rownames(cellMix)
  cellMix <- ScaleData(cellMix, features = all.genes)
  npcs = ifelse(ncol(counts)>50, 50, ncol(counts)-1)
  cellMix <- RunPCA(cellMix, features = VariableFeatures(object = cellMix), npcs = npcs)
  dim = ifelse(dim >= dim(cellMix@reductions$pca)[2], dim(cellMix@reductions$pca)[2], dim)
  cellMix <- FindNeighbors(cellMix, dims = 1:dim)
  cellMix <- FindClusters(cellMix, resolution = resolution)
  cellMix <- RunUMAP(cellMix, dims = 1:dim)
  
  return(cellMix)
}