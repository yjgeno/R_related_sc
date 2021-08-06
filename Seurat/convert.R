library(Seurat)
library(SeuratDisk)

setwd('./data')

#h5ad to seurat
Convert("hca_heart_immune_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
seuratObj <- LoadH5Seurat("hca_heart_immune_raw.h5seurat")

#view counts/data/scale.data
mat <- GetAssayData(object = heart, slot = "data")
df <- as.data.frame(as.matrix(mat)) 

#convert cell annotations from factors to chara
seuratObj[['cell_states']] <- lapply(seuratObj[['cell_states']], as.character)

seuratObj <- NormalizeData(seuratObj) 
seuratObj <- FindVariableFeatures(seuratObj)

#seurat to h5ad
SaveH5Seurat(seuratObj, filename = "hca_heart_immune.h5Seurat")
Convert("hca_heart_immune.h5Seurat", dest = "h5ad", assays = list(SCT = c("counts", "data")))


