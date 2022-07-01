rm(list = ls())

library(Seurat)
library(scTenifoldNet)
library(Matrix)

setwd("data/")
seurat.1 <- readRDS("WT_P27_merged/seurat_p27_processed3.rds")
seurat.2 <- readRDS("WT_P28_merged/seurat_p28_processed3.rds")
seurat_comb <- merge(seurat.1, y = seurat.2, 
                     add.cell.ids = c("p27", "p28"), project = "merge")
View(seurat_comb[[]])
print(dim(seurat_comb))
seurat_comb <- NormalizeData(seurat_comb)
seurat_comb <- FindVariableFeatures(seurat_comb, selection.method = "vst", nfeatures = 3000)
hvg <- seurat_comb@assays[["RNA"]]@var.features

seurat.1 <- seurat.1[hvg, seurat.1$orig.ident=="CSC"]
print(dim(seurat.1))
seurat.2 <- seurat.2[hvg, seurat.2$orig.ident=="CSC"]
print(dim(seurat.2))
exp.1 <- as.matrix(GetAssayData(object = seurat.1, slot = "counts"))
exp.2 <- as.matrix(GetAssayData(object = seurat.2, slot = "counts"))
# saveRDS(exp.1, "exp1.rds")
# saveRDS(exp.2, "exp2.rds")
rm(seurat.1)
rm(seurat.2)
rm(seurat_comb)


# for Matlab
writeMM(obj = exp, file = 'GRNs_p27p28/matrix.mtx')
writeLines(text = rownames(exp), con = 'GRNs_p27p28/features.tsv')
writeLines(text = seurat_comb$new.ident, con = 'GRNs_p27p28/barcodes.tsv')
# metadata <- seurat_comb@meta.data
# write.csv(x = metadata, file = 'metadata.csv', quote = FALSE)


output <- scTenifoldNet(X = exp.2, Y = exp.1,
                        nc_nNet = 10, nc_nCells = 500,
                        td_K = 3, qc_minLibSize = 30)
saveRDS(output, file = "output.rds")

