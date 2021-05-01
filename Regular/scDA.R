library(Matrix)
library(Seurat)
library(ggplot2)
library(harmony)
library(fgsea)
library(pbapply)
library(ggrepel)
library(enrichR)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
markerGenes <- gmtPathways('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/mmuPanglaoDB.gmt')

set.seed(1)
WT <- Read10X_h5('liver_WT/filtered_feature_bc_matrix.h5')
Tg <- Read10X_h5('liver_Tg/filtered_feature_bc_matrix.h5')
WT <- CreateSeuratObject(WT, project = 'WT')
Tg <- CreateSeuratObject(Tg, project = 'Tg')
ALL <- merge(WT, Tg)
rm(WT)
rm(Tg)
#gc()

ALL <- scQC(ALL)
ALL <- NormalizeData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- ScaleData(ALL)
ALL <- RunPCA(ALL)
#ALL <- RunTSNE(ALL)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunTSNE(ALL, reduction = 'harmony')
Idents(ALL) <- ALL$orig.ident

#png('sampleCells.png', width = 1500, height = 1000, res = 300)
#TSNEPlot(ALL) + xlab('t-SNE 1') + ylab('t-SNE 2') + theme_bw()
#dev.off()

ALL_orig <- ALL
ALL <- FindNeighbors(ALL, reduction = 'tsne', dims = 1:2)
ALL <- FindClusters(ALL, resolution = 0.001)
aDE <- FindAllMarkers(ALL)

cellTypes <- pbsapply(levels(aDE$cluster), function(X){
  tDE <- aDE[aDE$cluster %in% X,]
  FC <- tDE$avg_log2FC
  names(FC) <- tDE$gene
  E <- fgsea(markerGenes, FC, eps = 0)
  E <- E[E$padj < 0.05 & E$NES > 0,]
  E <- E[order(E$NES, decreasing = TRUE),]
  E$pathway[1]
})
levels(Idents(ALL)) <- cellTypes

#png('cellTypes.png', width = 1600, height = 1000, res = 300)
#TSNEPlot(ALL) + xlab('t-SNE 1') + ylab('t-SNE 2') + theme_bw()
#dev.off()

Idents(ALL) <- paste0(ALL$orig.ident, '_', Idents(ALL))
