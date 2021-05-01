library(rPanglaoDB)
library(Seurat)
library(harmony)
library(Nebulosa)
library(Matrix)

Tcells_list <- getMarkers(include = c('CD4', 'IL2RA')) # markers of clusters in DB
counts_Tcells <- getSamples(srs = unique(Tcells_list$SRS), celltype = 'T cells', specie = 'Homo sapiens') # raw counts

set.seed(1)
counts_Tcells <- NormalizeData(counts_Tcells)
counts_Tcells <- FindVariableFeatures(counts_Tcells)
counts_Tcells <- ScaleData(counts_Tcells)
counts_Tcells <- RunPCA(counts_Tcells, verbose = FALSE)
counts_Tcells <- RunHarmony(counts_Tcells, group.by.vars = 'orig.ident')
counts_Tcells <- FindNeighbors(counts_Tcells, reduction = 'harmony')
counts_Tcells <- FindClusters(counts_Tcells)

counts_Tcells <- RunTSNE(counts_Tcells, reduction = 'harmony')
TSNEPlot(counts_Tcells)
plot_density(counts_Tcells, features = c('CD4', 'IL2RA', 'IL2', 'CTLA4'), joint = TRUE)


Tcells_selected <- subset(counts_Tcells, idents = '9')
TSNEPlot(counts_Tcells, cells.highlight=Cells(Tcells_selected))

# check batch effect
counts_Tcells.c <- counts_Tcells
Idents(counts_Tcells.c) <- counts_Tcells$orig.ident
TSNEPlot(counts_Tcells.c)
