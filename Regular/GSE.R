# download GSE packed files first

library(Seurat)
library(readr)
library(Matrix)
library(Nebulosa)

counts <- readMM("./GSE167976/matrix.mtx.gz")
genes_id <- read_tsv("./GSE167976/features.tsv.gz", col_names = FALSE)$X1
cells_id <- read_tsv("./GSE167976/barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- genes_id
colnames(counts) <- cells_id

gc(verbose=T)
rm(genes_id)
rm(cells_id)

Tcells <- CreateSeuratObject(counts = counts, project = "Tcells")
#save(Tcells, file = 'Tcells.RData')

meta <- read_tsv("./GSE167976/metadata.tsv.gz", col_names = TRUE)$genotype
Tcells <- AddMetaData(object = Tcells, metadata = meta, col.name = 'genotype')
# select type 'HD' 
Idents(Tcells) <- Tcells@meta.data[["genotype"]]
Tcells <- subset(Tcells, idents = 'HD')

set.seed(1)
Tcells <- NormalizeData(Tcells)
Tcells <- FindVariableFeatures(Tcells)
Tcells <- ScaleData(Tcells)
Tcells <- RunPCA(Tcells, verbose = FALSE)
Tcells <- FindNeighbors(Tcells)
Tcells <- FindClusters(Tcells, resolution = 1.6)

Tcells <- RunTSNE(Tcells)
TSNEPlot(Tcells)

png('density_HD.png', width = 6000, height = 4000, res = 300)
plot_density(Tcells, features = c('CD4', 'IL2RA', 'IL2', 'CTLA4'), joint = TRUE, reduction = 'tsne')
dev.off()

Tcells_selected <- subset(Tcells, idents = '9')
TSNEPlot(Tcells, cells.highlight=Cells(Tcells_selected)) # check 

# check batch effect
Tcells.c <- Tcells
Idents(Tcells.c) <- Tcells$orig.ident # batch ident
TSNEPlot(Tcells.c)
