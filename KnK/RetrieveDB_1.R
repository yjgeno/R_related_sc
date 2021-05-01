library(rPanglaoDB)
library(Seurat)
library(harmony)
library(Matrix)
library(Nebulosa)

#BiocManager::install("harmony")
#library(devtools)
#install_github("immunogenomics/harmony")

csv_names <- list.files('./CSV_source', full.names=TRUE) 
dbQuery = lapply(csv_names,read.csv, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dbQuery = rbind (dbQuery[[1]], dbQuery[[2]], dbQuery[[3]])


#dbQuery <- read.csv('./CSV_source/CD4+.csv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
dbQuery <- dbQuery[grepl('Mus', dbQuery[,1]),] # grep Mus in first column

SRS <- dbQuery[,3]
SRS <- unique(unlist(lapply(strsplit(SRS, '\\_'), function(X){X[2]})))
View(SRS)


downloadedData <- lapply(SRS, function(X){
  try(getSamples(srs = X, celltype = 'T cells'))
})

# filter out Seurat objects
downloadedData <- downloadedData[unlist(lapply(downloadedData, class)) %in% 'Seurat']

# sum all to [[1]], [-1]: loop over seq except index 1
for(i in seq_along(downloadedData)[-1]){
  downloadedData[[1]] <- merge(downloadedData[[1]], downloadedData[[i]])
}

#for(i in seq_along(downloadedData)[-1]) print(i)
Tcells <- downloadedData[[1]]

#Tcells$orig.ident[1] # SRS No. labeled with Cell Barcodes

Tcells <- NormalizeData(Tcells)
Tcells <- FindVariableFeatures(Tcells)
Tcells <- ScaleData(Tcells)
Tcells <- RunPCA(Tcells) # which genes driving the new PCs
Tcells <- RunHarmony(Tcells, group.by.vars = 'orig.ident') # group.by.vars - effect of var. you want to remove, here by "SRS No."
Tcells <- RunTSNE(Tcells, reduction = 'harmony') # dimensional reduction to use for the tSNE. Default is PCA
TSNEPlot(Tcells)
FeaturePlot(Tcells, features = c('CD4', 'PDCD1', 'PTPRC'), order = TRUE)
png('./markers_Tcells.png', width = 3000, height = 2000, res = 300)
plot_density(Tcells, features =  c('CD4', 'PDCD1', 'PTPRC'), joint = TRUE) #nebulosa

dev.off() # the plot isn't written to the file until you call dev.off
Tcells <- FindNeighbors(Tcells, reduction = 'tsne', dims = 1:2) #dims of reduction to use to build SNN
Tcells <- FindClusters(Tcells, resolution = 0.01)
TSNEPlot(Tcells)
Tcells_selected <- subset(Tcells, idents = '1')

countMatrix <- Tcells_selected@assays$RNA@counts
writeMM(countMatrix, 'Tcells.mtx')
writeLines(rownames(countMatrix), 'genes.Tcells.mtx')
writeLines(colnames(countMatrix), 'barcodes.Tcells.mtx')
