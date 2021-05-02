library(Seurat)
pbmc <- pbmc_small # class Seurat

# Basic Attributes
dim(pbmc)
nrow(pbmc)
ncol(pbmc)
rownames(pbmc)
colnames(pbmc)
names(pbmc) # OO objects which start by $: "RNA","RNA_snn","pca","tsne" 

Idents(pbmc) # a col in cell-level metadata
pbmc$saved.idents <- Idents(pbmc) # saved to meta.data

HVFInfo(pbmc) # dispersion info for all genes, dataframe
VariableFeatures(pbmc) # HVG


# access data
rna <- pbmc[['RNA']] # == pbmc@assays[["RNA"]] == pbmc$RNA, class Assay
counts <- GetAssayData(pbmc, slot = 'scale.data') # == GetAssayData(rna, slot = 'scale.data'), matrix


# gene-level metadata
meta.genes <- rna[[]]
colnames(meta.genes)

# cell-level metadata
meta <- pbmc[[]] # dataframe
colnames(meta)
meta$groups == meta[['groups']] # == meta[[6]] # access cols
groups <- pbmc[['groups', drop = TRUE]] # drop: when only one col left after subset, turn dataframe into vector


# subset
# 1: by cluster
pbmc_selected1 <- subset(pbmc, idents = '1')
# 2: by cells (metadata)
cells.use <- colnames(pbmc)[which(meta['groups'] == 'g1')]
pbmc_selected2 <- subset(pbmc, cells = cells.use)
# 3: by genes (metadata)
HVG <- rownames(pbmc)[which(meta.genes$vst.variable == TRUE)]
pbmc_selected3 <- subset(pbmc, features = HVG)
# 4: by gene expression
Key(rna) # "rna_"
rna_CD9 <- FetchData(object = pbmc, vars = 'rna_CD9') # single gene count, dataframe
pbmc_selected4 <- pbmc[, which(rna_CD9 > 2 & rna_CD9 < 6)] #or

pbmc_selected4 <- subset(pbmc, subset = CD9 > 2 & CD9 < 6)


