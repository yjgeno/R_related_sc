rm(list = ls())

library(Seurat)
library(stringr)

setwd("data/")
########### p21
exp <- read.table(file = "WT_P21_merged/X.txt", 
                  header = F,
                  sep = ",",
                  as.is = TRUE) # raw
gene <- read.table(file = "WT_P21_merged/g.txt", 
                   header = F,
                   sep = ",",
                   as.is = TRUE)$V1
cell.type <- read.table(file = "WT_P21_merged/c_cell_type_tx.txt", 
                        header = F,
                        sep = ",",
                        as.is = TRUE)
cell.type[c('cell_type', 'cell_cluster')] <- str_split_fixed(cell.type$V1, '_', 2)
cell.type <- cell.type$cell_type

# Create Seurat object
seurat_p21 <- CreateSeuratObject(counts = exp, row.names = gene)
seurat_p21$orig.ident <- cell.type
rm(exp)

########### p23
exp <- read.table(file = "WT_P23_merged/X.txt", 
                  header = F,
                  sep = ",",
                  as.is = TRUE) # raw
gene <- read.table(file = "WT_P23_merged/g.txt", 
                   header = F,
                   sep = ",",
                   as.is = TRUE)$V1
cell.type <- read.table(file = "WT_P23_merged/c_cell_type_tx.txt", 
                        header = F,
                        sep = ",",
                        as.is = TRUE)
cell.type[c('cell_type', 'cell_cluster')] <- str_split_fixed(cell.type$V1, '_', 2)
cell.type <- cell.type$cell_type

# Create Seurat object
seurat_p23 <- CreateSeuratObject(counts = exp, row.names = gene)
seurat_p23$orig.ident <- cell.type
rm(exp)

###########
genes_touse_p21p23 <- intersect(rownames(seurat_p21), rownames(seurat_p23))
# write.table(genes_touse_p23p27, file = "genes_touse_p23p27.txt", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = F)
# filtering ['.', 'Rik', 'Gm', 'Rps', 'Rpl', 'mt-']

genes_touse <- read.table(file = "genes_touse_filtered_p21p23.txt", 
                        header = F,
                        sep = "\t",
                        as.is = TRUE)$V1
seurat_p21 <- seurat_p21[genes_touse, ]
seurat_p23 <- seurat_p23[genes_touse, ]
seurat_comb <- merge(seurat_p21, y = seurat_p23, 
                     add.cell.ids = c("p21", "p23"), project = "merge")
print(dim(seurat_comb))

seurat_comb <- NormalizeData(seurat_comb)
seurat_comb <- FindVariableFeatures(seurat_comb, selection.method = "vst", nfeatures = 5000)
var_genes_p21p23 <- VariableFeatures(seurat_comb)
# write.table(var_genes_p21p23, file = "var_genes_p21p23.txt", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = F)
seurat_comb <- ScaleData(seurat_comb, features = rownames(seurat_comb))

seurat_p21 <- seurat_comb[, 1:dim(seurat_p21)[2]]
seurat_p23 <- seurat_comb[, 1+dim(seurat_p21)[2]:dim(seurat_comb)[2]]
print(dim(seurat_p21))
print(dim(seurat_p23))
rm(seurat_comb)

saveRDS(seurat_p21, file = "seurat_p21_processed.rds")
saveRDS(seurat_p23, file = "seurat_p23_processed.rds")




