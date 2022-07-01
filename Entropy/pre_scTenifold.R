library(Seurat);
library(stringr);

setwd("data/")
########### p27
exp <- read.table(file = "WT_P27_merged/X.txt", 
                  header = F,
                  sep = ",",
                  as.is = TRUE) # raw
gene <- read.table(file = "WT_P27_merged/g.txt", 
                   header = F,
                   sep = ",",
                   as.is = TRUE)$V1
cell.type <- read.table(file = "WT_P27_merged/c_cell_type_tx.txt", 
                        header = F,
                        sep = ",",
                        as.is = TRUE)
cell.type[c('cell_type', 'cell_cluster')] <- str_split_fixed(cell.type$V1, '_', 2)
cell.type <- cell.type$cell_type

# Create Seurat object
seurat_p27 <- CreateSeuratObject(counts = exp, row.names = gene)
seurat_p27$orig.ident <- cell.type
rm(exp)

########### p28
exp <- read.table(file = "WT_P28_merged/X.txt", 
                  header = F,
                  sep = ",",
                  as.is = TRUE) # raw
gene <- read.table(file = "WT_P28_merged/g.txt", 
                   header = F,
                   sep = ",",
                   as.is = TRUE)$V1
cell.type <- read.table(file = "WT_P28_merged/c_cell_type_tx.txt", 
                        header = F,
                        sep = ",",
                        as.is = TRUE)
cell.type[c('cell_type', 'cell_cluster')] <- str_split_fixed(cell.type$V1, '_', 2)
cell.type <- cell.type$cell_type

# Create Seurat object
seurat_p28 <- CreateSeuratObject(counts = exp, row.names = gene)
seurat_p28$orig.ident <- cell.type
rm(exp)

# preprocessing
genes_touse_p27p28 <- intersect(rownames(seurat_p27), rownames(seurat_p28))
# write.table(genes_touse_p27p28, file = "genes_touse_p27p28.txt", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = F)
# filtering ['.', 'Rik', 'Gm', 'Rps', 'Rpl', 'mt-']
genes_touse <- read.table(file = "genes_touse_filtered_p27p28.txt", 
                          header = F,
                          sep = "\t",
                          as.is = TRUE)$V1
seurat_p27 <- seurat_p27[genes_touse, ]
seurat_p28 <- seurat_p28[genes_touse, ]
print(dim(seurat_p27))
print(dim(seurat_p28))

saveRDS(seurat_p27, file = "seurat_p27_processed3.rds")
saveRDS(seurat_p28, file = "seurat_p28_processed3.rds")
