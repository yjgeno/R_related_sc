rm(list = ls())

library(Seurat)
library(stringr)

setwd("data/")
seurat_processed <- readRDS("WT_P21_merged/seurat_p21_processed.rds")

# Cell cycle
# A list of cell cycle markers
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_processed <- CellCycleScoring(seurat_processed,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     set.ident = TRUE)
seurat_processed <- seurat_processed[ ,seurat_processed$Phase=="G1"]
print(dim(seurat_processed))
#
exp <- as.matrix(GetAssayData(object = seurat_processed, slot = "counts"))
print(dim(exp))
cell_type <- as.vector(seurat_processed$orig.ident)
var_genes_p21p23 <- toupper(row.names(seurat_processed))
rm(seurat_processed)

# write.table(cell_type, file = "celltype_p21_G1.txt", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = F)

library(org.Hs.eg.db);
hs <- org.Hs.eg.db
id <- select(hs, 
       keys = var_genes_p21p23,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL",
       multiVals = "first")
dup <- as.numeric(rownames(id[duplicated(id$SYMBOL),])) # if duplicate
print(dup)
id <- id[!rownames(id) %in% dup, ]$ENTREZID
rownames(exp) <- id
colnames(exp) <- cell_type
exp <- exp[!is.na(rownames(exp)), ]
print(dim(exp))


# Entropy
library(SCENT)
# library(devtools)
# remotes::install_github("aet21/SCENT",
#                         INSTALL_opts = c("--no-multiarch")
# )

# PPI network
data(net13Jun12);
print(dim(net13Jun12.m));
# write.table(net13Jun12.m, file = "net13Jun12.txt", sep = ",",
#             row.names = T, col.names = T, quote = F)

# net <- as.matrix(read.table(file = "net13Jun12.txt", 
#                             header = T,
#                        sep = ",",
#                        as.is = TRUE))
# colnames(net) <- sub("X", "", colnames(net))


# example
data(dataChu)
print(dim(scChu.m)); # genes by cells
print(summary(factor(phenoChu.v))); # cells
# for PPI
lscChu.m <- log2(scChu.m+1.1);
range(lscChu.m);
# for CCAT
lscChu0.m <- log2(scChu.m+1);
range(lscChu0.m);
# PPI
integ.l <- DoIntegPPI(exp.m = lscChu.m, ppiA.m = net13Jun12.m)
str(integ.l)
sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 1)
boxplot(srChu.v ~ phenoChu.v, main = "SR potency estimates", 
        xlab = "Cell Type", 
        ylab = "SR")
# CCAT
ccat.v <- CompCCAT(exp = lscChu0.m, ppiA = net13Jun12.m)
boxplot(ccat.v ~ phenoChu.v, main = "SR potency estimates", 
        xlab = "Cell Type", ylab = "SR")


# real dataset
print(dim(exp));
# PPI
exp.m <- log2(exp+1.1);
range(exp.m);
integ.l <- DoIntegPPI(exp.m = exp.m, ppiA.m = net13Jun12.m)
# str(integ.l)
sr.o <- CompSRana(integ.l, local = FALSE, mc.cores = 1)
saveRDS(sr.o, file = "sr.o.rds")

boxplot(sr.o$SR ~ cell.type, main = "SR potency estimates", 
        xlab = "Cell Type", 
        ylab = "SR") # signaling entropy rate
write.table(sr.o$SR, file = "SR.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE)

# CCAT
exp0.m <- as.matrix(log2(exp+1));
range(exp0.m);

ccat <- CompCCAT(exp = exp0.m, ppiA = net13Jun12.m);
boxplot(ccat ~ cell_type, main = "SR potency estimates", 
        xlab = "Cell Type", ylab = "SR")
write.table(ccat, file = "ccat_p23_G0.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = F)










