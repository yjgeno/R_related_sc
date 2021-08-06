library(Seurat)

write.table(as.matrix(GetAssayData(object = ALL[["RNA"]], slot = "scale.data")), 
            './cellphoneDB/cpdb_counts.txt', sep = '\t', quote = F)

# generating meta file
Idents(ALL) <- sub('.*_\\s*', '', Idents(ALL))
labels <- Idents(ALL)
meta <- data.frame(cell = names(labels), cell_type = labels) 
write.table(meta, './cellphoneDB/cpdb_meta.txt', sep='\t', quote=F, row.names=F)
