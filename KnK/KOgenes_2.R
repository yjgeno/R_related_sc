library(scTenifoldKnk)
library(scTenifoldNet)
library(Matrix)

Tcells <- readMM('Tcells.mtx')
rownames(Tcells) <- readLines('genes.Tcells.mtx')
colnames(Tcells) <- readLines('barcodes.Tcells.mtx')  

Tcells <- Tcells[rowMeans(Tcells != 0) >= 0.05,] # filter out low-express genes: select rows mean >= 0.05

WT <- scTenifoldKnk(Tcells, gKO = c('TCF7', 'LEF1'))
save(WT, file = 'Tcells_KO.RData')

#
geneList <- rownames(WT$tensorNetworks$WT)

sapply(geneList, function(gene){
  X <- WT$tensorNetworks$WT # WT GRN
  Y <- X
  Y[gene,] <- 0 # KO each gene including 'TCF7', 'LEF1'
  MA <- manifoldAlignment(X = X,Y = Y, d = 2)
  DR <- scTenifoldKnk:::dRegulation(MA, gKO = gene)
  write.csv(DR, file = paste0('./koProfiles/', gene, '.csv'), quote = FALSE, row.names = FALSE)  
})

# for 'TCF7+LEF1'
# X <- WT$tensorNetworks$WT # WT GRN
# Y <- WT$tensorNetworks$KO
# MA <- manifoldAlignment(X = X,Y = Y, d = 2)
# DR <- scTenifoldKnk:::dRegulation(MA, gKO = c('TCF7', 'LEF1'))
# write.csv(DR, file = paste0('./koProfiles/', 'TCF7+LEF1', '.csv'), quote = FALSE, row.names = FALSE)
