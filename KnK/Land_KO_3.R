library(RSpectra)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(GSVA)
library(fgsea)

GOBP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')

koList <- list.files(path = 'koProfiles/', full.names = TRUE)
koProfiles <- lapply(koList, function(X){
  koProfile <- read.csv(X)
  Z <- koProfile$Z # Z-score
  names(Z) <- koProfile$gene # gene names
  return(Z)  
}) # list

allGenes <- unique(unlist(lapply(koProfiles, names)))

koProfiles <- sapply(koProfiles, function(X){X[allGenes]}) # to matrix

colnames(koProfiles) <- gsub('.csv','', basename(koList))
rownames(koProfiles) <- toupper(allGenes)
write.csv(koProfiles, 'allKOProfiles.csv')

################################################################
E <- gsva(koProfiles, gset.idx.list = GOBP, method = 'ssgsea')
E <- round(E, 2)
write.csv(E, 'gobpEnrichment.csv')
SS <- E[grepl('shear', rownames(E)),]
write.csv(SS, 'ssEnrichment.csv')
################################################################

# tsne on PCA 50
covKO <- cov(koProfiles)
pcaKO <- eigs(covKO, 50)$vectors
rownames(pcaKO) <- colnames(koProfiles)

set.seed(42)
tsneKO <- Rtsne(pcaKO)$Y # 2d coordinates
rownames(tsneKO) <- colnames(koProfiles)

# keeping original data
write.csv(tsneKO, 'tSNEcoordinates.csv')
# tsneKO_orig <- tsneKO
# tsneKO <- tsneKO_orig

tsneKO <- as.data.frame(tsneKO)
# Creating k-means 
fit_cluster_kmeans=kmeans(tsneKO, 30)  
tsneKO$kcluster <-  factor(fit_cluster_kmeans$cluster)

tsneKO$gene <- colnames(koProfiles)
#write.csv(tsneKO, 'tSNEcoordinates.csv')
tsneKO$gene[!tsneKO$gene %in% c('TCF7', 'LEF1', 'TCF7+LEF1', 'LAG3', 'CTLA4')] <- NA # only gene labels, others NA

png('landscape_k11.png', width = 2000, height = 2000, res = 300)
ggplot(tsneKO, aes(V1,V2, label = gene)) + 
  geom_point(alpha=ifelse(is.na(tsneKO$gene) == TRUE, 0.2, 1),
             color=ifelse(tsneKO$kcluster == tsneKO ['TCF7+LEF1',]$kcluster, 'red', 'black'),
             size=ifelse(is.na(tsneKO$gene) == TRUE, 1, 2)) + 
  theme_bw() + 
  xlab('t-SNE 1') + 
  ylab('t-SNE 2') + 
  geom_text_repel(nudge_x = 5, nudge_y = 5) + 
  labs(title = 'TCF7+LEF1') +
  theme(plot.title = element_text(face = 2))
dev.off()

# get cluster genes
target_set <- rownames(tsneKO[tsneKO$kcluster == tsneKO['TCF7+LEF1',]$kcluster, ])
write.table(target_set, 'target_geneset.csv', col.names = FALSE, row.names = FALSE)
