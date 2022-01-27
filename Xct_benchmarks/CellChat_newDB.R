library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(Seurat)
options(stringsAsFactors = FALSE)
setwd("")

seuratObject <- readRDS('humanSkin.rds')
HVG <- VariableFeatures(seuratObject)
s <- s[HVG,]

# cellchat object
cellchat <- createCellChat(object = s, group.by = "ident") # group.by to define cell groups
#groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# curate database
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
# write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")

interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) # subset to signaling genes for saving computation cost
dim(cellchat@data.signaling)

# future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat, 
                                       thresh.fc = -100,
                                       thresh.p = 100) # identify over-expressed ligands or receptors in one cell group
DE <- cellchat@var.features$features.info#[[paste0(features.name, ".info")]] # DE result

# or assign own DE genes instead of identifyOverExpressedGenes
cellchat@var.features[["features"]] <- rep(rownames(seuratObject), 2)

# vst.info <- cellchat@var.features$features.info
# write.csv(vst.info, file = "Z:/Cailab/scTenifoldXct/QianEdit/vst_CellChat.csv")
# write.csv(rownames(seuratObject), file = "Z:/Cailab/scTenifoldXct/QianEdit/genes_CellChat.csv")
# write.csv(rownames(seuratObject), file = "genes_CellChat.csv")
vst.info1 <- read.csv(file = 'vst_CellChat_processed.csv', row.names = 1)
cellchat@var.features$features.info <- vst.info1


length(cellchat@var.features$features) # total DE genes
length(unique(cellchat@var.features$features))
sum(cellchat@var.features$features.info$clusters == 'Inflam. FIB') # cluster-specific DE genes
sum(cellchat@var.features$features.info$clusters == 'Inflam. DC')

cellchat <- identifyOverExpressedInteractions(cellchat) # then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed
dim(cellchat@LR$LRsig)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat) # for pairs in @LR 
dim(cellchat@net$prob)

# cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat, thresh = 1000) # all the inferred cell-cell communications at the level of ligands/receptors
dim(df.net)[1]

# cellchat <- computeCommunProbPathway(cellchat)
df.net <- df.net[df.net$source == 'Inflam. FIB' & df.net$target =='Inflam. DC',]
dim(df.net)[1]

length(which(cellchat@net$prob[1:1,2:2,] != 0)) # probs (LR scores of FIB-DC), others are 0

# write.csv(df.net, file = "CellChat_result_0114.csv") # 44 using all genes
