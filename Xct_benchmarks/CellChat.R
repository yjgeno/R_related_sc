library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

load(url("https://ndownloader.figshare.com/files/25950872")) # This is a combined data from two biological conditions: normal and diseases
data.input = data_humanSkin$data # input normalized data matrix
meta = data_humanSkin$meta # a dataframe of cell metadata
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels") # group.by to define cell groups
#groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset to signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat) # identify over-expressed ligands or receptors in one cell group
cellchat <- identifyOverExpressedInteractions(cellchat) # then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed
cellchat <- projectData(cellchat, PPI.human) #project gene expression data onto protein-protein interaction (PPI) network

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat) # all the inferred cell-cell communications at the level of ligands/receptors
cellchat <- computeCommunProbPathway(cellchat)
#groupSize <- as.numeric(table(cellchat@idents))


# visualization
# aggregated 
cellchat <- aggregateNet(cellchat)
par(mfrow = c(1,2), xpd=TRUE) # graphical parameters
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# for specific pathway/L-R pair
pathways.show <- c("CXCL") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # slot 'netP' for pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# communication pattern
selectK(cellchat, pattern = "outgoing") # 3 drop 
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing") # river plot
