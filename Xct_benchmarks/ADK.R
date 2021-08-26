library(Seurat)
library(CellChat)
library(ggplot2)

ADK <- readRDS('./ADK.rds')
ADK$saved.idents <- Idents(ADK)
cell.use <- c("B cells", "T cells", "Endothelial cells", "Neutrophil", "Hepatocyte" , "Monocyte" , "Macrophage")

#ADK[['labels']] <- as.vector(Idents(ADK))
ADK <- subset(ADK, idents = cell.use)
ADK$saved.idents <- droplevels(ADK$saved.idents)
WT <- subset(ADK, orig.ident == 'WT')
TG <- subset(ADK, orig.ident == 'TG')
rm(ADK)

#cellchat
cellchat_WT <- createCellChat(object = WT, group.by = 'saved.idents')
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat_WT@DB <- CellChatDB.use 
cellchat_WT <- subsetData(cellchat_WT) 
future::plan("multiprocess", workers = 4) 
cellchat_WT <- identifyOverExpressedGenes(cellchat_WT) 
cellchat_WT <- identifyOverExpressedInteractions(cellchat_WT) 
cellchat_WT <- computeCommunProb(cellchat_WT, raw.use = TRUE)
cellchat_WT <- filterCommunication(cellchat_WT, min.cells = 10)
df.net_WT <- subsetCommunication(cellchat_WT) # all the inferred cell-cell communications at the level of ligands/receptors
#write.csv(x=df.net_WT, file="cellchat_WT.csv")
cellchat_WT <- computeCommunProbPathway(cellchat_WT)

cellchat_WT <- netAnalysis_computeCentrality(cellchat_WT, slot.name = "netP")
cellchat_WT <- aggregateNet(cellchat_WT)

groupSize_WT <- as.numeric(table(cellchat_WT@idents))
png('WT_all_str.png', width = 4000, height = 4000, res = 300)
par(mfrow = c(1,2)) 
netVisual_circle(cellchat_WT@net$count, vertex.weight = groupSize_WT, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", vertex.label.cex = 1.5)
netVisual_circle(cellchat_WT@net$weight, vertex.weight = groupSize_WT, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 1.5)
dev.off()

#TG
cellchat_TG <- createCellChat(object = TG, group.by = "saved.idents")
cellchat_TG@DB <- CellChatDB.use 
cellchat_TG <- subsetData(cellchat_TG) 
cellchat_TG <- identifyOverExpressedGenes(cellchat_TG) 
cellchat_TG <- identifyOverExpressedInteractions(cellchat_TG) 
cellchat_TG <- computeCommunProb(cellchat_TG, raw.use = TRUE)
cellchat_TG <- filterCommunication(cellchat_TG, min.cells = 10)
df.net_TG <- subsetCommunication(cellchat_TG)
#write.csv(x=df.net_TG, file="cellchat_TG.csv")
cellchat_TG <- computeCommunProbPathway(cellchat_TG)

cellchat_TG <- netAnalysis_computeCentrality(cellchat_TG, slot.name = "netP")
cellchat_TG <- aggregateNet(cellchat_TG)

groupSize_TG <- as.numeric(table(cellchat_TG@idents))
png('TG_all_str.png', width = 4000, height = 4000, res = 300)
par(mfrow = c(1,2), xpd=TRUE) 
netVisual_circle(cellchat_TG@net$count, vertex.weight = groupSize_TG, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 1.5)
netVisual_circle(cellchat_TG@net$weight, vertex.weight = groupSize_TG, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 1.5)
dev.off()

rm(WT)
rm(TG)

#pathway in single cellchat object
pathways.show <- c("AGT") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_TG, signaling = pathways.show, layout = "circle",
                    edge.width.max = 10, arrow.size=1)


#merge cellchat
object.list <- list(WT = cellchat_WT, TG = cellchat_TG)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

png('interaction_numbers.png', width = 4000, height = 2000, res = 300)
compareInteractions(cellchat, show.legend = F, group = c(1,2), size.text = 20)
#gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
dev.off()


#cell-level interactions: circle plot
weight.max <- getMaxWeight(object.list[1:2], slot.name = c("net"), attribute = 'weight')
cell.max <- getMaxWeight(object.list[1:2], slot.name = 'idents', attribute = 'idents')
png('groups.png', width = 4000, height = 2000, res = 300)
par(mfrow = c(1,2)) 
for (i in 1:length(object.list)) {
  netVisual_circle(cellchat_TG@net$count, vertex.weight = as.numeric(table(object.list[[i]]@idents)), weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[1], edge.width.max = 0.4, vertex.weight.max = cell.max[1],
                   title.name = "Number of interactions", vertex.label.cex = 2.5)
}
dev.off()

for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, vertex.weight = as.numeric(table(object.list[[i]]@idents)), weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[1], vertex.weight.max = cell.max[1],
                   title.name = "Interaction weights/strength", vertex.label.cex = 2.5)
}

#cell-level interactions: alternative heatmap: cell vs cell
png('heatmap.png', width = 3000, height = 3000, res = 300)
netVisual_heatmap(cellchat, measure = "weight", font.size = 16, font.size.title = 16)
dev.off()

#diff interaction numbers and strength
png('interaction_diff_.png', width = 4000, height = 4000, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, vertex.label.cex = 2.5, arrow.size = 1)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", vertex.label.cex = 2.5, arrow.size = 1)
dev.off()

#bubble plot: cell-cell vs L-R
png('cell_bubble_macro.png', width = 3000, height = 2000, res = 300)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(7, 3, 5, 6),  comparison = c(1, 2), angle.x = 45,
                 font.size = 16, font.size.title = 16)
netVisual_bubble(cellchat, sources.use = c(5, 6, 7), targets.use = 6,  comparison = c(1, 2), angle.x = 45,
                 font.size = 16, font.size.title = 16)
dev.off()

#diff of in/out signaling of each cell type
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, font.size = 16,
                                               font.size.title = 16)
}
png('interaction_diff_str_wrap.png', width = 4000, height = 2000, res = 300)
patchwork::wrap_plots(plots = gg)
dev.off()


#pathway-level
png('pathways_diff.png', width = 3000, height = 5000, res = 300)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, font.size = 20, color.use = c('#377EB8', '#E41A1C'))
dev.off()

#for a specific pathway diff
pathways.show <- c("AGT")
weight.max <- getMaxWeight(object.list[1:2], slot.name = c("netP"), attribute = pathways.show)

png('PARs.png', width = 4000, height = 2000, res = 300)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, arrow.size=1,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
} #color.use = color.use.list[[i]]
dev.off()





