require(liana)
require(tibble)
require(purrr)
library(Seurat)
library(readxl)
library(VennDiagram)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplot2)

setwd("")
seuratObject <- readRDS('humanSkin.rds')
HVG <- VariableFeatures(seuratObject)
s <- seuratObject[HVG,]

dim(s)
# 
# squidpy <- liana_wrap(seuratObject,resource = 'OmniPath',method = 'squidpy')

# check default setting
liana_defaults(call_italk.params = list(.DE = F))$call_italk


setwd("/liana")
xct_omni <- read.csv('omnipath_intercell_toUse_forLIANA.csv') # xct omnipath DB

# natmi 
call_natmi <- liana_wrap(s, resource = 'OmniPath', method = 'call_natmi')
call_natmi <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'call_natmi')
call_natmi <- call_natmi[call_natmi$source == 'Inflam. FIB' & call_natmi$target =='Inflam. DC',] # 102
write.csv(call_natmi,file = 'call_natmi_s_xctDB2.csv')

# cytotalk
cytotalk <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'cytotalk')
cytotalk <- cytotalk[cytotalk$source == 'Inflam. FIB' & cytotalk$target =='Inflam. DC',]
write.csv(cytotalk,file = 'cytotalk_s_xctDB2.csv')

# sca
call_sca <- liana_wrap(s,resource = 'OmniPath',method = 'call_sca')
call_sca <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'call_sca')
call_sca <- call_sca[call_sca$source == 'Inflam. FIB' & call_sca$target =='Inflam. DC',] # 80
write.csv(call_sca,file = 'call_sca_s_xctDB2.csv')


# cellchat 
cellchat <- read.csv("CellChat_result_0114_HVG.csv")

cellchat <- liana_wrap(s, resource = 'OmniPath', method = 'cellchat',
                       cellchat.params = 
                         list(
                           nboot = 10,
                           exclude_anns = NULL,
                           thresh = 10000, # permutation p-val
                           de_thresh = 100000, # DE p-val
                           assay = "RNA",
                           .normalize = FALSE,
                           .do_parallel = FALSE,
                           .raw_use = TRUE,
                           sources.use = 'Inflam. FIB',
                           target.use = 'Inflam. DC'
                         ))

cellchat <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'cellchat', 
                       cellchat.params = 
                         list(
                           nboot = 100,
                           .format = T,
                           thresh = 10, # permutation p-val
                           de_thresh = 10 # DE p-val
                         ))
cellchat <- cellchat[cellchat$source == 'Inflam. FIB' & cellchat$target =='Inflam. DC',]
write.csv(cellchat,file = 'cellchat.csv')


# italk
call_italk <- liana_wrap(s, resource = 'OmniPath', method = 'call_italk')

call_italk <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'call_italk', 
                         call_italk.params = list(assay = "RNA", .DE = TRUE, logfc.threshold = 0))
call_italk <- call_italk[call_italk$source == 'Inflam. FIB' & call_italk$target =='Inflam. DC',]


# alternative italk
call_italk <- liana_pipe(
  s,
  op_resource = xct_omni,
  decomplexify = F,
  test.type = "wilcox",
  pval.type = "all",
  trim = 0,
  assay = "RNA",
  assay.type = "logcounts"
)

# temp <- call_italk
# call_italk <- temp

call_italk <- call_italk[call_italk$ligand.log2FC>0 & call_italk$receptor.log2FC>0, ]
call_italk <- call_italk[call_italk$source == 'Inflam. FIB' & call_italk$target =='Inflam. DC',]
call_italk$logfc_comb <- call_italk$ligand.log2FC * call_italk$receptor.log2FC
write.csv(call_italk,file = 'call_italk_s_xctDB2.csv')

# connectome
call_connectome <- liana_wrap(s,resource = 'OmniPath',method = 'call_connectome')
call_connectome <- liana_wrap(s, resource='custom', external_resource = xct_omni, method = 'call_connectome')

call_connectome <- call_connectome[call_connectome$source == 'Inflam. FIB' & call_connectome$target =='Inflam. DC',] # 388
# call_connectome <- call_connectome[call_connectome$p_val_adj.lig < 0.05 & call_connectome$p_val_adj.rec < 0.05,]
write.csv(call_connectome,file = 'call_connectome_s_xctDB2.csv')


# nichenet
nichienet <- read.csv('nichienet.csv')
nichienet$cat <- paste(nichienet$from,'&',nichienet$to)


# Xct
xct_res <- read.csv("LS_skin_FIB2DC_full_grace.csv") # full genes
xct_res <- read.csv("LS_skin_FIB2DC_full.csv") # 160; HVG 2.5% default
xct_res <- read.csv("LS_skin_FIB2DC_5.csv") # HVG 5%

xct_res$cat <- paste(xct_res$ligand,'-', xct_res$receptor)

# analysis
call_sca <- read.csv('call_sca_s_xctDB2.csv')
cellchat <- read.csv('CellChat_result_0114_HVG.csv')
call_italk <- read.csv('call_italk_s_xctDB2.csv')
call_connectome <- read.csv('call_connectome_s_xctDB2.csv') 
call_natmi <- read.csv('call_natmi_s_xctDB2.csv')
cytotalk <- read.csv('cytotalk_s_xctDB2.csv')

call_sca$cat <- paste(call_sca$ligand,'-',call_sca$receptor)
cellchat$cat <- paste(cellchat$ligand,'-',cellchat$receptor)
call_italk$cat <- paste(call_italk$ligand,'-',call_italk$receptor)
call_connectome$cat <- paste(call_connectome$ligand,'-',call_connectome$receptor)
call_natmi$cat <- paste(call_natmi$ligand,'-',call_natmi$receptor)
cytotalk$cat <- paste(cytotalk$ligand,'-',cytotalk$receptor)


n = 30
call_connectome_trim <- call_connectome[order(call_connectome$weight_sc, decreasing = T),][1:n,]
call_natmi_trim <- call_natmi[order(call_natmi$edge_specificity, decreasing = T),][1:n,]
call_sca_trim <- call_sca[order(call_sca$LRscore, decreasing = T),][1:n,]
call_italk_trim <- call_italk[order(call_italk$logfc_comb, decreasing = T),][1:n,]
cytotalk_trim <- cytotalk[order(cytotalk$crosstalk_score, decreasing = T),][1:n,]
if (n > length(cellchat)) {cellchat_trim <- cellchat[order(cellchat$prob, decreasing = T),][1:n,]} else {cellchat_trim <- cellchat}

############ upset plot ###########
library('UpSetR')
require(ggplot2); require(plyr); require(gridExtra); require(grid);

listInput <- list(XCT = xct_res$cat[1:n],
                  Cellchat = cellchat_trim$cat,
                  italk=call_italk_trim$cat,
                  SCA = call_sca_trim$cat,
                  CONNECTOME = call_connectome_trim$cat,
                  NATMI = call_natmi_trim$cat,
                  cytotalk = cytotalk_trim$cat)

upset(fromList(listInput), order.by = "freq", nsets = 7)

pdf("upsetplot_0116_top30.pdf",         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")          # Paper size

dev.off()


############ venn diagram ###########
set1 <- xct_res$cat[1:n]
set2 <- call_sca_trim$cat
set3 <- cellchat_trim$cat
set4 <- call_italk_trim$cat
set5 <- call_connectome_trim$cat
set6 <- call_natmi_trim$cat
set7 <- cytotalk$cat

x <- list(scTenifoldXct = set1, CytoTalk = set7)
png('venn_chat_xct_cytotalk.png', res = 300, height = 1500, width = 2000)
ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")
dev.off()


cytotalk_only <- set7[!set7 %in% intersect(set1, set7)]
xct_only <- set1[!set1 %in% intersect(set1, set7)]

# output pairs
zz <- file("xct_only.txt", "wb")
writeBin(paste(xct_only, collapse="\n"), zz) 
close(zz)



