library(Seurat)
library(CytoTalk)
library(reticulate)
use_python("path_to_Miniconda/envs/py/python.exe", required = T)

setwd("path")

seuratObject <- readRDS('humanSkin.rds')
# s <- subset(x = seuratObject, ident = c("Inflam. FIB", "Inflam. DC"))
HVG <- VariableFeatures(seuratObject)
s <- seuratObject[HVG,]

counts <- GetAssayData(s, slot = 'counts')
s_cyto <- list(mat = counts, cell_types = s$ident)

xct_omni <- read.csv('omnipath_intercell_toUse_forLIANA.csv') # xct omnipath DB
lrp_human_xct <- data.frame(ligand = xct_omni$source_genesymbol, 
                      receptor = xct_omni$target_genesymbol)

results <- CytoTalk::run_cytotalk(s_cyto, "Inflam. FIB", "Inflam. DC",
                                  pcg = CytoTalk::pcg_human, # protein-coding genes
                                  lrp = lrp_human_xct, # L-R database
                                  cutoff_a=0.1, cutoff_b=0.1, # gene cutoffs
                                  beta_max=500, # tune beta: a parameter for balancing the edge costs and node prizes
                                  omega_min=0.1, omega_max=1.5, # tune w: cost of all artificial edges (same)
                                  depth=3, ntrial=1000,
                                  cores=NULL, echo=TRUE, 
                                  dir_out = './CytoTalk/')
saveRDS(results, 'CytoTalk/results.rds')
results$pathways


# example
load("C:/Users/yjyang027/Downloads/test/data/scrna_cyto.rda")
scrna_cyto$mat <- scrna_cyto$mat[1:2000, ]
# run CytoTalk process
results <- CytoTalk::run_cytotalk(scrna_cyto, "EndothelialCells", "Macrophages")


