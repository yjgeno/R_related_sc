library(fgsea)
library(GSVA)
library(data.table)
library(ggplot2)

rm(list = ls())

### example
data(examplePathways)
data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  eps = 0.0
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")


# GSEA
res <- read.table("./gene_list_KL_GSEA.txt") # rnk-like table
BIOP <- gmtPathways("./GO_Biological_Process_2021.txt")

set.seed(123)
gene_list <- as.numeric(res$V2)
names(gene_list) <- toupper(as.vector(res$V1))
gene_list <- sort(gene_list, decreasing = TRUE)

U <- gene_list
U <- U[!grepl('RPL|RPS|RP[[:digit:]]+|MT-', names(U))]
#U <- U[!grepl('GM|RIK', names(U))]
eD <- fgseaMultilevel(BIOP, U, scoreType = 'pos')
eD$leadingEdge <- unlist(lapply(eD$leadingEdge, function(X){paste0(X,collapse = ';')}))
eD <- eD[eD$padj < 0.05,,drop=FALSE]
# write.csv(eD, file = paste0('GSEA_BIOP.csv')) 


# Plot the enrichment plot of the specific pathway
png('GSEA_cholesterol_efflux.png', width = 1000, height = 800, res = 260)
gSet <- 'positive regulation of cholesterol efflux (GO:0010875)'
pTitle <- 'positive regulation of cholesterol efflux'
plotEnrichment(BIOP[[gSet]], U, ticksSize = 0.001) + 
  labs(title = pTitle) +
  theme_bw() + 
  ylab('Enrichment score') +
  xlab('Gene rank') +
  theme(plot.title = element_text(face = 2, size = 12))
dev.off()


