library(fgsea)
library(GSVA)
library(data.table)
library(ggplot2)

rm(list = ls())

### example
data(examplePathways)
data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, # List of gene sets to check
                  stats    = exampleRanks, # Named vector of gene-level stats. Names should be the same as in 'pathways'
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
eD <- fgseaMultilevel(BIOP, U, scoreType = 'pos') # scoreType: https://github.com/ctlab/fgsea/issues/27
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

### Appendix: R script for obtain gnk
library(GEOquery)
library(limma)
library(org.Mm.eg.db)
library(data.table)
# for collapseBy
source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")

gse14308 <- getGEO("GSE14308")[[1]]

pData(gse14308)$condition <- sub("-.*$", "", gse14308$title)

es <- collapseBy(gse14308, fData(gse14308)$ENTREZ_GENE_ID, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)+1), method="quantile")

es.design <- model.matrix(~0+condition, data=pData(es))

fit <- lmFit(es, es.design)

fit2 <- contrasts.fit(fit, makeContrasts(conditionTh1-conditionNaive,
                                         levels=es.design))
fit2 <- eBayes(fit2)
de <- data.table(topTable(fit2, adjust.method="BH", number=12000, sort.by = "B"), keep.rownames = TRUE)


ranks <- de[order(t), list(rn, t)]
write.table(ranks, "inst/extdata/naive.vs.th1.rnk", sep="\t", quote = FALSE, row.names = FALSE)

