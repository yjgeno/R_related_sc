library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

seqdata <- read.delim("GSE146425_gene_counts.txt", header = TRUE, stringsAsFactors = FALSE)
seqdata <- seqdata %>% column_to_rownames('ID')
X <- factor(c('WT', 'WT', 'WT', 'WT', 'KO', 'KO', 'KO'), levels = c('KO', 'WT'))
Categories <- data.frame(X)
rownames(Categories) <- colnames(seqdata)

#ncol(seqdata) == nrow(Categories)
seqDE <- DESeqDataSetFromMatrix(seqdata, Categories, design = ~ X)
seqDE <- DESeq(seqDE)
seqDE <- results(seqDE)
seqDE <- as.data.frame(seqDE) # gene with NA
seqDE <- seqDE[complete.cases(seqDE),] # complete genes list

rownames(seqDE) <- sub("\\..*", "", rownames(seqDE)) # gene names
write.csv(seqDE, file = 'seqDE.csv')

# add gene symbols
gene_Symbols <- read.delim("gene_symbols.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_Symbols) <- 'gene'
seqDE <- cbind(seqDE, new_col = gene_Symbols)

seqDE$gene[!(abs(seqDE$log2FoldChange) > 1 & seqDE$padj < 0.05)] <- NA

seqDE_sign <- seqDE[!(is.na(seqDE$gene)),] # genes that significant DE
write.csv(seqDE_sign, file = 'seqDE_sign.csv')

png('seqDE.png', width = 1500, height = 1500, res = 300)
ggplot(seqDE, aes(log2FoldChange, -log10(padj), label = gene)) + 
  geom_point(color = ifelse(is.na(seqDE$gene), 'black', 'red'), size = 1, alpha = 0.5) + 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "darkgreen", size=0.5) +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", 
             color = "darkgreen", size=0.5) + 
  theme_minimal() + 
  labs(title = 'GSE146425', subtitle = 'WT - KO') +
  xlab(log[2]~(Fold~Change)) + 
  ylab(-log[10]~(adj_P-value)) + 
  geom_text_repel(cex=2.5, segment.size = 0.1) + 
  theme(plot.title = element_text(face = 2, size = 15))
dev.off()  
