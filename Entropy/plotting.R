rm(list = ls())
library(SCENT)
library(ggplot2)

setwd("data/")
ccat_p21 <- read.table(file = "ccat_p21_G1.txt", 
                   header = F,
                   sep = "\t",
                   as.is = TRUE)$V1
ccat_p23 <- read.table(file = "ccat_p23_G1.txt", 
                       header = F,
                       sep = "\t",
                       as.is = TRUE)$V1
celltype_p21 <- read.table(file = "celltype_p21_G1.txt", 
                       header = F,
                       sep = "\t",
                       as.is = TRUE)$V1
celltype_p23 <- read.table(file = "celltype_p23_G1.txt", 
                       header = F,
                       sep = "\t",
                       as.is = TRUE)$V1

ccat_df <- data.frame(
  ccat = c(ccat_p21, ccat_p23),
  cell.type = c(celltype_p21, celltype_p23),
  batch = c(rep("p21", length(celltype_p21)), rep("p23", length(celltype_p23))),
  stringsAsFactors = T
)

png('ccat_G1.png', width = 3500, height = 2000, res = 200)
ggplot(ccat_df, aes(cell.type, ccat, fill=factor(batch))) +
  geom_boxplot() +
  theme_bw() + 
  theme(text = element_text(size = 20)) +
  xlab('Cell type') + 
  ylab('CCAT') + 
  guides(fill=guide_legend(title="Batch")) +
  theme(plot.title = element_text(face = 2))
dev.off()

