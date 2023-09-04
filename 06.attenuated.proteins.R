
#----- library -----#

library(mclust)
library(limma)
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(clusterProfiler)


#----- Figure 2A and 2B -----#

#----- CNA matrix -----#

#load data
cna.mat <- read.table("cna.mat.txt", row.names = 1, header = TRUE, sep = "\t")

#----- Protein matrix -----#

#load data
protein.mat <- read.table("protein.mat.txt", header = TRUE, row.names = 1)

#filter
i <- 1
cols.f <- c()
for (i in 1:ncol(protein.mat)) {
  len <- length(which(protein.mat[,i] >= 100))
  cols.f <- c(cols.f,len)
}
protein.mat <- protein.mat[,which(cols.f > (nrow(protein.mat)*0.6))]

j <- 1
rows.f <- c()
for (j in 1:nrow(protein.mat)) {
  len <- length(which(protein.mat[j,] >= 100))
  rows.f <- c(rows.f, len)
}

protein.mat <- protein.mat[which(rows.f > (ncol(protein.mat)*0.6)),]
protein.mat <- normalizeQuantiles(protein.mat)

IDgene.p <- read.table("IDgene.p.txt", row.names = 1, header = TRUE)

protein.mat <- aggregate(x = protein.mat, by = list(IDgene.p[rownames(protein.mat),1]), FUN = mean)
rownames(protein.mat) <- protein.mat$Group.1
protein.mat <- as.matrix(protein.mat[,-1])
#protein.mat <- log2(protein.mat+1)


#----- mRNA matrix -----#

#load data
RNA.mat.count <- read.table("RNA.mat.count.txt", header = TRUE, row.names = 1)
RNA.mat.fpkm <- read.table("RNA.mat.fpkm.txt", header = TRUE, row.names = 1)

rows.f <- which(rowMeans(RNA.mat.fpkm) < 1)
RNA.mat.fpkm <- RNA.mat.fpkm[-rows.f,]
RNA.mat.count <- RNA.mat.count[-rows.f,]

#fpkm-UQ
RNA.mat.fpkm_UQ <- RNA.mat.fpkm

i <- 1
for(i in 1:ncol(RNA.mat.fpkm_UQ))
{
  RNA.mat.fpkm_UQ[,i] <- RNA.mat.fpkm[,i]*(sum(RNA.mat.countt[,i])/quantile(data.count[,i], 0.75))
}

IDgene.rna <- read.table("IDgene.rna.txt", row.names = 1, header = TRUE)

RNA.mat.fpkm_UQ <- aggregate(x = RNA.mat.fpkm_UQ, by = list(IDgene[rownames(RNA.mat.fpkm_UQ),2]), FUN = mean)
rownames(RNA.mat.fpkm_UQ) <- RNA.mat.fpkm_UQ$Group.1
RNA.mat.fpkm_UQ <- as.matrix(RNA.mat.fpkm_UQ[,-1])
#RNA.mat.fpkm_UQ <- log2(RNA.mat.fpkm_UQ+1)


#----- correlation analysis -----#

overlap.Gene <- intersect(intersect(rownames(cna.mat), rownames(protein.mat)), 
                          rownames(RNA.mat.fpkm_UQ))
overlap.Sample <- intersect(intersect(colnames(cna.mat), colnames(protein.mat)), 
                            colnames(RNA.mat.fpkm_UQ))

cna.mat <- cna.mat[overlap.Gene, overlap.Sample]
protein.mat <- protein.mat[overlap.Gene, overlap.Sample]
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[overlap.Gene, overlap.Sample]


cna_protein_cor <- c()
i <- 1
for (i in 1:nrow(cna.mat)) {
  cna_protein_cor <- c(cna_protein_cor, cor(as.numeric(cna.mat[i,]), 
                                            as.numeric(protein.mat[i,]), method = "spearman"))
}

cna_rna_cor <- c()
i <- 1
for (i in 1:nrow(cna.mat)) {
  cna_rna_cor <- c(cna_ran_cor, cor(as.numeric(cna.mat[i,]), 
                                    as.numeric(RNA.mat.fpkm_UQ[i,]), method = "spearman"))
}


#attenuated proteins

set.seed(10)
mclust.mat <- data.frame(cna_rna_cor, cna_protein_cor)
mod <- Mclust(mclust.mat, G = 2)

cols <- mod$classification
cols[which(cols == 1)] <- "red"
cols[which(cols == 2)] <- "grey"

plotData <- data.frame(x = cna_rna_cor, y = cna_protein_cor, class = as.factor(cols))

p1 <- ggplot(plotData, aes(x = cna_rna_cor, y = cna_protein_cor, color = class)) + 
  geom_point(color = cols, cex=1, alpha = 0.6, pch = 20) + 
  scale_x_continuous(breaks = seq(-0.4, 0.8, 0.2),limits = c(-0.45, 0.8)) + 
  scale_y_continuous(breaks = seq(-0.4, 0.8, 0.2),limits = c(-0.45, 0.8)) + 
  geom_abline(intercept = 0, slope = 1) + stat_density2d() + 
  xlab("CNA-mRNA correlation") + ylab("CNA-protein correlation") + 
  theme_bw() + guides(color = FALSE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())

xdens <- axis_canvas(p1, axis = "x") + geom_density(data = plotData, aes(x = cna_rna_cor, fill = class))

ydens <- axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
  geom_density(data = plotData, aes(x = cna_protein_cor, fill = class)) + coord_flip()

p2 <- insert_xaxis_grob(p1, xdens, grid::unit(.2, "null"), position = "top")
p3 <- insert_yaxis_grob(p2, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p3)


#----- GSEA analysis -----#

ratio <- cna_rna_cor-cna_protein_cor

#ranked gene set
rank.geneSet <- ratio
names(rank.geneSet) <- rownames(cna.mat)
rank.geneSet <- sort(rank.geneSet, decreasing = TRUE)

go.mt <- read.gmt("allComplexes.txt")
go.mt <- gomt[which(nchar(gomt$gene) != 0),]

gsea.res <- GSEA(rank.geneSet, TERM2GENE = go.mt, pvalueCutoff = 0.05, pAdjustMethod = "BH")

#plot GSEA
NES <- gsea.res@result$NES[order(gsea.res@result$NES, decreasing = TRUE)][1:9]
fdr <- gogsea@result$p.adjust[order(gogsea@result$NES, decreasing = TRUE)][1:9]
pathway <- gogsea@result$Description[order(gogsea@result$NES, decreasing = TRUE)][1:9]

plotData <- data.frame(x = pathway, y = fdr)

ggplot(plotData, aes(x, y)) + geom_bar(stat = "identity", width = 0.5)+
  coord_flip() + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=10)) + 
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4))

