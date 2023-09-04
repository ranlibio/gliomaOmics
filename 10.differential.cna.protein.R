
#----- library -----#

library(stringr)
library(limma)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(survival)
library(survminer)
library(survMisc)


#----- CNA matrix -----#

#load data

#mutation matrix
cna.mat <- read.table("cna.mat.txt", row.names = 1, header = TRUE, sep = "\t")


#----- mRNA matrix -----#

#load data
RNA.mat.count <- read.table("RNA.mat.count.txt", header = TRUE, row.names = 1)
RNA.mat.fpkm <- read.table("RNA.mat.fpkm.txt", header = TRUE, row.names = 1)

IDgene.rna <- read.table("IDgene.rna.txt", row.names = 1, header = TRUE)

remove.f <- apply(RNA.mat.fpkm, 1, function(x) length(which(x < 0.5)) > length(x)*0.8)

RNA.mat.fpkm <- RNA.mat.fpkm[-remove.f, ]
RNA.mat.count <- RNA.mat.count[-remove.f, ]

#fpkm-UQ
RNA.mat.fpkm_UQ <- RNA.mat.fpkm

i <- 1
for(i in 1:ncol(RNA.mat.fpkm_UQ))
{
  RNA.mat.fpkm_UQ[,i] <- RNA.mat.fpkm[,i]*(sum(RNA.mat.countt[,i])/quantile(data.count[,i], 0.75))
}

RNA.mat.fpkm_UQ <- aggregate(x = RNA.mat.fpkm_UQ, by = list(IDgene.rna[rownames(RNA.mat.fpkm_UQ),2]), FUN = mean)
rownames(RNA.mat.fpkm_UQ) <- RNA.mat.fpkm_UQ$Group.1
RNA.mat.fpkm_UQ <- as.matrix(RNA.mat.fpkm_UQ[,-1])

RNA.mat.fpkm_UQ <- log2(RNA.mat.fpkm_UQ + 0.01)
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[which(rownames(RNA.mat.fpkm_UQ) %in% 
                                           IDgene.rna$V2[which(IDgene.rna$V3 == "protein_coding")]),]


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
protein.mat <- protein.mat[,which(cols.f > (nrow(protein.mat)*0.4))]

j <- 1
rows.f <- c()
for (j in 1:nrow(protein.mat)) {
  len <- length(which(protein.mat[j,] >= 100))
  rows.f <- c(rows.f, len)
}

protein.mat <- protein.mat[which(rows.f > (ncol(protein.mat)*0.4)),]
protein.mat <- normalizeQuantiles(protein.mat)

IDgene.p <- read.table("IDgene.p.txt", row.names = 1, header = TRUE)

protein.mat <- aggregate(x = protein.mat, by = list(IDgene.p[rownames(protein.mat),1]), FUN = mean)
rownames(protein.mat) <- protein.mat$Group.1
protein.mat <- as.matrix(protein.mat[,-1])
protein.mat <- log2(protein.mat + 1)


#----- Correlation analysis -----#

overlap.Gene <- intersect(intersect(rownames(cna.mat), rownames(protein.mat)), 
                          rownames(RNA.mat.fpkm_UQ))
overlap.Sample <- intersect(intersect(colnames(cna.mat), colnames(protein.mat)), 
                            colnames(RNA.mat.fpkm_UQ))

cna.mat <- cna.mat[overlap.Gene, overlap.Sample]
protein.mat <- protein.mat[overlap.Gene, overlap.Sample]
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[overlap.Gene, overlap.Sample]

#CNA-mRNA
cna_rna_cor <- matrix(0, nrow(cna.mat), 2)
rownames(cna_rna_cor) <- rownames(cna.mat)
i <- 1
for (i in 1:nrow(cna_rna_cor)) {
  
  cor.value <- cor.test(as.numeric(cna.mat[i,]), as.numeric(RNA.mat.fpkm_UQ[i,]), method = "spearman")
  cna_rna_cor[i,1] <- cor.value$estimate
  cna_rna_cor[i,2] <- cor.value$p.value
}
cna_rna_cor[,2] <- p.adjust(cna_rna_cor[,2], method = "BH")

#CNA-protein
cna_protein_cor <- matrix(0, nrow(cna.mat), 2)
rownames(cna_protein_cor) <- rownames(cna.mat)
i <- 1
for (i in 1:nrow(cna_protein_cor)) {
  
  cor.value <- cor.test(as.numeric(cna.mat[i,]), as.numeric(protein.mat[i,]), method = "spearman")
  cna_protein_cor[i,1] <- cor.value$estimate
  cna_protein_cor[i,2] <- cor.value$p.value
}
cna_protein_cor[,2] <- p.adjust(cna_protein_cor[,2], method = "BH")


#----- ANOVA analysis -----#

clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, row.names = 1, sep = "\t")

#CNA
cna.aov <- matrix(0, nrow(cna.mat), 2)
rownames(cna.aov) <- rownames(cna.mat)
i <- 1
for (i in 1:nrow(cna.aov)) {
  
  tmp <- data.frame(cna = as.numeric(cna.mat[rownames(cna.aov)[i],]), 
                    class = clinicalData[colnames(cna.mat), "WHO_2021_subtype"])
  aov.res <- aov(cna~class, data = tmp)
  aov.res <- summary(aov.res)
  cna.aov [i,1] <- aov.res[[1]][1,5]
}
cna.aov[,1] <- p.adjust(cna.aov[,1], method = "BH")

#Protein
i <- 1
for (i in 1:nrow(cna.aov)) {
  
  tmp <- data.frame(cna = as.numeric(protein.mat[rownames(cna.aov)[i],]), 
                    class = clinicalData[colnames(protein.mat), "WHO_2021_subtype"])
  aov.res <- aov(cnv~class, data = tmp)
  aov.res <- summary(aov.res)
  cna.aov[i,2] <- aov.res[[1]][1,5]
}
cna.aov[,2] <- p.adjust(cna.aov[,2],method = "BH")


#----- Figure 2F -----#

#Overlap
sign.gene <- rownames(cna.aov)[which(cna.aov[,1] < 0.01 & cna.aov[,2] < 0.01)]


#----- Chromosome analysis -----#

annot.chr <- read.table("./gencode.v22.annot.chration.gene.probeMap", header = TRUE, row.names = 1)

chr7.gene <- annot.chr$gene[annot.chr$chrom == "chr7"] #chr1, chr7, chr10, chr19
chr7.gene <- intersect(sign.gene,chr7.gene)

IDHmut_codel <- intersect(rownames(clinicalData)[which(clinicalData$WHO_2021_subtype == 'IDHmut-codel')], 
                          colnames(protein.mat))
IDHmut_noncodel <- intersect(rownames(clinicalData)[which(clinicalData$sampleGroup == 'IDHmut-noncodel')], 
                             colnames(protein.mat))
IDHwt <- intersect(rownames(clinicalData)[which(clinicalData$sampleGroup == 'IDHwt')], colnames(protein.mat))


#----- Figure S3G -----#

#Plot
sub1.cna <- rowMeans(cna.mat[chr7.gene, IDHmut_codel])
sub2.cna <- rowMeans(cna.mat[chr7.gene, IDHmut_noncodel])
sub3.cna <- rowMeans(cna.mat[chr7.gene, IDHwt])

sub1.protein <- rowMeans(protein.mat[chr7.gene, IDHmut_codel])
sub2.protein <- rowMeans(protein.mat[chr7.gene, IDHmut_noncodel])
sub3.protein <- rowMeans(protein.mat[chr7.gene, IDHwt])

plotData <- data.frame(x = c(sub1.cna, sub2.cna, sub3.cna), y = c(sub1.protein, sub2.protein, sub3.protein), 
                       group = c(rep(c("IDHmut-codel","IDHmut-noncodel","IDHwt"), each = length(chr7.gene))))

ggplot(plotData, aes(x, y, fill = group)) + geom_point(size = 0.8) + 
  scale_fill_manual(values = rev(brewer.pal(3, "Set1"))) + 
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1,1))+
  scale_y_continuous(breaks = seq(15,27,2), limits = c(15,27.5))+
  xlab("CNAs") + ylab("Protein abundance") + 
  theme_bw() + theme(panel.grid = element_blank())


#----- Figure 2F -----#

#22 genes
cna_rna_protein_cor <- matrix(0, nrow(cna.mat), 6)
rownames(cna_rna_protein_cor) <- rownames(cna.mat)

i <- 1
for (i in 1:nrow(cna_rna_protein_cor)) {
  
  cor.value <- cor.test(as.numeric(cna.mat[rownames(cna.mat)[i],]), as.numeric(RNA.mat.fpkm_UQ[rownames(cna.mat)[i],]), 
                        method = "spearman")
  cna_rna_protein_cor[i,1] <- cor.value$p.value
  cna_rna_protein_cor[i,4] <- cor.value$estimate
  
  cor.value <- cor.test(as.numeric(cna.mat[rownames(cna.mat)[i],]), as.numeric(protein.mat[rownames(cna.mat)[i],]), 
                        method = "spearman")
  cna_rna_protein_cor[i,2] <- cor.value$p.value
  cna_rna_protein_cor[i,5] <- cor.value$estimate
  
  cor.value <- cor.test(as.numeric(RNA.mat.fpkm_UQ[rownames(cna.mat)[i],]), as.numeric(protein.mat[rownames(cna.mat)[i],]), 
                        method = "spearman")
  cna_rna_protein_cor[i,3] <- cor.value$p.value
  cna_rna_protein_cor[i,6] <- cor.value$estimate
}

cna_rna_protein_cor[,1] <- p.adjust(cna_rna_protein_cor[,1], method = "BH")
cna_rna_protein_cor[,2] <- p.adjust(cna_rna_protein_cor[,2], method = "BH")
cna_rna_protein_cor[,3] <- p.adjust(cna_rna_protein_cor[,3], method = "BH")

cna_rna_protein_cor <- cna_rna_protein_cor[sign.gene,]

sign.gene.2 <- rownames(cna_rna_protein_cor)[which(cna_rna_protein_cor[,1] < 0.0001 & 
                                                     cna_rna_protein_cor[,2] < 0.0001 & 
                                                     cna_rna_protein_cor[,3] < 0.01)]

#----- Figure 2I -----#

#Plot
cna.1 <- cna.mat[sign.gene.2, IDHmut_codel]
cna.2 <- cna.mat[sign.gene.2, IDHmut_noncodel]
cna.3 <- cna.mat[sign.gene.2, IDHwt]

plotData <- data.frame(x = c(rep(sign.gene.2, each = 29), rep(sign.gene.2, each = 24), rep(sign.gene.2, each = 50)), 
                       y = c(c(t(cna.1)), c(t(cna.2)), c(t(cna.3))), 
                       group = rep(c("IDHmut-codel","IDHmut-noncodel","IDHwt"), time = c(22*29, 22*24, 22*50)))

plotData$x <- fct_inorder(plotData$x)

ggplot(plotData, aes(x = x, y = y, color = group)) + geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(breaks = seq(-1.5,2,0.5), limits = c(-1.5,2))+
  scale_color_manual(values = rev(brewer.pal(3, "Set1"))) + 
  ylab("") + xlab("") + guides(color=FALSE) + 
  theme_bw() + theme(panel.grid = element_blank())


#Survival analysis for 22 genes
sur.p <- c()

i <- 1
for (i in 1:length(sign.gene.2)) {
  
  Sur <- data.frame(clinicalData[colnames(cna.mat), c("OS", "Death")])
  colnames(Sur) <- c("time", "status")
  Sur <- Sur[which(!is.na(Sur$time) & !is.na(Sur$status)),]
  
  Sur$exp <- as.numeric(cna.mat[sign.gene.2[i], rownames(Sur)])
  cut.off <- surv_cutpoint(Sur, time = "time", event = "status", variables = "exp", minprop = 0.3)
  cut.off <- cut.off$cutpoint$cutpoint
  
  Sur$lable <- "H"
  Sur[which(Sur$exp < cut.off), "lable"] <- "L"
  
  Sur.fit <- survfit(Surv(time,status) ~ lable, data = Sur)
  p.value <- surv_pvalue(Sur.fit)
  sur.p <- c(sur.p, p.value$pval)
  
}


plotData <- data.frame(x = c(sign.gene.2), y = -log10(sur.p), color = -log10(sur.p))
plotData$x <- fct_inorder(plotData$x)

ggplot(plotData,aes(x = x, y = y, color = color)) + geom_point(cex = 2) + 
  scale_colour_gradient(low = "#FEE0D2", high = "#DE2D26") + 
  scale_y_continuous(breaks = seq(1, 7, 1),limits = c(1, 7)) + 
  ylab("-log10(P value)")+xlab("")+
  theme_bw() + theme(panel.grid = element_blank())


#----- Figure 2F, 2G and S3I-----#

#cancer-associated genes
cancerGene <- read.delim("./cancerGeneList.tsv", row.names = 1, header = TRUE, sep = "\t")
cancerGene <- intersect(sign.gene, rownames(cancerGene))

heatmapData <- rbind(cna.mat[cancerGene, c(IDHmut_codel,IDHmut_noncodel,IDHwt)], 
                     protein.mat[cancerGene, c(IDHmut_codel,IDHmut_noncodel,IDHwt)])

annot_col <- data.frame(
  groups = factor(clinicalData[colnames(heatmapData), "WHO_2021_subtype"])
)
rownames(annot_col) <- colnames(heatmapData)

ann_colors <- list(groups = c(IDHmut_codel = "#E41A1C", IDHmut_noncodel = "#377EB8", IDHwt = "#4DAF4A"))

bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(heatmapData, scale = "row", show_rownames = TRUE, show_colnames = FALSE, 
         cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, 
         fontsize_row = 8, treeheight_row = 0, 
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2), 
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
         legend_breaks = seq(-2, 2, 2), breaks = bk, 
         annotation_col = annot_col, annotation_colors = ann_colors)

