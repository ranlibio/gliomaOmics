
#----- library -----#

library(estimate)
library(limma)
library(ConsensusClusterPlus)
library(ggplot2)
library(estimate)
library(ggnewscale)


#----- mRNA matrix -----#

#load data
RNA.mat.count <- read.table("RNA.mat.count.txt", header = TRUE, row.names = 1)
RNA.mat.fpkm <- read.table("RNA.mat.fpkm.txt", header = TRUE, row.names = 1)

clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, 
                           row.names = 1, sep = "\t")

IDgene.rna <- read.table("IDgene.rna.txt", row.names = 1, header = TRUE)

remove.f <- apply(RNA.mat.fpkm, 1, function(x) length(which(x < 0.5)) > length(x)*0.5)

RNA.mat.fpkm <- RNA.mat.fpkm[-remove.f, ]
RNA.mat.count <- RNA.mat.count[-remove.f, ]

#fpkm-UQ
RNA.mat.fpkm_UQ <- RNA.mat.fpkm

i <- 1
for(i in 1:ncol(RNA.mat.fpkm_UQ))
{
  RNA.mat.fpkm_UQ[,i] <- RNA.mat.fpkm[,i]*(sum(RNA.mat.countt[,i])/quantile(data.count[,i], 0.75))
}

RNA.mat.fpkm_UQ <- aggregate(x = RNA.mat.fpkm_UQ, 
                             by = list(IDgene.rna[rownames(RNA.mat.fpkm_UQ), 2]), FUN = mean)
rownames(RNA.mat.fpkm_UQ) <- RNA.mat.fpkm_UQ$Group.1
RNA.mat.fpkm_UQ <- as.matrix(RNA.mat.fpkm_UQ[,-1])

#RNA.mat.fpkm_UQ <- log2(RNA.mat.fpkm_UQ + 1)
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[which(rownames(RNA.mat.fpkm_UQ) %in% 
                                           IDgene.rna$V2[which(IDgene.rna$V3 == "protein_coding")]),]



#ssGSEA analysis

microglia.markers <- c("P2RY12", "TMEM119", "SLC2A5", "TGFBR1", "GPR34", "SALL1", "GAS6", "MERTK", 
                       "C1QA", "PROS1", "CD68", "ADGRE1", "AIF1", "CX3CR1", "TREM2", "ITGAM")

go.mt <- list()
go.mt$microglia <- microglia.markers

GSVA.rna.mat <- as.matrix(log2(RNA.mat.fpkm_UQ + 1))
ssgsea.rna.res <- gsva(GSVA.rna.mat, go.mt, method = "ssgsea")


#xCell
xCell.p <- read.table("./xCell.pvals.txt", row.names = 1, header = TRUE, sep = "\t")
xCell.score <- read.table("./xCell.rawScore.txt", row.names = 1, header = TRUE, sep = "\t")

xCell.p <- xCell.p[1:64, which(colnames(xCell.score) %in% rownames(clinicalData))]
xCell.score <- xCell.score[1:64, which(colnames(xCell.score) %in% rownames(clinicalData))]

sign.celltype <- c()

i <- 1
for (i in 1:nrow(xCell.p)) {
  
  tmp <- which(xCell.p[i,] < 0.01)
  sign.celltype <- c(sign.celltype, length(tmp))
  
}

xCell.score <- xCell.score[which(sign.celltype > 6),]
xCell.score <- rbind(xCell.score, ssgsea.rna.res[,colnames(xCell.score)])
rownames(xCell.score)[nrow(xCell.score)] <- "microglia"


#Consensus clustering
cc.res <- ConsensusClusterPlus(d = as.matrix(xCell.score), maxK = 6, reps = 1000, 
                               plot = "pdf", clusterAlg = "km", distance = "euclidean", 
                               title = "Immune subtypes")

#sil.sim <- silhouette_SimilarityMatrix(cc.res[[4]]$consensusClass, cc.res[[4]]$consensusMatrix)
#plot(sil.sim)


#estimate

estimateCount <- RNA.mat.count[which(rownames(RNA.mat.count) %in% 
                                       IDgene.rna$V2[which(IDgene.rna$V3 == "protein_coding")]),]

estimateCount <- aggregate(estimateCount, by = list(IDgene.rna[rownames(estimateCount), 2]), 
                           FUN = mean)
rownames(estimateCount) <- estimateCount$Group.1
estimateCount <- estimateCount[,-1]

estimate.cpm <- as.data.frame(log2(edgeR::cpm(estimateCount) + 0.01))

#write.table(estimate.cpm, file = "estimate.cpm.txt", row.names = TRUE, 
#            col.names = TRUE, sep = "\t", quote = FALSE)

filterCommonGenes(input.f = "estimate.cpm.txt", output.f = "estimate.immune.gct", 
                  id = "GeneSymbol")

estimateScore(input.ds = "estimate.immune.gct", output.ds = "estimateScore.gct", 
              platform = "illumina")

estimateScore  <- read.table("estimateScore.gct", skip = 2, header = T)
rownames(estimateScore) <- estimateScore[,1]
estimateScore <- estimateScore[,3:ncol(estimateScore)]


#Contour plot

im1 <- names(cc.res[[4]]$consensusClass)[which(cc.res[[4]]$consensusClass == 1)]
im2 <- names(cc.res[[4]]$consensusClass)[which(cc.res[[4]]$consensusClass == 2)]
im3 <- names(cc.res[[4]]$consensusClass)[which(cc.res[[4]]$consensusClass == 3)]
im4 <- names(cc.res[[4]]$consensusClass)[which(cc.res[[4]]$consensusClass == 4)]

plotData <- data.frame(x = scale(as.numeric(estimateScore[2, c(im1, im2, im3, im4)])), 
                       y = scale(as.numeric(estimateScore[1, c(im1, im2, im3, im4)])), 
                       group = rep(c("im1", "im2", "im3", "im4"), 
                                   time = c(length(im1), length(im2), length(im3), length(im4))))

rownames(plotData) <- c(im1, im2, im3, im4)

plotData.1 <- plotData[im1,]
plotData.2 <- plotData[im2,]
plotData.3 <- plotData[im3,]
plotData.4 <- plotData[im4,]

#plot
p <- ggplot() + scale_x_continuous(limits = c(-2.5, 3)) + scale_y_continuous(limits = c(-2.5, 3)) + 
  ylab("StromalScore")+xlab("ImmuneScore") + theme_bw() + theme(panel.grid = element_blank())

p1 <- p + stat_density_2d(data = plotData.1, aes(x, y, fill = ..level..), geom = "polygon", 
                          colour = NA, alpha = 0.4, show.legend = FALSE) + 
  scale_fill_gradient(low = "#ed5e5f", high = "#E41A1C")

p2 <- p1 + new_scale("fill") + stat_density_2d(data = plotData.2, aes(x, y, fill = ..level..), 
                                               geom = "polygon", colour = NA, alpha = 0.4, 
                                               show.legend = FALSE) + 
  scale_fill_gradient(low = "#69a3d2", high = "#377EB8")

p3 <- p2 + new_scale("fill") + stat_density_2d(data = plotData.3, aes(x, y, fill = ..level..), 
                                               geom = "polygon", colour = NA, alpha = 0.4, 
                                               show.legend = FALSE) + 
  scale_fill_gradient(low = "#80c87d", high = "#4DAF4A")

p3 + new_scale("fill") + stat_density_2d(data = plotData.4, aes(x, y, fill = ..level..), 
                                         geom = "polygon", 
                                         colour = NA, alpha = 0.1, show.legend = FALSE) + 
  scale_fill_gradient(low = "#b87dc1", high = "#984EA3")


