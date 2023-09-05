
#----- library -----#

library(TCGAbiolinks)
library(EDASeq)
library(ggplot2)
library(forcats)


#----- CIN -----#

#load data
seg.score <- read.table("./Glioma.Reference.seg", header = TRUE)
seg.score <- seg.score[which(seg.score$chrom != "chrX"),]
seg.score$ID <- sapply(seg.score$ID, function(x) unlist(stringr::str_split(x, "\\."))[1])

seg.score$weightScore <- abs(seg.score$seg.mean)*(seg.score$loc.end - seg.score$loc.start)
sampleID <- unique(seg.score$ID)

sumScores <- c()

i <- 1
for (i in 1:length(sampleID)) {
  sumScores <- c(sumScores, sum(seg.score[which(seg.score$ID %in% sampleID[i]), 7]))
}

sumScores <- sumScores/100000000
names(sumScores) <- sampleID

sampleGroup <- read.table("./sampleGroup.txt", row.names = 1, header = TRUE)

#plot CIN
plotData <- data.frame(x = sampleGroup$sampleGroup, 
                       y = sumScores[rownames(sampleGroup)])
plotData$x <- fct_inorder(plotData$x)

ggplot(plotData, aes(x = x, y = y, fill = x)) + 
  geom_boxplot(outlier.color = "black", outlier.size = 0.8) + 
  scale_fill_manual(values = brewer.pal(5, "Set1")) + 
  scale_y_continuous(breaks = seq(0, 10, 2), limits = c(0, 10)) + 
  xlab("") + ylab("Chromosome instability index") + 
  guides(fill = FALSE) + theme_bw() + theme(panel.grid = element_blank())


#----- stemness -----#

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

RNA.mat.fpkm_UQ <- log2(RNA.mat.fpkm_UQ + 1)
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[which(rownames(RNA.mat.fpkm_UQ) %in% 
                                           IDgene.rna$V2[which(IDgene.rna$V3 == "protein_coding")]),]


#stemness score
stemness.score <- TCGAanalyze_Stemness(stemSig = PCBC_stemSig, dataGE = RNA.mat.fpkm_UQ)
