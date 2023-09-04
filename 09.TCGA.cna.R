
#----- library -----#

library(maftools)
library(stringr)
library(ggplot2)
library(RColorBrewer)

#----- load data -----#

load("mutLGG.rda")
load("mutGBM.rda")

##LGG and GBM mutation

gliomaMut <- rbind(LGGmut, GBMmut)
glioma.maf <- read.maf(maf = gliomaMut, isTCGA = TRUE)

tcgaIDHmut <- unique(glioma.maf@data$Tumor_Sample_Barcode[which(glioma.maf@data$Hugo_Symbol == "IDH1" | 
                                                                  glioma.maf@data$Hugo_Symbol == "IDH2")])
tcgaIDHwt <- setdiff(unique(glioma.maf@data$Tumor_Sample_Barcode), tcgaIDHmut)

tcgaIDHmut <- unique(glioma.maf@data$Tumor_Sample_Barcode[which(glioma.maf@data$Hugo_Symbol == "IDH1" | 
                                                                  glioma.maf@data$Hugo_Symbol == "IDH2")])
tcgaIDHwt <- setdiff(unique(glioma.maf@data$Tumor_Sample_Barcode), tcgaIDHmut)

tcgaIDHmut <- sapply(tcgaIDHmut, function(x) paste0(unlist(str_split(x, "-"))[1:4], collapse = "-"))
tcgaIDHmut <- sapply(tcgaIDHmut, function(x) substr(x, 1, nchar(x)-1))
tcgaIDHmut <- unique(tcgaIDHmut)

tcgaIDHwt <- sapply(tcgaIDHwt, function(x) paste0(unlist(str_split(x, "-"))[1:4], collapse = "-"))
tcgaIDHwt <- sapply(tcgaIDHwt, function(x) substr(x, 1, nchar(x)-1))
tcgaIDHwt <- unique(tcgaIDHwt)


#LGG and GBM CNA

lgg.cna <- read.delim("./lgg_tcga/data_linear_CNA.txt", header = TRUE, check.names = FALSE)
lgg.cna <- lgg.cna[which(duplicated(lgg.cna$Hugo_Symbol) == FALSE),]
rownames(lgg.cna) <- lgg.cna$Hugo_Symbol
lgg.cna <- lgg.cna[,-c(1,2)]

gbm.cna <- read.delim("./gbm_tcga/data_linear_CNA.txt", header = TRUE, check.names = FALSE)
gbm.cna <- gbm.cna[which(duplicated(gbm.cna$Hugo_Symbol) == FALSE),]
rownames(gbm.cna) <- gbm.cna$Hugo_Symbol
gbm.cna <- gbm.cna[,-c(1,2)]

glioma.cna <- cbind(lgg.cna, gbm.cna[rownames(lgg.cna),])


#LGG and CNA mRNA expression

lgg.rna <- read.delim("./TCGA-LGG.htseq_fpkm.tsv", header = TRUE, check.names = FALSE, row.names = 1)
gbm.rna <- read.delim("./TCGA-GBM.htseq_fpkm.tsv", header = TRUE, check.names = FALSE, row.names = 1)
glioma.rna <- cbind(lgg.rna, gbm.rna[rownames(lgg.rna),])

annot <- read.table("./gencode.v22.annotation.gene.probeMap", header = TRUE, row.names = 1)

glioma.rna <- aggregate(x = glioma.rna, by = list(annot[rownames(glioma.rna), "gene"]), FUN = mean)
rownames(glioma.rna) <- glioma.rna$Group.1
glioma.rna <- as.matrix(glioma.rna[,-1])

colnames(glioma.rna) <- sapply(colnames(glioma.rna), function(x) substr(x, 1, nchar(x)-1))

remove.f <- apply(glioma.rna, 1, function(x) length(which(x > 0.5)) > length(x)/2)
glioma.rna <- glioma.rna[remove.f,]


#overlap
overlap.Gene <- intersect(rownames(glioma.cna), rownames(glioma.rna))

tcgaIDHmut <- intersect(intersect(tcgaIDHmut, colnames(glioma.cna)), colnames(glioma.rna))
tcgaIDHwt <- intersect(intersect(tcgaIDHwt, colnames(glioma.cna)), colnames(glioma.rna))

tcgaIDHmut.cna <- glioma.cna[overlap.Gene, tcgaIDHmut]
tcgaIDHmut.rna <- glioma.rna[overlap.Gene, tcgaIDHmut]

tcgaIDHwt.cna <- glioma.cna[overlap.Gene, tcgaIDHwt]
tcgaIDHwt.rna <- glioma.rna[overlap.Gene, tcgaIDHwt]


#correlation analysis
cor.matrix <- matrix(0, nrow(tcgaIDHmut.cna), 2)

i <- 1
for (i in 1:nrow(cor.matrix)){
  cor.matrix[i,1] <- cor(as.numeric(tcgaIDHmut.cna[i,]), as.numeric(tcgaIDHmut.rna[i,]), method = "spearman")
  cor.matrix[i,2] <- cor(as.numeric(tcgaIDHwt.cna[i,]), as.numeric(tcgaIDHwt.rna[i,]), method = "spearman")
}


#----- Figure S3B -----#

#plot
plotData <- data.frame(y = c(cor.matrix[,1],cor.matrix[,2]),
                       x = rep(c("tcga_IDHmut","tcga_IDHwt"), each = nrow(cor.matrix)))

ggplot(plotData, aes(x = x, y = y, fill = x)) + 
  geom_boxplot(outlier.colour = NA, width = 0.3) + 
  scale_y_continuous(limits = c(-0.4,0.8), breaks = seq(-0.4,0.8,0.2)) + 
  scale_fill_manual(values = brewer.pal(3, "Set1")) + 
  ylab("Spearman's correlation") + xlab("") + 
  guides(fill = FALSE) + 
  theme_bw() + theme(panel.grid = element_blank())


#----- Figure S3F -----#

#plot
ratio.1 <- apply(tcgaIDHmut.cna, 2, function(x) length(which(x > (0.2)))/length(x)) # x > 0.2 or x < (-0.2)
ratio.2 <- apply(tcgaIDHwt.cna, 2, function(x) length(which(x > (0.2)))/length(x))

plotData <- data.frame(x = rep(c("tcgaIDHmut", "tcgaIDHwt"), time = c(ncol(tcgaIDHmut.cna), ncol(tcgaIDHwt.cna))), 
                       y = c(ratio.1, ratio.2)*100)

ggplot(plotData, aes(x = x, y = y, fill = x)) + 
  geom_boxplot(outlier.colour = NA, width = 0.3) + 
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) + 
  scale_fill_manual(values = brewer.pal(3,"Set1")) + 
  ylab("Percentage of log2(cn/2)<-0.2") + xlab("") + 
  guides(fill = FALSE) + 
  theme_bw() + theme(panel.grid = element_blank())

