
#----- library -----#

library(limma)
library(survival)
library(survminer)
library(survMisc)
library(ConsensusClusterPlus)
library(CancerSubtypes)


#load data

protein.mat <- read.table("protein.mat.txt", header = TRUE, row.names = 1)
IDgene.p <- read.table("IDgene.p.txt", row.names = 1, header = TRUE)
clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, 
                           row.names = 1, sep = "\t")

IDHwt <- intersect(rownames(clinicalData)[which(clinicalData$WHO_2021_subtype == "IDHwt")], 
                   colnames(protein.mat))


#Consensus clustering

#pre-processing
retain.cols <- apply(protein.mat, 2, function(x) length(which(x >= 100)) > length(x)*0.5)
protein.mat <- protein.mat[,retain.cols]
retain.rows <- apply(protein.mat, 1, function(x) length(which(x >= 100)) > length(x)*0.5)
protein.mat <- protein.mat[retain.rows,]

protein.mat <- normalizeQuantiles(protein.mat)
protein.mat <- log2(protein.mat + 0.01)
protein.mat <- protein.mat[,intersect(colnames(protein.mat), IDHwt)]

mad.res <- apply(protein.mat, 1, mad)
protein.mat <- protein.mat[order(mad.res, decreasing = TRUE)[1:913],]


#clustering
cc.result <- ConsensusClusterPlus(d = as.matrix(protein.mat), maxK = 6, pItem = 0.8, pFeature = 0.8, 
                                  reps = 1000, clusterAlg = "hc", distance = "spearman", 
                                  title = "IDHwt subtypes", plot = "pdf", seed = 6)

#sil <- silhouette_SimilarityMatrix(cc.result[[3]]$consensusClass, 
#                                   cc.result[[3]]$consensusMatrix)
#plot(sil)


#survival analysis

Sur <- data.frame(clinicalData[colnames(protein.mat), "OS"], 
                  clinicalData[colnames(protein.mat), "Death"])

colnames(Sur) <- c("time","status")
rownames(Sur) <- colnames(protein.mat)

Sur$label <- ""
Sur[names(cc.result[[3]]$consensusClass), "label"] <- cc.result[[3]]$consensusClass

Sur <- Sur[which(!is.na(Sur$time) & !is.na(Sur$status)),]

Sur.fit <- survfit(Surv(time,status) ~ label, data = Sur)

ggsurv <- ggsurvplot(Sur.fit,pval = TRUE,risk.table = TRUE,pval.method = TRUE)

ggsurv$plot+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  scale_x_continuous(limits = c(0,35),breaks = seq(0,35,5))


