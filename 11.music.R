
#----- library -----#

library(MuSiC)
library(Biobase)

#single-cell RNA-seq
#load data
sce.mat <- read.table("./sce.counts.txt", row.names = 1, header = TRUE, check.names = FALSE)
sce.clinical <- read.table("./sce.metadata.txt", row.names = 1, header = TRUE, sep = "\t", 
                           check.names = FALSE)
#sce.clinical <- sce.clinical[which(sce.clinical[,6] == "Tumor"),]

sce.clinical <- sce.clinical[,c(6, 11)]
colnames(sce.clinical) <- c("sampleID", "celltype")
sce.clinical <- sce.clinical[rownames(sce.clinical) %in% colnames(sce.mat),]

sce.mat <- sce.mat[,rownames(sce.clinical)]
sce.meta <- data.frame(labelDescription = c("sampleID", "celltype"), 
                       row.names = c("sampleID", "celltype"))
sce.eset <- ExpressionSet(assayData = as.matrix(sce.mat), 
                          phenoData =  new("AnnotatedDataFrame", data = sce.clinical, 
                                           varMetadata = sce.meta))


#bulk RNA-seq
#load data
RNA.mat.count <- read.table("RNA.mat.count.txt", header = TRUE, row.names = 1)
RNA.mat.fpkm <- read.table("RNA.mat.fpkm.txt", header = TRUE, row.names = 1)

clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, row.names = 1, sep = "\t")
IDgene.rna <- read.table("IDgene.rna.txt", row.names = 1, header = TRUE)

remove.f <- apply(RNA.mat.fpkm, 1, function(x) length(which(x < 0.5)) > length(x)*0.5)
RNA.mat.fpkm <- RNA.mat.fpkm[-remove.f, ]
RNA.mat.count <- RNA.mat.count[-remove.f, ]

RNA.mat.count <- aggregate(x = RNA.mat.count, by = list(IDgene.rna[rownames(RNA.mat.fpkm_UQ),2]), 
                           FUN = mean)
rownames(RNA.mat.count) <- RNA.mat.count$Group.1
RNA.mat.count <- RNA.mat.count[,-1]

RNA.mat.count <- RNA.mat.count[which(rownames(RNA.mat.count) %in% 
                                           IDgene.rna$V2[which(IDgene.rna$V3 == "protein_coding")]),]


#MuSiC analysis
bulk.meta <- data.frame(labelDescription = c("sampleID","subjectID"), 
                        row.names = c("sampleID","subjectID"))

bulk.clinical <- data.frame(sampleID = colnames(RNA.mat.count), 
                            subjectID = clinicalData[colnames(RNA.mat.count), "WHO_2021_subtype"])
rownames(bulk.clinical) <- colnames(RNA.mat.count)

bulk.eset <- ExpressionSet(assayData = as.matrix(RNA.mat.count), 
                           phenoData =  new("AnnotatedDataFrame", data = bulk.clinical, varMetadata = bulk.meta))

music.res <- music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'celltype', verbose = F)
music.res <- music.res$Est.prop.allgene

