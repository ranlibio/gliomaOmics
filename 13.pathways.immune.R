
#----- library -----#

library(limma)
library(edgeR)
library(clusterProfiler)


#function for differential analysis
DE.fun <- function(s1, s2, protein.mat, IDgene){
  
  #pre-processing 1
  retain.f <- apply(protein.mat, 2, function(x) length(which(x >= 100)) > length(x)*0.4)
  protein.mat <- protein.mat[,retain.f]
  
  s1 <- intersect(s1, colnames(protein.mat))
  s2 <- intersect(s2, colnames(protein.mat))
  protein.mat <- protein.mat[,c(s1,s2)]
  
  #pre-processing 2
  i <- 1
  rows_f <- c()
  for (i in 1:nrow(protein.mat)) {
    
    t1 <- which(protein.mat[i, s1] > 100)
    t2 <- which(protein.mat[i, s2] > 100)
    
    if(length(t1) > length(s1)*0.3 & length(t2) > length(s2)*0.3){
      
      rows_f <- c(rows_f, i)
      
    }
    
  }
  
  protein.mat <- protein.mat[rows_f,]
  protein.mat <- normalizeQuantiles(protein.mat)
  protein.mat <- log2(protein.mat + 1)
  
  #limma
  group_list <- factor(c(rep("n1", length(s1)), rep("n2", length(s2))))
  design <- model.matrix(~0 + group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(protein.mat)
  
  contrast.matrix <- makeContrasts("n1-n2", levels = design)
  fit <- lmFit(protein.mat, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit, trend = TRUE)
  
  de.res <- topTable(fit, n=Inf, adjust.method = "BH")
  de.res$Gene <- IDgene[rownames(de.res), 1]
  return(de.res)
  
}


imCluster <- read.table("./immuneCluster.txt", row.names = 1, header = TRUE)
im1 <- rownames(immuneCluster)[which(immuneCluster$immuneCluster == 1)]
im2 <- rownames(immuneCluster)[which(immuneCluster$immuneCluster == 2)]
im3 <- rownames(immuneCluster)[which(immuneCluster$immuneCluster == 3)]
im4 <- rownames(immuneCluster)[which(immuneCluster$immuneCluster == 4)]


#load data
protein.mat <- read.table("protein.mat.txt", header = TRUE, row.names = 1)
IDgene.p <- read.table("IDgene.p.txt", row.names = 1, header = TRUE)

#differential analysis
#
n1.sample <- im1
n2.sample <- c(im2, im3, im4)
proteinDE_im1 <- DE.fun(n1.sample, n2.sample, protein.mat, IDgene.p)

#
n1.sample <- im2
n2.sample <- c(im1, im3, im4)
proteinDE_im2 <- DE.fun(n1.sample, n2.sample, protein.mat, IDgene.p)

#
n1.sample <- im3
n2.sample <- c(im1, im2, im4)
proteinDE_im3 <- DE.fun(n1.sample, n2.sample, protein.mat, IDgene.p)

#
n1.sample <- im4
n2.sample <- c(im1, im2, im3)
proteinDE_im4 <- DE.fun(n1.sample,n2.sample,protein.mat,IDgene.p)


#GSEA
#proteinDE_im1, proteinDE_im2, proteinDE_im3 and proteinDE_im4
rank_geneSet <- proteinDE_im1$logFC 
names(rank_geneSet) <- proteinDE_im1$Gene
rank_geneSet <- sort(rank_geneSet, decreasing = TRUE)

go.mt <- read.gmt("./c2.cp.v7.4.symbols.gmt")

gsea_res <- GSEA(rank_geneSet, TERM2GENE = go.mt, pvalueCutoff = 0.01, eps = 0, 
                    minGSSize = 10, maxGSSize = 300, seed = TRUE)
gsea_res <- gsea_res@result

#write.table(gsea_res, file = "./gsea_res_im1.txt", row.names = TRUE, 
#            col.names = TRUE, sep = "\t", quote = FALSE)






