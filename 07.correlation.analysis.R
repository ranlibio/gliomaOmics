
#----- library -----#

library(limma)
library(forcats)
library(RColorBrewer)
library(ggplot2)

#----- Figure 2C and 2D -----#

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

RNA.mat.fpkm_UQ <- aggregate(x = RNA.mat.fpkm_UQ, by = list(IDgene.rna[rownames(RNA.mat.fpkm_UQ),2]), FUN = mean)
rownames(RNA.mat.fpkm_UQ) <- RNA.mat.fpkm_UQ$Group.1
RNA.mat.fpkm_UQ <- as.matrix(RNA.mat.fpkm_UQ[,-1])
#RNA.mat.fpkm_UQ <- log2(RNA.mat.fpkm_UQ+1)


#----- correlation analysis -----#

clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, row.names = 1, sep = "\t")

overlap.Gene <- intersect(intersect(rownames(cna.mat), rownames(protein.mat)), 
                          rownames(RNA.mat.fpkm_UQ))
overlap.Sample <- intersect(intersect(colnames(cna.mat), colnames(protein.mat)), 
                            colnames(RNA.mat.fpkm_UQ))

cna.mat <- cna.mat[overlap.Gene, overlap.Sample]
protein.mat <- protein.mat[overlap.Gene, overlap.Sample]
RNA.mat.fpkm_UQ <- RNA.mat.fpkm_UQ[overlap.Gene, overlap.Sample]

#Three WHO subtypes of adult-type diffuse gliomas
IDHmut_codel <- intersect(rownames(clinicalData)[which(clinicalData$WHO_2021_subtype == "IDHmut-codel")], overlap.Sample)
IDHmut_noncodel <- intersect(rownames(clinicalData)[which(clinicalData$WHO_2021_subtype == "IDHmut-noncodel")], overlap.Sample)
IDHwt <- intersect(rownames(clinicalData)[which(clinicalData$WHO_2021_subtype == "IDHwt")], overlap.Sample)

#correlation
corFun <- function(mat.1, mat.2, s1, s2, s3){
  
  cor.mat <- matrix(NA, nrow = nrow(mat.1), 3)
  
  i <- 1
  for (i in 1:nrow(mat.1)) {
    cor.mat[i,1] <- cor(as.numeric(mat.1[i, IDHmut_codel]), as.numeric(mat.2[i, IDHmut_codel]), method = "spearman")
    cor.mat[i,2] <- cor(as.numeric(mat.1[i, IDHmut_noncodel]), as.numeric(mat.2[i, IDHmut_noncodel]), method = "spearman")
    cor.mat[i,3] <- cor(as.numeric(mat.1[i, IDHwt]), as.numeric(mat.2[i, IDHwt]), method = "spearman")
  }
  
  return(cor.mat)
  
}

cna_rna_cor <- corFun(cna.mat, RNA.mat.fpkm_UQ, IDHmut_codel, IDHmut_noncodel, IDHwt)
rna_protein_cor <- corFun(RNA.mat.fpkm_UQ, protein.mat, IDHmut_codel, IDHmut_noncodel, IDHwt)
cna_protein_cor <- corFun(cna.mat, protein.mat, IDHmut_codel, IDHmut_noncodel, IDHwt)

#plot
plotData <- data.frame(x = rep(rep(c("IDHmut_codel","IDHmut_noncodel","IDHwt"), each = nrow(cna_rna_cor)), time = 3), 
                       y = c(cna_rna_cor, rna_protein_cor, cna_protein_cor), 
                       group = rep(c("cna_rna","rna_protein","cna_protein"),each = nrow(cna_rna_cor)*3))
plotData$x <- fct_inorder(plotData$x)


ggplot(plotData, aes(x = x, y = y, fill = group)) + 
  geom_boxplot(outlier.colour = NA) + 
  scale_fill_manual(values = rev(brewer.pal(3,"Set1"))) + 
  scale_y_continuous(breaks = seq(-0.8,0.8,0.4),limits = c(-0.8,0.8)) + 
  facet_wrap(~x, scales = "free_x",ncol = 3) + guides(fill = FALSE) + 
  theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("")




