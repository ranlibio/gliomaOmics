
#----- library -----#

library(maftools)
library(TCGAbiolinks)
library(ggplot2)
library(RColorBrewer)
library(forcats)

#----------#

top10 <- c("TP53", "TTN", "CIC", "PIK3R1", "ATRX", "NF1", "EGFR", "FUBP1", "PIK3CA", "PTEN")
mut.class <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
               "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")


#----- TCGA -----#

#load data
load("mutLGG.rda")
load("mutGBM.rda")

gliomaMut <- rbind(LGGmut, GBMmut)
gliomaMut <- gliomaMut[which(gliomaMut$Variant_Classification %in% mut.class),]

tcgaIDHmut <- unique(gliomaMut$Tumor_Sample_Barcode[which(gliomaMut$Hugo_Symbol == "IDH1" | 
                                                            gliomaMut$Hugo_Symbol == "IDH2")])
tcgaIDHwt <- setdiff(unique(gliomaMut$Tumor_Sample_Barcode), tcgaIDHmut)

tcgaIDHmut <- gliomaMut[which(gliomaMut$Tumor_Sample_Barcode %in% tcgaIDHmut),]
tcgaIDHwt <- gliomaMut[which(gliomaMut$Tumor_Sample_Barcode %in% tcgaIDHwt),]

#mutation frequency
tcga.mut.fre <- matrix(0, 10, 2)
rownames(tcga.mut.fre) <- top10
colnames(tcga.mut.fre) <- c("IDHmut", "IDHwt")

i <- 1
for (i in 1:nrow(tcga.mut.fre)) {
  
  n1 <- tcgaIDHmut$Tumor_Sample_Barcode[which(tcgaIDHmut$Hugo_Symbol == rownames(tcga.mut.fre)[i])]
  n2 <- length(unique(tcgaIDHmut$Tumor_Sample_Barcode))
  tcga.mut.fre[i,1] <- length(unique(n1))/n2
  
  n3 <- tcgaIDHwt$Tumor_Sample_Barcode[which(tcgaIDHwt$Hugo_Symbol == rownames(tcga.mut.fre)[i])]
  n4 <- length(unique(tcgaIDHwt$Tumor_Sample_Barcode))
  tcga.mut.fre[i,2] <- length(unique(n3))/n4
}


#----- CGGA -----#

#load data
cggaMut <- read.table("./CGGA.WEseq_286.20200506.txt", row.names = 1, header = TRUE, sep = "\t")
cggaMut <- cggaMut[,-which(colSums(is.na(cggaMut)) == nrow(cggaMut))]

cggaIDHmut <- unique(colnames(cggaMut)[which(cggaMut["IDH1",] != "" | cggaMut["IDH2",] != "")])
cggaIDHwt <- setdiff(unique(colnames(cggaMut)), cggaIDHmut)
cggaIDHmut <- cggaMut[,cggaIDHmut]
cggaIDHwt <- cggaMut[,cggaIDHwt]

#mutation frequency
cgga.mut.fre <- matrix(0, 10, 2)
rownames(cgga.mut.fre) <- top10
colnames(cgga.mut.fre) <- c("IDHmut","IDHwt")

i <- 1
for (i in 1:nrow(cgga.mut.fre)) {
  
  n1 <- which(cggaIDHmut[rownames(cgga.mut.fre)[i],] != "")
  n2 <- ncol(cggaIDHmut)
  cgga.mut.fre[i,1] <- length(unique(n1))/n2
  
  n3 <- which(cggaIDHwt[rownames(cgga.mut.fre)[i],] != "")
  n4 <- ncol(cggaIDHwt)
  cgga.mut.fre[i,2] <- length(unique(n3))/n4
}


#----- Cohort 1 -----#

cohort1.mut <- read.maf(maf = "./cohort1.mutation", clinicalData = "./cohort1.clinicalInfor.txt", 
                        vc_nonSyn = mut.class)

cohort1.TMB <- cohort1.mut@variants.per.sample$Variants/38
names(cohort1.TMB) <- as.character(cohort1.mut@variants.per.sample$Tumor_Sample_Barcode)
cohort1.mut@clinical.data$cohort1.TMB <- cohort1.TMB[as.character(cohort1.mut@clinical.data$Tumor_Sample_Barcode)]
hypermut.cohort1 <- names(cohort1.TMB)[which(cohort1.TMB > 50)]

cohort1IDHmut <- cohort1.mut@clinical.data$Tumor_Sample_Barcode[which(cohort1.mut@clinical.data$WHO_subtype != "IDHwt")]
cohort1IDHmut <- setdiff(cohort1IDHmut,hypermut.cohort1)
cohort1IDHwt <- cohort1.mut@clinical.data$Tumor_Sample_Barcode[which(cohort1.mut@clinical.data$WHO_subtype == "IDHwt")]
cohort1IDHwt <- setdiff(cohort1IDHwt,hypermut.cohort1)

cohort1IDHmut <- subsetMaf(maf = cohort1.mut, tsb = cohort1IDHmut)
cohort1IDHwt <- subsetMaf(maf = cohort1.mut, tsb = cohort1IDHwt)

#mutation frequency
cohort1.mut.fre <- matrix(0, 10, 2)
rownames(cohort1.mut.fre) <- top10
colnames(cohort1.mut.fre) <- c("IDHmut", "IDHwt")

i <- 1
for (i in 1:nrow(cohort1.mut.fre)) {
  
  n1 <- cohort1IDHmut@data$Tumor_Sample_Barcode[which(cohort1IDHmut@data$Hugo_Symbol == 
                                                        rownames(cohort1.mut.fre)[i])]
  n2 <- length(unique(cohort1IDHmut@data$Tumor_Sample_Barcode))
  cohort1.mut.fre[i,1] <- length(unique(n1))/n2
  
  n3 <- cohort1IDHwt@data$Tumor_Sample_Barcode[which(cohort1IDHwt@data$Hugo_Symbol == 
                                                       rownames(cohort1.mut.fre)[i])]
  n4 <- length(unique(cohort1IDHwt@data$Tumor_Sample_Barcode))
  cohort1.mut.fre[i,2] <- length(unique(n3))/n4
}


#----- Cohort 2 -----#

cohort2.mut <- read.maf(maf = "./cohort2.mutation", clinicalData = "./cohort2.clinicalInfor.txt", 
                        vc_nonSyn = mut.class)

cohort2.TMB <- cohort2.mut@variants.per.sample$Variants/38
names(cohort2.TMB) <- as.character(cohort2.mut@variants.per.sample$Tumor_Sample_Barcode)
cohort2.mut@clinical.data$cohort2.TMB <- cohort2.TMB[as.character(cohort2.mut@clinical.data$Tumor_Sample_Barcode)]
hypermut.cohort2 <- names(cohort2.TMB)[which(cohort2.TMB > 50)]

cohort2IDHmut <- cohort2.mut@clinical.data$Tumor_Sample_Barcode[which(cohort2.mut@clinical.data$WHO_subtype != "IDHwt")]
cohort2IDHmut <- setdiff(cohort2IDHmut,hypermut.cohort2)
cohort2IDHwt <- cohort2.mut@clinical.data$Tumor_Sample_Barcode[which(cohort2.mut@clinical.data$WHO_subtype == "IDHwt")]
cohort2IDHwt <- setdiff(cohort2IDHwt,hypermut.cohort2)

cohort2IDHmut <- subsetMaf(maf = cohort2.mut, tsb = cohort2IDHmut)
cohort2IDHwt <- subsetMaf(maf = cohort2.mut, tsb = cohort2IDHwt)

cohort2.mut.fre <- matrix(0, 10, 2)
rownames(cohort2.mut.fre) <- top10
colnames(cohort2.mut.fre) <- c("IDHmut", "IDHwt")

i <- 1
for (i in 1:nrow(cohort2.mut.fre)) {
  
  n1 <- cohort2IDHmut@data$Tumor_Sample_Barcode[which(cohort2IDHmut@data$Hugo_Symbol == 
                                                        rownames(cohort2.mut.fre)[i])]
  n2 <- length(unique(cohort2IDHmut@data$Tumor_Sample_Barcode))
  cohort2.mut.fre[i,1] <- length(unique(n1))/n2
  
  n3 <- cohort2IDHwt@data$Tumor_Sample_Barcode[which(cohort2IDHwt@data$Hugo_Symbol == 
                                                       rownames(cohort2.mut.fre)[i])]
  n4 <- length(unique(cohort2IDHwt@data$Tumor_Sample_Barcode))
  cohort2.mut.fre[i,2] <- length(unique(n3))/n4
}


#----- Figure 1E -----#

#Plot barplot
yValue <- c(cohort1.mut.fre[,1], cohort2.mut.fre[,1], cgga.mut.fre[,1], tcga.mut.fre[,1]) #cohort1.mut.fre[,1:]
xValue <- rep(top10, time = 4)

plotData <- data.frame(x = xValue, y = yValue, group = rep(c("Cohort1", "Cohort2", "CGGA", "TCGA"), each = 10))
plotData$x <- fct_inorder(plotData$x)
plotData$group <- fct_inorder(plotData$group)

ggplot(plotData, aes(x = x, y = y, fill = group)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.5) + 
  facet_wrap(~x, scales = "free_x", ncol = 5) + 
  scale_fill_manual(values = brewer.pal(4,"Set1"))+ 
  labs(x = '', y = 'Percentage (%)') + guides(fill = FALSE) + 
  theme_bw() + theme(panel.grid = element_blank())


#Fisher-test