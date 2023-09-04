#----- library -----#

library(maftools)
library(ggplot2)
library(NMF)
library(stringr)
library(RColorBrewer)
library(SomaticSignatures)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)

#----- load data -----#

#MAF file and clinical information from two Chinese cohorts were merged
cohort.mut <- read.maf(maf ="cohort.mutation", clinicalData = "cohort.clinicalInfor.txt",
                    vc_nonSyn = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
                                  "Nonsense_Mutation", "Nonstop_Mutation", "Multi_Hit", "In_Frame_Del", 
                                  "In_Frame_Ins", "Missense_Mutation", "Silent"))

#TMB
TMB <- cohort.mut@variants.per.sample$Variants/38
names(TMB) <- as.character(cohort.mut@variants.per.sample$Tumor_Sample_Barcode)
cohort.mut@clinical.data$TMB <- TMB[as.character(cohort.mut@clinical.data$Tumor_Sample_Barcode)]
hypermut <- names(TMB)[which(TMB>50)]


#----- Figure S2H -----#

#Extract mutational signatures
glioma.tnm <- trinucleotideMatrix(maf = cohort.mut, add = TRUE,
                                  ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

#glioma.sign <- estimateSignatures(mat = glioma.tnm, nTry = 15, nrun = 10)
glioma.sign <- extractSignatures(mat = glioma.tnm, n = 3)
glioma.cos <- compareSignatures(nmfRes = glioma.sign, sig_db = "SBS")

pheatmap::pheatmap(mat = glioma.cos$cosine_similarities, 
                   cluster_rows = FALSE, main = "cosine similarity against validated signatures",
                   cellwidth = 5,cellheight = 12,fontsize_row = 6,
                   fontsize_col = 6,fontsize = 8,treeheight_col = 20)


#----- Figure S2G -----#

#Plot mutation signatures
prob.mat <- glioma.sign$signatures
plotData <- as.data.frame(t(prob.mat))
nsigs <- nrow(plotData)

col.mat <- c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
cols.mat <- rep(col.mat , each=16)

par(mfrow = c(1,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,2.5,0) + 0.1, las = 1, 
    tcl = -.25, font.main = 4, xpd = NA)

ae <- glioma.cos$best_match[[3]][1] #best_match[[1:3]]
ae.bm <- glioma.cos$best_match[[3]][2]
ae <- paste0(ae, " \n Aetiology: ", ae.bm)
prob.data <- as.matrix(plotData[3,]) #plotData[1:3,]

bh <- 0.2
colnames(prob.data) <- paste(paste(str_sub(colnames(prob.data),1,1), 
                                   str_sub(colnames(prob.data),5,5), sep = ""), 
                             str_sub(colnames(prob.data),7,7),sep = "")

barplot(prob.data, yaxt = "n", col = cols.mat, beside = TRUE, ylim = c(-0.1, bh),
        border = NA, las = "2", cex.names = 0.1)

title(main = ae, cex.main = 0.8, line = 0, font.main = 1)
axis(side = 2, at = seq(0, bh, 0.05), pos = -2, las = 2, lwd = 2, hadj = 1.1,
     font = 1, cex.axis = 1)
rect(xleft = seq(0, 192, 32), ybottom = -0.05, xright = 192, ytop = -0.02, 
     col = col.mat, border = 'gray70')
text(labels = c("C>A","C>G","C>T","T>A","T>C","T>G"),
     y = rep(-0.1,6),x = seq(0, 192, 32)[2:7]-16, cex = 0.5,
     font = 0.2, font.lab = 1, pos = 3)


#----- Figure S2J -----#

#Contribution
ContributionRate <- glioma.sign$contributions
plotData <- data.frame(y = ContributionRate[1,], #ContributionRate[1:3,]
                       group = cohort.mut@clinical.data$WHO_2021_subtype[match(colnames(ContributionRate), 
                                                                               cohort.mut@clinical.data$Tumor_Sample_Barcode)])

ggplot(plotData, aes(x = group,y = y,fill = group))+
  geom_boxplot(outlier.colour = NA, width = 0.3)+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  geom_jitter(width = 0.05,shape = 16,size=0.8,fill="black")+
  scale_fill_manual(values = rev(brewer.pal(3,"Set1")))+
  ylab("Contribution")+xlab("")+
  guides(fill=FALSE)+
  stat_compare_means(method = "aov")+
  theme_bw()+theme(panel.grid = element_blank())


#----- Figure S2I -----#

glioma.se <- signatureEnrichment(maf = cohort.mut, sig_res = glioma.sign)
sig.assign <- glioma.se$Signature_Assignment
subCount <- cohort.mut@clinical.data
subCount$signature <- sig.assign$Signature[match(subCount$Tumor_Sample_Barcode, sig.assign$Tumor_Sample_Barcode)]

y0 <- table(subCount$WHO_2021_subtype)/nrow(subCount)
y1 <- table(subCount[,c(4,8)])[,1]/table(subCount$signature)[1]
y2 <- table(subCount[,c(4,8)])[,2]/table(subCount$signature)[2]
y3 <- table(subCount[,c(4,8)])[,3]/table(subCount$signature)[3]

plotData <- data.frame(y = c(y0, y1, y2, y3)*100, fill = rep(c("1","2","3"), time = 4),
                       x = rep(c("Control","SBS1","SBS55","SBS15"), each = 3))
plotData$x <- fct_inorder(plotData$x)

ggplot(plotData, aes(x, y, fill = fill))+
  geom_bar(position = "stack", stat = "identity", width = 0.5)+
  scale_fill_manual(values = rev(brewer.pal(3, 'Set1')))+
  scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20))+
  ylab("Percentage (%)")+xlab("")+guides(fill=FALSE)+
  theme_bw()+theme(panel.grid = element_blank())


#----- Figure S2K -----#

#non-synonymous and synonymous mutations from two Chinese cohorts
TMB.all <- read.csv("./TMB.all.csv")
Xvalue <- sig.assign$Signature
Yvalue <- TMB.all$TMB_nonsilent[match(sig.assign$Tumor_Sample_Barcode, TMB.all$Tumor_Sample_Barcode)] #TMB.all$TMB_silent

plotData <- data.frame(x = Xvalue, y = Yvalue)

ggplot(plotData, aes(x = x, y = y, fill = x))+
  geom_boxplot(outlier.colour = NA, width = 0.3)+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,5,1))+
  geom_jitter(width = 0.05, shape = 16, size = 0.8, fill = "black")+
  scale_fill_manual(values = c("#4DAF4A", "#E41A1C", "#377EB8")) + 
  ylab("Mutations / Mb") + xlab("")+
  guides(fill=FALSE) +
  theme_bw() + theme(panel.grid = element_blank())

#T-test
t.test(plotData$y[which(plotData$x == "Signature_3")], plotData$y[which(plotData$x == "Signature_2")])


