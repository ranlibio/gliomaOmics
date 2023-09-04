
#----- library -----#

library(maftools)
library(ggplot2)

#----- load data -----#

#MAF file and clinical information from Chinese cohort
cohort.mut <- read.maf(maf ="cohort.mutation", clinicalData = "cohort.clinicalInfor.txt",
                       vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", 
                                   "Nonsense_Mutation", "Nonstop_Mutation", "Multi_Hit", "In_Frame_Del", 
                                   "In_Frame_Ins", "Missense_Mutation"))

#TMB
TMB <- cohort.mut@variants.per.sample$Variants/38
names(TMB) <- as.character(cohort.mut@variants.per.sample$Tumor_Sample_Barcode)
cohort.mut@clinical.data$TMB <- TMB[as.character(cohort.mut@clinical.data$Tumor_Sample_Barcode)]
hypermut <- names(TMB)[which(TMB>50)]


#Three WHO subtypes of adult-type diffuse gliomas
IDHmut_codel <- cohort.mut@clinical.data$Tumor_Sample_Barcode[which(cohort.mut@clinical.data$WHO_2021_subtype=="IDHmut-codel")]
IDHmut_noncodel <- cohort.mut@clinical.data$Tumor_Sample_Barcode[which(cohort.mut@clinical.data$WHO_2021_subtype=="IDHmut-noncodel")]
IDHwt <- cohort.mut@clinical.data$Tumor_Sample_Barcode[which(cohort.mut@clinical.data$WHO_2021_subtype=="IDHwt")]

IDHmut_codel.mut <- subsetMaf(maf = cohort.mut, tsb = setdiff(IDHmut_codel,hypermut))
IDHmut_noncodel.mut <- subsetMaf(maf = cohort.mut, tsb = setdiff(IDHmut_noncodel,hypermut))
IDHwt.mut <- subsetMaf(maf = cohort.mut, tsb = setdiff(IDHwt,hypermut))

#Matrix
mut.gene <- c("IDH1", "TP53", "CIC", "ATRX", "FUBP1", "EGFR", "NF1")
oddRatio <- matrix(1, length(mut.gene), 3)
oddRatio.p <- matrix(1, length(mut.gene), 3)

rownames(oddRatio) <- mut.gene
colnames(oddRatio) <- c("IDHmut-codel vs. IDHmut-noncodel","IDHmut-codel vs. IDHwt",
                        "IDHmut-noncodel vs. IDHwt")
rownames(oddRatio.p) <- mut.gene
colnames(oddRatio.p) <- colnames(oddRatio)

#----- Figure 1F and S2E -----#

#Compare mutation frequency
mafCom <- mafCompare(m1 = sub2.mut, m2 = sub3.mut, m1Name = "IDHmut-noncodel", #m1 and m2 were two subtypes of gliomas
                     m2Name = "IDHwt", minMut = 5)

i <- 1
for (i in 1:length(mut.gene)){
  tmp <- which(mafCom$results$Hugo_Symbol == mut.gene[i])
  
  if(length(tmp)>0){
    if(mafCom$results$or[tmp] == "Inf"){
      oddRatio[mut.gene[i],1] <- 500
    }else if(mafCom$results$or[tmp] == 0){
      oddRatio[mut.gene[i],1] <- 0.001
    }else{
      oddRatio[mut.gene[i],1] <- mafCom$results$or[tmp]
    }
    oddRatio.p[mut.gene[i],1] <- mafCom$results$pval[tmp]
  }
}

#
oddRatio <- log2(oddRatio)
oddRatio.p <- -log10(oddRatio.p)

oddRatio <- oddRatio[,c(3,2,1)]
oddRatio.p <- oddRatio.p[,c(3,2,1)]

plot.oddRatio <- c(t(oddRatio))
plot.oddRatio.p <- c(t(oddRatio.p))
plot.x <- rep(mut.gene,each = 3)
plot.y <- rep(colnames(oddRatio),time = 7)

plotData <- data.frame(x = plot.x, y = plot.y, size = plot.oddRatio.p, color = plot.oddRatio,
                       shape = as.character(sign(plot.oddRatio)))
plotData$x <- fct_inorder(plotData$x)
plotData$y <- fct_inorder(plotData$y)

ggplot(plotData, aes(x = x, y = y, size = size, color = color)) + 
  geom_point() + scale_size_continuous(breaks = seq(1,16,3)) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red") + xlab("") + 
  ylab("") + labs(size = "-log10 (P value)", color = "log2 (odds ratio)") + theme_bw()

