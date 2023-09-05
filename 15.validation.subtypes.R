
#----- library -----#

library(CMScaller)
library(CancerSubtypes)
library(survival)
library(survminer)
library(survMisc)

#Validation of proteomic subtypes

#pre-processing
Oh.mat <- read.delim("./protein.mat.Oh.txt", header = TRUE)
Oh.mat$Symbol <- gsub("'", "", Oh.mat$Symbol)
colnames(Oh.mat) <- gsub("X.." ,"", colnames(Oh.mat))
colnames(Oh.mat) <- gsub("\\.", "-", colnames(Oh.mat))

#load signature genes
signature.gene <- read.table("./signature.gene.txt", header = TRUE)

Oh.mat <- Oh.mat[which(Oh.mat$Symbol %in% signature.gene$probe),]
Oh.mat <- Oh.mat[!duplicated(Oh.mat$Symbol),]
rownames(Oh.mat) <- Oh.mat$Symbol
Oh.mat <- Oh.mat[,-c(1)]

#
Oh.mat <- data.imputation(Oh.mat, fun = "median")

#Nearest Template Prediction
ntp.res <- ntp(Oh.mat, signature.gene, distance = "cosine", doPlot = FALSE, 
           nPerm = 1000, seed = 123)

sampleID <- read.table("./sampleID.txt", header = TRUE)
ntp.res <- ntp.res[sampleID$Sample_name,]


#survival
clinicalData <- read.table("./clinicalData.Oh.txt", row.names = 1, header = TRUE)

Sur <- data.frame(clinicalData[rownames(ntp.res), c(2, 1)])
colnames(Sur) <- c("time", "status")
Sur$label <- ntp.res[rownames(Sur), "prediction"]
Sur <- Sur[which(!is.na(Sur$time) & !is.na(Sur$status)),]
Sur$time <- Sur$time/31

Sur.fit <- survfit(Surv(time, status) ~ label, data = Sur)
ggsurv <- ggsurvplot(Sur.fit, pval = TRUE, risk.table = TRUE, pval.method = TRUE)

ggsurv$plot+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
  scale_x_continuous(limits = c(0,80),breaks = seq(0,80,20))



