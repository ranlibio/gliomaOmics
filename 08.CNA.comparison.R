
#----- library -----#

library(ggplot2)

#----- Figure 2E and S3G -----#

#----- CNA matrix -----#

#load data

#mutation matrix from Chinese cohort 1 or cohort 2
cna.mat <- read.table("cna.mat.txt", row.names = 1, header = TRUE, sep = "\t")

clinicalData <- read.table("./cohort.clinicalInfor.txt", header = TRUE, row.names = 1, sep = "\t")

ratio <- matrix(0, 2, ncol(cna.mat))

i <- 1
for (i in 1:ncol(ratio)) {
  ratio[1,i] <- length(which(cna.mat[,i] > 0.2))/nrow(cna.mat)
  ratio[2,i] <- length(which(cna.mat[,i] < (-0.2)))/nrow(cna.mat)
}

plotData <- data.frame(x = clinicalData[colnames(cna.mat),"WHO_2021_subtype"], 
                       y = ratio[1,]*100) #ratio[1:2,]

ggplot(plotData, aes(x = x, y = y, fill = x)) + 
  geom_boxplot(outlier.colour = NA, width=0.5) + 
  geom_jitter(shape = 16, size = 0.8, position = position_jitterdodge(jitter.width = 0.04, dodge.width = 0.5)) + 
  scale_fill_manual(values = rev(brewer.pal(3,"Set1"))) + 
  scale_y_continuous(breaks = seq(0,60,20), limits = c(0,60)) + guides(fill = FALSE) + 
  theme_bw() + theme(panel.grid = element_blank()) + xlab("") + ylab("")



