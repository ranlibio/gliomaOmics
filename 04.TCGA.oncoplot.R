
#----- library -----#

library(maftools)
library(TCGAbiolinks)

#----- load data -----#

load("mutLGG.rda")
load("mutGBM.rda")

#Merge
gliomaMut <- rbind(LGGmut, GBMmut)
glioma.maf <- read.maf(maf = gliomaMut, isTCGA = TRUE)

tcgaIDHmut <- unique(glioma.maf@data$Tumor_Sample_Barcode[which(glioma.maf@data$Hugo_Symbol == "IDH1" | 
                                                                  glioma.maf@data$Hugo_Symbol == "IDH2")])
tcgaIDHwt <- setdiff(unique(glioma.maf@data$Tumor_Sample_Barcode), tcgaIDHmut)

IDHmut.maf <- subsetMaf(maf = glioma.maf, tsb = tcgaIDHmut)
IDHwt.maf <- subsetMaf(maf = glioma.maf, tsb = tcgaIDHwt)

#----- Figure S2B-----#

top20 <- c("IDH1","TP53","TTN","CIC","MUC16","PIK3R1","OBSCN","NF1","ATRX",
           "USH2A","EGFR","FUBP1","HMCN2","CCDC168","PLEC","APOB","DNAH3","DNAH9",
           "PIK3CA","PTEN")

#Oncoplot
vc_cols <- RColorBrewer::brewer.pal(n = 9, name = 'Set1')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'In_Frame_Del',
  'Splice_Site'
)

oncoplot(maf = IDHmut.maf, fontSize = 0.6, showTumorSampleBarcodes = F, #IDHmut.maf and IDHwt.maf
         genes = top20, sortByAnnotation = TRUE, colors = vc_cols, bgCol = "grey80", 
         annoBorderCol = NA, showTitle=FALSE, keepGeneOrder = TRUE, 
         legendFontSize = 1, draw_titv=TRUE, includeColBarCN = FALSE, 
         drawRowBar = FALSE, sepwd_genes = 0.5, sepwd_samples = 0.5, 
         drawColBar = FALSE, barcode_mar = 5, removeNonMutated = FALSE, 
         annotationFontSize = 0.8, anno_height = 3, groupAnnotationBySize = FALSE)

