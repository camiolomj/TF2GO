library("limma")

###Compare across disease severity

###dsev_SARP_BEC is a factor containing disease severity classifications where HC = healthy controls, MMA = mild to moderate asthma and SA = severe asthma
###SARP_sex_BEC is a factor containing biologic sex assigned at birth
###SARP_ICS_BEC is a factor containing ICS use information

#Create a contrast matrix
design <- model.matrix(~0+dsev_SARP_BEC+SARP_sex_BEC+SARP_ICS_BEC)

contr.matrix <- makeContrasts(
  HCvsMMA = dsev_SARP_BECMMA - dsev_SARP_BECHC, 
  MMAvsSA = dsev_SARP_BECSA - dsev_SARP_BECMMA, 
  HCvsSA = dsev_SARP_BECSA - dsev_SARP_BECHC, 
  levels = colnames(design))

#Perform DEG analysis
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast, proportion = 0.01)

summary(decideTests(fit, method="separate", adjust.method = "fdr"))
DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")

#Create a matrix of differentially expressed genes and their expression values across participants
asthma.var.SARP.BEC <- as.matrix(SARP_BEC_data[which(DEG_index[,1]!=0|DEG_index[,2]!=0|DEG_index[,3]!=0),])

#Output the Venn diagram illustrated in Figure 1A
png("DEG venn.png", height = 2000, width = 2000, res = 300)
vennDiagram(vennCounts(DEG_index, include="both"), include="both")
dev.off()
