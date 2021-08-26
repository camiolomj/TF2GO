##############################SARP Participant Clustering#######################################

library(multiClust)
library(NbClust)

#Create a scaled gene expression data set from the analysis detailed in 'Disease Severity DEG'
sig.BEC.scaled<-scale(t(asthma.var.SARP.BEC))

#Estimate number of clusters to assume for k-means
num_clust<-as.numeric(NbClust(sig.BEC.scaled, diss = NULL, min.nc = 2, max.nc = 6, distance = "euclidean",
        method = "kmeans", index = "dunn")$Best.nc[1])

#Call the multiclust version of Kmeans that allows for alternate distance inputs (e.g. pearson)
BEC.kmeans<-Kmeans(sig.BEC.scaled, num_clust, method = "pearson", iter.max = 500)
SARP.clusters<-factor(BEC.kmeans$cluster)

#Create a contrast matrix
design <- model.matrix(~0+SARP.clusters+SARP_gender_BEC+SARP_ICS_BEC)

contr.matrix <- makeContrasts(
  PC1_vs_PC2 = SARP.clusters1 - SARP.clusters2,
  PC1_vs_PC3 = SARP.clusters1 - SARP.clusters3,
  PC1_vs_PC4 = SARP.clusters1 - SARP.clusters4,
  PC2_vs_PC3 = SARP.clusters2 - SARP.clusters3,
  PC2_vs_PC4 = SARP.clusters2 - SARP.clusters4,
  PC3_vs_PC4 = SARP.clusters3 - SARP.clusters4,
  levels = colnames(design))

#Perform DEG analysis across patient clusters
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast, proportion = 0.01)

summary(decideTests(fit, method="separate", adjust.method = "fdr"))
DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")

###Create lists of genes that are up or down in one specific group by at least 
###two comparisons to other participant clusters.  

pc1.down.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,1]==-1&DEG_index[,2]==-1)|
                                                (DEG_index[,2]==-1&DEG_index[,3]==-1)|
                                                (DEG_index[,1]==-1&DEG_index[,3]==-1)),])

pc1.up.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,1]==1&DEG_index[,2]==1)|
                                              (DEG_index[,2]==1&DEG_index[,3]==1)|
                                              (DEG_index[,1]==1&DEG_index[,3]==1)),])

pc2.down.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,1]==1&DEG_index[,4]==-1)|
                                                (DEG_index[,1]==1&DEG_index[,5]==-1)|
                                                (DEG_index[,4]==-1&DEG_index[,5]==-1)),])

pc2.up.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,1]==-1&DEG_index[,4]==1)|
                                              (DEG_index[,1]==-1&DEG_index[,5]==1)|
                                              (DEG_index[,4]==1&DEG_index[,5]==1)),])


pc3.down.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,2]==1&DEG_index[,4]==1)|
                                                (DEG_index[,2]==1&DEG_index[,6]==-1)|
                                                (DEG_index[,4]==1&DEG_index[,6]==-1)),])

pc3.up.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,2]==-1&DEG_index[,4]==-1)|
                                              (DEG_index[,2]==-1&DEG_index[,6]==1)|
                                              (DEG_index[,4]==-1&DEG_index[,6]==1)),])

pc4.down.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,3]==1&DEG_index[,5]==1)|
                                                (DEG_index[,3]==1&DEG_index[,6]==1)|
                                                (DEG_index[,5]==1&DEG_index[,6]==1)),])

pc4.up.genes<-as.matrix(SARP_BEC_data[which((DEG_index[,3]==-1&DEG_index[,5]==-1)|
                                              (DEG_index[,3]==-1&DEG_index[,6]==-1)|
                                              (DEG_index[,5]==-1&DEG_index[,6]==-1)),])

###Create lists of genes that are up or down in one specific group compared to the remaining 
###pool of participants.  Includes fold change measurements for GSEA

SARP_PC1_v_rest<-factor(ifelse(as.character(SARP.clusters)==1, "1", "rest"))
SARP_PC2_v_rest<-factor(ifelse(as.character(SARP.clusters)==2, "2", "rest"))
SARP_PC3_v_rest<-factor(ifelse(as.character(SARP.clusters)==3, "3", "rest"))
SARP_PC4_v_rest<-factor(ifelse(as.character(SARP.clusters)==4, "4", "rest"))

design <- model.matrix(~0+SARP_PC1_v_rest+SARP_gender_BEC+SARP_ICS_BEC)
contr.matrix <- makeContrasts(
  PC1_vs_rest = SARP_PC1_v_rest1 - SARP_PC1_v_restrest,
  levels = colnames(design))
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast)
GSEA_start<-topTable(fit, p.value = 0.05, number = 5000)
GSEA_start<-GSEA_start[order(GSEA_start$logFC, decreasing = T),]
GOI_fc_list_pc1<-GSEA_start$logFC
names(GOI_fc_list_pc1)<-rownames(GSEA_start)

DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")
pc1.v.rest.genes<-rownames(as.matrix(SARP_BEC_data[which(DEG_index[,1]==1),]))

design <- model.matrix(~0+SARP_PC2_v_rest+SARP_gender_BEC+SARP_ICS_BEC)
contr.matrix <- makeContrasts(
  PC2_vs_rest = SARP_PC2_v_rest2 - SARP_PC2_v_restrest,
  levels = colnames(design))
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast)
GSEA_start<-topTable(fit, p.value = 0.05, number = 5000)
GSEA_start<-GSEA_start[order(GSEA_start$logFC, decreasing = T),]
GOI_fc_list_pc2<-GSEA_start$logFC
names(GOI_fc_list_pc2)<-rownames(GSEA_start)
DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")
pc2.v.rest.genes<-rownames(as.matrix(SARP_BEC_data[which(DEG_index[,1]==1),]))

design <- model.matrix(~0+SARP_PC3_v_rest+SARP_gender_BEC+SARP_ICS_BEC)
contr.matrix <- makeContrasts(
  PC3_vs_rest = SARP_PC3_v_rest3 - SARP_PC3_v_restrest,
  levels = colnames(design))
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast)
GSEA_start<-topTable(fit, p.value = 0.05, number = 5000)
GSEA_start<-GSEA_start[order(GSEA_start$logFC, decreasing = T),]
GOI_fc_list_pc3<-GSEA_start$logFC
names(GOI_fc_list_pc3)<-rownames(GSEA_start)
DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")
pc3.v.rest.genes<-rownames(as.matrix(SARP_BEC_data[which(DEG_index[,1]==1),]))

design <- model.matrix(~0+SARP_PC4_v_rest+SARP_gender_BEC+SARP_ICS_BEC)
contr.matrix <- makeContrasts(
  PC4_vs_rest = SARP_PC4_v_rest4 - SARP_PC4_v_restrest,
  levels = colnames(design))
lm<-lmFit(SARP_BEC_data, design = design, robust = F)
lm.contrast<-contrasts.fit(lm, contrasts=contr.matrix)
fit <- eBayes(lm.contrast)
GSEA_start<-topTable(fit, p.value = 0.05, number = 5000)
GSEA_start<-GSEA_start[order(GSEA_start$logFC, decreasing = T),]
GOI_fc_list_pc4<-GSEA_start$logFC
names(GOI_fc_list_pc4)<-rownames(GSEA_start)
DEG_index <- decideTests(fit, method="separate",adjust.method = "fdr")
pc4.v.rest.genes<-rownames(as.matrix(SARP_BEC_data[which(DEG_index[,1]==1),]))

