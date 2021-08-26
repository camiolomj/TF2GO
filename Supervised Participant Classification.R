###top_eigen_genes generates gene lists based on input selection parameters from PCA of an input gene matrix

top_eigen_genes<-function(df, percent_var, num_genes){
  
  if(percent_var>99){
    
    print("Input percent variance lowered to maximum of 100%")
    
    percent_var<-99
    
  }
  
  theoretic_min<-as.numeric(min(FactoMineR::PCA(scale(t(df)), ncp = ncol(df),
                                       graph = F, scale.unit = F)$eig[,3], na.rm = T))
  
  if(percent_var<theoretic_min){
    
    print(paste0("Input percent variance increased to minimum of ", ceiling(theoretic_min), "%"))
    
    percent_var<-ceiling(theoretic_min)
    
  }
  
  
  pick_ncomp<-as.numeric(min(which(FactoMineR::PCA(scale(t(df)), ncp = ncol(df),
                                                   graph = F, scale.unit = F)$eig[,3]>percent_var), na.rm = T))
  
  filtered_eigenvectors<-FactoMineR::PCA(scale(t(df)), ncp = ncol(df),
                                         graph = F, scale.unit = F)$var$contrib[,1:pick_ncomp]
  
  
  top_contrib<-list()
  
  for(c in 1:ncol(filtered_eigenvectors)){
    
    top_contrib[[c]]<-filtered_eigenvectors[order(filtered_eigenvectors[,c], decreasing = T),c][1:num_genes]
    
  }
  
  return(unique(rownames(as.matrix(unlist(top_contrib)))))
}

#Input test parameters for gene list generation
variance_to_test<-c(30,40,50,60,70,80)
gene_number_to_test<-c(3,4,5,6)

#Create combinations of parameters for generating gene of interest sets
input_parameters<-as.list(tidyr::crossing(variance_to_test, gene_number_to_test))

#Filter genes used for initial SARP participant clustering based on pressence in external (IMSA) RNA sequencing data
filtered_SARP_genes<-asthma.var.SARP.BEC[rownames(asthma.var.SARP.BEC)%in%rownames(IMSA_BEC_data),]

#Create a list of input genes for testing based on input parameters
genes_to_test<-list()

for(t in 1:length(input_parameters[[1]])){

  genes_to_test[[t]]<-top_eigen_genes(filtered_SARP_genes, input_parameters[[1]][t],
                                      input_parameters[[2]][t])
}  
  
###test_splsda_BER accepts an input data frame (df), number of components to test during model tuning (max_ncomps),
###gene list for testing (gene) and a class assignment for model training (train_assignment)

test_splsda_BER<-function(df, max_ncomps, genes, train_assignment){
  
  features_for_PLSDA<-as.matrix(scale(t(df[genes,])))
  
  if(ncol(features_for_PLSDA)<=10){
  
    list.keepX <- c(1:ncol(features_for_PLSDA))
  
  }else{
    
    list.keepX <- c(1:10, seq(10, ncol(features_for_PLSDA), 5))
    
  }
  
  if(max_ncomps>ncol(features_for_PLSDA)){
    
    max_ncomps<-as.numeric(ncol(features_for_PLSDA))
    
  }
  
  set.seed(40)
  tune.splsda.BEC <- mixOmics::tune.splsda(features_for_PLSDA, train_assignment, ncomp = max_ncomps, 
                                           validation = 'Mfold', folds = 5, scale = F,
                                           progressBar = TRUE, dist = 'centroids.dist', measure = "BER",
                                           test.keepX = list.keepX, nrepeat = 100)
  
  error <- tune.splsda.BEC$error.rate  # error rate per component for the keepX grid
  choice.ncomp <- tune.splsda.BEC$choice.ncomp$ncomp # optimal number of components based on t-tests
  select.keepX <- tune.splsda.BEC$choice.keepX[1:choice.ncomp]  # optimal number of variables to select
  
  set.seed(40)
  splsda.res <- mixOmics::splsda(features_for_PLSDA, train_assignment, ncomp = choice.ncomp, keepX = select.keepX)
  
  set.seed(40) 
  perf.srbct <- mixOmics::perf(splsda.res, validation = "Mfold", folds = 5,
                     dist = 'centroids.dist', nrepeat = 100,
                     progressBar = T) 
  
  return(min(perf.srbct$error.rate$BER))
  
}

test_num<-length(unique(genes_to_test))

#Calculate model error rates based on the gene lists generated above.  In this case, the train_assignment variable
#is the initial participant cluster identified from the R script 'Participant Clustering and DEG'

error_rates<-mapply(test_splsda_BER, replicate(test_num, list(filtered_SARP_genes)), as.list(rep(5, test_num)),
       unique(genes_to_test), replicate(test_num, list(SARP.clusters)))

###Plot model error versus input parameters from gene list generation

png("BER vs gene number.png", height = 1600, width = 2100, res = 300)
par(mar=c(5.1, 4.1, 4.1, 14.1))
  
  plot(error_rates~input_parameters$gene_number_to_test,
       pch = 19, cex = 2, xaxt="n", 
       ylab = "Model Error Rate",
       xlab = "Genes per Eigenvector",
       col = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(input_parameters$variance_to_test))]
       [match(input_parameters$variance_to_test,
              c(unique(input_parameters$variance_to_test)))],
       xlim = c(min(input_parameters$gene_number_to_test)-0.5,
                max(input_parameters$gene_number_to_test)+0.5)
  )
  axis(1, at = seq(min(input_parameters$gene_number_to_test), 
                   max(input_parameters$gene_number_to_test), by = 1), las=2)
  
  legend("topright", inset=c(-0.65,0.2), legend = paste0("Cumulative Variance = ", unique(input_parameters$variance_to_test)),
         fill = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(input_parameters$variance_to_test))],
         cex = 0.8, xpd = T)
dev.off()



png("BER vs variance.png", height = 1600, width = 2100, res = 300)
par(mar=c(5.1, 4.1, 4.1, 14.1))

plot(error_rates~input_parameters$variance_to_test,
     pch = 19, cex = 2, xaxt="n", 
     ylab = "Model Error Rate",
     xlab = "Genes per Eigenvector",
     col = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(input_parameters$gene_number_to_test))]
     [match(input_parameters$gene_number_to_test,
            c(unique(input_parameters$gene_number_to_test)))],
     xlim = c(min(input_parameters$variance_to_test)-0.5,
              max(input_parameters$variance_to_test)+0.5)
)
axis(1, at = seq(min(input_parameters$variance_to_test), 
                 max(input_parameters$variance_to_test), by = 10), las=2)

legend("topright", inset=c(-0.65,0.2), legend = paste0("Genes per Eigenvector = ", unique(input_parameters$gene_number_to_test)),
       fill = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(input_parameters$gene_number_to_test))],
       cex = 0.8, xpd = T)
dev.off()

###Plot model error versus total number of genes used in class prediction
###Optimal solution should be identified using the "elbow" method based on LOESS line superimposed on plot area

png("total genes vs error rate.png",height = 1600, width = 2100, res = 300)
lw1 <- loess(error_rates~c(unlist(lapply(genes_to_test, length))))
plot(error_rates~c(unlist(lapply(genes_to_test, length))),
     pch = 19, cex = 1.25,
     ylab = "Model Error Rate",
     xlab = "Total Genes Used for Prediction")
j <- order(c(unlist(lapply(genes_to_test, length))))
lines(c(unlist(lapply(genes_to_test, length)))[j],lw1$fitted[j],col="red",lwd=3)
dev.off()

###The selection variable corresponds to the solution number that gives best error minimization at lowest gene number

selection<-11
best_gene_set<-unique(genes_to_test)[selection] 
features_for_PLSDA<-as.matrix(scale(t(SARP_BEC_data[best_gene_set[[1]],])))

###Recapitulate model error results and create a sparse partial least square regression discriminant analysis
###object for class prediction in the validation cohort


#Model tuning and parameter selection
list.keepX <- c(1:ncol(features_for_PLSDA))

set.seed(40)
tune.splsda.BEC <- mixOmics::tune.splsda(features_for_PLSDA, SARP.clusters, ncomp = 5, 
                                         validation = 'Mfold', folds = 5, scale = F,
                                         progressBar = TRUE, dist = 'centroids.dist', measure = "BER",
                                         test.keepX = list.keepX, nrepeat = 100)

error <- tune.splsda.BEC$error.rate  # error rate per component for the keepX grid
choice.ncomp <- tune.splsda.BEC$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX <- tune.splsda.BEC$choice.keepX[1:choice.ncomp]  # optimal number of variables to select

###Create splsda object and measure performance
set.seed(40)
splsda.res <- mixOmics::splsda(features_for_PLSDA, SARP.clusters, ncomp = choice.ncomp, keepX = select.keepX)

set.seed(40) 
perf.splsda <- mixOmics::perf(splsda.res, validation = "Mfold", folds = 5,
                             dist = 'centroids.dist', nrepeat = 100,
                             progressBar = T) 

###Display balanced error rate (BER)
min(perf.splsda$error.rate$BER)

###Generate and print ROC curve
auc.BEC<- mixOmics::auroc(splsda.res, outcome.test = SARP.clusters,
                          multilevel = NULL, plot = TRUE, roc.comp = as.numeric(choice.ncomp))

png("auroc training SARP BEC clusters.png", height = 5*300, width = 5*300, res = 300, pointsize = 16)
auc.BEC$graph.Comp3+ ggplot2::scale_color_manual(values=c(topo.colors(4)))
dev.off()

###Test class prediction on training set
SARP_prediction<-predict(splsda.res, features_for_PLSDA, dist = "all")
SARP.predicted.clusters <- SARP_prediction$MajorityVote$centroids.dist[,length(choice.keepX)] 

table(SARP.predicted.clusters, SARP.clusters)


###Filter the IMSA cohort according to feature selection from training set
IMSA_BEC_data_for_pls<-scale(t(IMSA_BEC_data[rownames(IMSA_BEC_data)%in%colnames(features_for_PLSDA),]))

IMSA_BEC_data_for_pls<-IMSA_BEC_data_for_pls[,match(colnames(features_for_PLSDA),
                                              colnames(IMSA_BEC_data_for_pls))]

###Test class prediction on validation set
IMSA_prediction<-predict(splsda.res, IMSA_BEC_data_for_pls, dist = "all")
IMSA.clusters <- IMSA_prediction$MajorityVote$max.dist[,length(choice.keepX)] 