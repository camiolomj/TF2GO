###Unpack loadings generated from splsrda as detailed in the R script 'Supervised Participant Classification"

loadings<-as.data.frame(abs(splsda.res$loadings$X))

genes_from_pls<-vector()

for(w in 1:ncol(loadings)){
  
  genes_from_pls<-c(genes_from_pls, rownames(loadings)[which(loadings[,w]>0)])
  
}

genes_from_pls<-unique(genes_from_pls)

library(glmnet); require(caret)

#####Elastic net modeling of FEV1 in SARP

FEV<-scale(SARP_fev1_BEC[is.na(SARP_fev1_BEC)==F])
scaled.elastic.data<-as.matrix(scale(t(SARP_BEC_data[,is.na(SARP_fev1_BEC)==F])))
scaled.elastic.data<-scaled.elastic.data[,colnames(scaled.elastic.data)%in%genes_from_pls]

fev.model<-cv.glmnet(x = scaled.elastic.data, y = FEV,
                     penalty.factor = penalty.factor, nfolds = 154, alpha = 1)

#Plot non-zero coefficients from model
coefficients<-as.matrix(coef(fev.model, s = "lambda.min"))
filter.coefficients<-as.matrix(coefficients[which(abs(coefficients)>0),])
coef.order<-order(filter.coefficients)

filter.coefficients<-as.matrix(filter.coefficients[coef.order,])
remove.int<-as.numeric(grep("(Intercept)",rownames(filter.coefficients)))
filter.coefficients<-as.matrix(filter.coefficients[-remove.int,])
filter.coefficients<-cbind(scale(filter.coefficients), 1:length(filter.coefficients))
colnames(filter.coefficients)<-c("Coefficient", "Rank")

png("FEV EN coef SARP.png", height = 5*450, width = 5*800, res = 300)
par(xpd=TRUE)
par(mar= c(5.1, 8.1, 4.1, 5.1))

plot(filter.coefficients[,2]~filter.coefficients[,1] , type = "n", xlab = "", ylab = "Predictor",
     cex.axis = "1", cex.lab = 3, xaxt='n', yaxt='n')

name.vector<-rownames(filter.coefficients)
FEV.pos.egraph<-matrix(ncol=length(name.vector))
FEV.neg.egraph<-matrix(ncol=length(name.vector))
for(i in 1:length(filter.coefficients)){
  if(filter.coefficients[i]>0){FEV.pos.egraph[i]<-name.vector[i]}
  else{FEV.pos.egraph[i]<- " "}}

for(i in 1:length(filter.coefficients)){
  if(filter.coefficients[i]<0){FEV.neg.egraph[i]<-name.vector[i]}
  else{FEV.neg.egraph[i]<- " "}}

text(filter.coefficients, FEV.pos.egraph, col = "red", cex = 1.4)
text(filter.coefficients, FEV.neg.egraph, col = "blue", cex = 1.4)
dev.off()

##SARP model performance evaluation using Leave-one-out 

predicted=vector()
x=scaled.elastic.data
for (i in seq(nrow(x))){
  new.model<-cv.glmnet(x = x[-i,], y = FEV[-i],  
                       family = "gaussian", nfolds = 154)
  predicted[i]=predict(new.model, x, "lambda.min")[i]
}

FEV_data<-as.data.frame(cbind(FEV, predicted))
colnames(FEV_data)<-c("measured", "predicted")

library(ggpubr)

png("EN correlation FEV SARP.png", height = 1000, width = 1000, res = 300)
ggscatter(FEV_data, x = "measured", y = "predicted", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Measured FEV1 %Predicted", ylab = "EN Predicted FEV1 %Predicted")
dev.off()

