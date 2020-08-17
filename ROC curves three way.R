##############################
#ROC CURVES cfRRBS CASE SUBTYPES
##############################

#PREPROCESSING
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")

cluster.light <- readRDS("cluster.light.Rds") 


#binning samples together --> need smaller bins
#btw binnen van de niet CpG eiland geclusterde files lukt niet --> te weinig RAM
library(prospectr)
meth_bins <- binning(t(cluster.light), 1000)

#annotatie van enkel case
library("readxl")
names_case<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=3)
annot <- c(names_case$type)

#create three way annotation
#SCLC
labelscase_SCLC<- annot
#transforming labels to binary division
for (i in 1:length(labelscase_SCLC)){
  if (labelscase_SCLC[i]!="SCLC"){
    labelscase_SCLC[i]="other"
  }
}
#LUAD
labelscase_LUAD<- annot
#transforming labels to binary division
for (i in 1:length(labelscase_LUAD)){
  if (labelscase_LUAD[i]!="LUAD"){
    labelscase_LUAD[i]="other"
  }
}
#LUSC
labelscase_LUSC<- annot
#transforming labels to binary division
for (i in 1:length(labelscase_LUSC)){
  if (labelscase_LUSC[i]!="LUSC"){
    labelscase_LUSC[i]="other"
  }
}

#selecteer enkel de case stalen
meth_bins_case <- meth_bins[match(c(names_case$case),c(rownames(meth_bins))),]
####################################################################################################

###MODELLING
############
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Modelling case")
###model 1:random forest with LOOV
#function for RF with LOOV
exec.RFLOOV <- function(samples, labels, j){
  library(randomForest)
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    rfmodel=randomForest(samples[-i,], factor(labels[-i]),mtry=j, ntree=500) #labels moeten factor zijn voor RF
    rf_prediction=predict(rfmodel, samples[i,], type="prob") 
    predictions <- c(predictions, as.character(rf_prediction)) 
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}
#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_case, labelscase_LUAD, 5)
result=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob[,2]), t1="other")
auc.result_LUAD=simple.auc(result)
print(auc.result_LUAD)
plot.ROC(result)

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_case, labelscase_LUSC, 5)
result=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob[,2]), t1="other")
auc.result_LUSC=simple.auc(result)
print(auc.result_LUSC)
plot.ROC(result)

#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_case, labelscase_SCLC, 5)
result=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob[,2]), t1="SCLC")
auc.result_SCLC=simple.auc(result)
print(auc.result_SCLC)
plot.ROC(result)

#calculate average AUC
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)

#####herhalen met top 20#################
correlations=sapply(1:(dim(meth_bins_case)[2]), function(x)(cor(meth_bins_case[,x],as.numeric(factor(annot)), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
features.top20=correlations.sorted$ix[1:20] 
meth_bins_top20=meth_bins_case[,features.top20]

#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_top20, labelscase_LUAD, 5)
result=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob[,2]), t1="other")
auc.result_LUAD=simple.auc(result)
print(auc.result_LUAD)
plot.ROC(result)

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_top20, labelscase_LUSC, 5)
result=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob[,2]), t1="other")
auc.result_LUSC=simple.auc(result)
print(auc.result_LUSC)
plot.ROC(result)

#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(meth_bins_top20, labelscase_SCLC, 5)
result=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob[,2]), t1="SCLC")
auc.result_SCLC=simple.auc(result)
print(auc.result_SCLC)
plot.ROC(result)

#calculate average AUC
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)


###model 2: BInomial (ridge) regression met LOOV
####################################################
exec.BRLOOV <- function(samples, labels){
  library(glmnet)
  #model
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    cvfit = cv.glmnet(samples[-i,], (factor(labels[-i])), family = "binomial",type = "class", alpha = 0) #labels moeten factor zijn voor RF
    cv_prediction = predict(cvfit, matrix(samples[i,], nrow=1), s = "lambda.min",type = "response")
    predictions<-rbind(predictions, as.character(cv_prediction))
  }
  return(predictions)
}

par(mfrow=c(3,1))
#execute BR LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(meth_bins_top20, labelscase_LUAD)
result=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob), t1="other")
auc.result_LUAD=simple.auc(result)
print(auc.result_LUAD)
plot.ROC(result)

#execute BR LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(meth_bins_top20, labelscase_LUSC)
result=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob), t1="other")
auc.result_LUSC=simple.auc(result)
print(auc.result_LUSC)
plot.ROC(result)

#execute BR LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(meth_bins_top20, labelscase_SCLC)
result=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob), t1="SCLC")
auc.result_SCLC=simple.auc(result)
print(auc.result_SCLC)
plot.ROC(result)

#calculate average AUC
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)

#model 3: LSVM met LOOV
set.seed(2020)
tuned_parameters <- tune.svm(meth_bins_top20,factor(annot), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear", probability= TRUE)

exec.LSVMLOOV <- function(samples, labels){
  library(e1071)
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    svmmodel=svm(samples[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.1) #labels moeten factor zijn voor RF
    svm_prediction_prob=predict(svmmodel, matrix(samples[i,], nrow=1), probability=TRUE)
    predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}

par(mfrow=c(3,1))
#execute LSVM LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob1=exec.LSVMLOOV(meth_bins_top20, labelscase_LUAD)
result=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob1[,2]), t1="other")
auc.result_LUAD=simple.auc(result)
print(auc.result_LUAD)
plot.ROC(result)

#execute LSVM LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob2=exec.LSVMLOOV(meth_bins_top20, labelscase_LUSC)
result=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob2), t1="LUSC")
auc.result_LUSC=simple.auc(result)
print(auc.result_LUSC)
plot.ROC(result)

#execute LSVM LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob3=exec.LSVMLOOV(meth_bins_top20, labelscase_SCLC)
result=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob3), t1="other")
auc.result_SCLC=simple.auc(result)
print(auc.result_SCLC)
plot.ROC(result)

#calculate average AUC
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)




