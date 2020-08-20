##########################
##ACC cfRRBS 12 aug############
#############################

#PREPROCESSING
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")

cluster.light <- readRDS("cluster.light.Rds") 


#binning samples together --> need smaller bins
#btw binnen van de niet CpG eiland geclusterde files lukt niet --> te weinig RAM
library(prospectr)
meth_bins1 <- binning(t(cluster.light), 1000)
meth_bins2 <- binning((cluster.light.subset), 1000)
meth_bins3 <- binning((cluster.light.subsetcancer), 1000)

#annotatie van enkel case
library("readxl")
names_case<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=3)

annot <- c(names_case$type)
#annot van NSCLC vs SCLC
for (i in 1:length(annot)){
  if (annot[i]!="SCLC"){
    annot[i]="NSCLC"
  }
}
meth_bins <-meth_bins2

#selecteer enkel de case stalen
meth_bins_case <- meth_bins[match(c(names_case$case),c(rownames(meth_bins))),]
meth_bins_case <- meth_bins_case[-25,]
annot <- annot[-25]






################
#functions
#model 1:random forest with LOOV  (two subtypes --> ROCR does not support non binary classes)
#function for RF with LOOV
exec.RFLOOV<- function(samples, labels, soort="prob"){
  library(randomForest)
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    rfmodel=randomForest(samples[-i,], factor(labels[-i]),mtry=5, ntree=500) #labels moeten factor zijn voor RF
    if (soort=="prob"){
      rf_prediction_prob=predict(rfmodel, samples[i,], type="prob") #probabilities van alle drie de typen
      predictions_all <-rbind(predictions_all, rf_prediction_prob)#probabilities bijhouden
    }
    else {
      rf_prediction=predict(rfmodel, samples[i,], type="class") #class van alle drie de typen
      predictions_all <-rbind(predictions_all, as.character(rf_prediction))#classes bijhouden
    }
  }
  return(predictions_all)
}
#model 2: BInomial (ridge) regression met LOOV
exec.BRLOOV<- function(samples, labels, a, soort="prob"){
  library(glmnet)
  #model
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    cvfit = cv.glmnet(samples[-i,], (factor(labels[-i])), family = "binomial", type.measure = "class",  alpha = a) #labels moeten factor zijn voor RF
    if (soort=="prob"){
      cv_prediction_prob = predict(cvfit, matrix(samples[i,], nrow=1), s = "lambda.min", type="response")
      predictions_all <-rbind(predictions_all, cv_prediction_prob)
    }
    else {
      cv_prediction= predict(cvfit, matrix(samples[i,], nrow=1), s = "lambda.min", type="class")
      predictions_all <-rbind(predictions_all, cv_prediction)
    }
  }
  return(predictions_all)
}

#model 3: LSVM met LOOV
exec.LSVMLOOV<- function(samples, labels, soort="prob"){
  library(e1071)
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    svmmodel=svm(samples[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.01) 
    if (soort=="prob"){
      svm_prediction_prob=predict(svmmodel, matrix(samples[i,], nrow=1), probability=TRUE)
      predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
    }
    else {
      svm_prediction=predict(svmmodel, matrix(samples[i,], nrow=1))
      predictions <- c(predictions, as.character(svm_prediction)) 
    }
  }
  
  return(predictions)
}

#######################
#TWO WAY
##########################

#execute RF LOOV 
prob=exec.RFLOOV(meth_bins_case, annot, "acc")
RF_accuracy=sum(prob==annot)/length(prob)
print(RF_accuracy)

#alpha =1
#execute BR LOOV 
prob=exec.BRLOOV(meth_bins_case, annot, 1, "acc")
BR_accuracy_1=sum(prob==annot)/length(prob)
print(BR_accuracy_1)

#alpha =0
#execute BR LOOV 
prob=exec.BRLOOV(meth_bins_case, annot, 0, "acc")
BR_accuracy_0=sum(prob==annot)/length(prob)
print(BR_accuracy_0)


#alpha =0.5
#execute BR LOOV 
prob=exec.BRLOOV(meth_bins_case, annot, 0.5, "acc")
BR_accuracy_05=sum(prob==annot)/length(prob)
print(BR_accuracy_05)


set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(meth_bins_case,factor(annot), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear")
print(tuned_parameters) #3fill in in model!!

#execute LSVM LOOV 

prob=exec.LSVMLOOV(meth_bins_case, annot, "acc")
LSVM_accuracy=sum(prob==annot)/length(prob)
print(LSVM_accuracy)

x<- data.frame(annot, prob)


Acc_all_two_way_zscore<- data.frame(LSVM_accuracy, BR_accuracy_0, BR_accuracy_05, BR_accuracy_1, RF_accuracy)
