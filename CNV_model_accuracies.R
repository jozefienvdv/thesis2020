############################################################
#modelling case 3 subtypes by ratio/z-score + three way ROC curves 03/07
###########################################################


###PREPROCESSING###
##################
#choose RATIO or Z-score or Z-score segments
#removing bins that don't have RATIO for every entry + convert factors to numeric
liquidcase_x<-na.omit(liquidcase) #!!!!
liquidcase_x <- sapply(liquidcase_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
liquidcase_tot <- data.frame(liquidcase_x)

#removing bins that don't have Z-score for every entry + convert factors to numeric
liquidcase_x<-na.omit(liquidcase_z) #!!!!
liquidcase_x <- sapply(liquidcase_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
liquidcase_tot <- data.frame(liquidcase_x)

#removing bins that don't have Z-score  SEGMENT for every entry + convert factors to numeric
liquidcase_x<-na.omit(liquidcase_seg) #!!!!
liquidcase_x <- sapply(liquidcase_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
liquidcase_tot <- data.frame(liquidcase_x)


#binning samples together --> need smaller bins
library(prospectr)
lcase_bins <- binning(t(liquidcase_tot), 1000)

#get annotation of sample subtype
library("readxl")
caseinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=2)
#remove samples without annotation (1 tot 44, want enkel liquid stalen)
lcase_bins_full <- lcase_bins[match(caseinfo$CNV.ID[1:44],rownames(lcase_bins)),]

#annotation
labelscase <- c(caseinfo$type[1:44])
#create three way annotation
#SCLC
labelscase_SCLC<- c(caseinfo$type[1:44])
#transforming labels to binary division
for (i in 1:length(labelscase_SCLC)){
  if (labelscase_SCLC[i]!="SCLC"){
    labelscase_SCLC[i]="other"
  }
}
#LUAD
labelscase_LUAD<- c(caseinfo$type[1:44])
#transforming labels to binary division
for (i in 1:length(labelscase_LUAD)){
  if (labelscase_LUAD[i]!="LUAD"){
    labelscase_LUAD[i]="other"
  }
}
#LUSC
labelscase_LUSC<- c(caseinfo$type[1:44])
#transforming labels to binary division
for (i in 1:length(labelscase_LUSC)){
  if (labelscase_LUSC[i]!="LUSC"){
    labelscase_LUSC[i]="other"
  }
}
###############################################################

###ROC CURVES###
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
    svmmodel=svm(samples[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.1) 
    if (soort=="prob"){
      svm_prediction_prob=predict(svmmodel, matrix(samples[i,], nrow=1), probability=TRUE)
      predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
    }
    else {
      svm_prediction=predict(svmmodel, matrix(samples[i,], nrow=1))
      predictions <- c(predictions, as.character(svm_prediction)) 
    }
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}

################################
#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_LUAD, "acc")
RF_accuracy_LUAD=sum(prob==labelscase_LUAD)/length(prob)
print(RF_accuracy_LUAD)

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_LUSC, "acc")
RF_accuracy_LUSC=sum(prob==labelscase_LUSC)/length(prob)
print(RF_accuracy_LUSC)
#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_SCLC, "acc")
RF_accuracy_SCLC=sum(prob==labelscase_SCLC)/length(prob)
print(RF_accuracy_SCLC)

#calculate average acc
avgacc_RF =(RF_accuracy_LUAD+RF_accuracy_LUSC+RF_accuracy_SCLC)/3
print(avgacc_RF)

######################################
#alpha =1
#execute BR LOOV for LUAD
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUAD, 1, "acc")
BR_accuracy_LUAD_1=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_1)
x<-matrix(c(labelscase_LUAD, prob), nrow = 44, ncol = 2)
#execute BR LOOV for LUSC
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUSC, 1, "acc")
BR_accuracy_LUSC_1=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_1)


#execute BR LOOV for SCLC
prob=exec.BRLOOV(lcase_bins_full, labelscase_SCLC, 1, "acc")
BR_accuracy_SCLC_1=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_1)

#calculate average acc
avgacc_BR_1 =(BR_accuracy_LUAD_1+BR_accuracy_LUSC_1+BR_accuracy_SCLC_1)/3
print(avgacc_BR_1)

#alpha =0
#execute BR LOOV for LUAD
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUAD, 0, "acc")
BR_accuracy_LUAD_0=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_0)

#execute BR LOOV for LUSC
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUSC, 0, "acc")
BR_accuracy_LUSC_0=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_0)


#execute BR LOOV for SCLC
prob=exec.BRLOOV(lcase_bins_full, labelscase_SCLC, 0, "acc")
BR_accuracy_SCLC_0=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_0)

#calculate average acc
avgacc_BR_0 =(BR_accuracy_LUAD_0+BR_accuracy_LUSC_0+BR_accuracy_SCLC_0)/3
print(avgacc_BR_0)

#alpha =0.5
#execute BR LOOV for LUAD
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUAD, 0.5, "acc")
BR_accuracy_LUAD_05=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_05)

#execute BR LOOV for LUSC
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUSC, 0.5, "acc")
BR_accuracy_LUSC_05=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_05)


#execute BR LOOV for SCLC
prob=exec.BRLOOV(lcase_bins_full, labelscase_SCLC, 0.5, "acc")
BR_accuracy_SCLC_05=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_05)

#calculate average acc
avgacc_BR_05=(BR_accuracy_LUAD_05+BR_accuracy_LUSC_05+BR_accuracy_SCLC_05)/3
print(avgacc_BR_05)
##########################################
set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(lcase_bins_full,factor(labelscase), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear")
print(tuned_parameters) #3fill in in model!!

#execute LSVM LOOV for LUAD + execute ROC curve and AUC (functions Lennart)

prob=exec.LSVMLOOV(lcase_bins_full, labelscase_LUAD, "acc")
LSVM_accuracy_LUAD=sum(prob==labelscase_LUAD)/length(prob)
print(LSVM_accuracy_LUAD)

#execute LSVM LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.LSVMLOOV(lcase_bins_full, labelscase_LUSC, "acc")
LSVM_accuracy_LUSC=sum(prob==labelscase_LUSC)/length(prob)
print(LSVM_accuracy_LUSC)

#execute LSVM LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.LSVMLOOV(lcase_bins_full, labelscase_SCLC, "acc")
LSVM_accuracy_SCLC=sum(prob==labelscase_SCLC)/length(prob)
print(LSVM_accuracy_SCLC)

#calculate average acc
avgacc_LSVM =(LSVM_accuracy_LUAD+LSVM_accuracy_LUSC+LSVM_accuracy_SCLC)/3
print(avgacc_LSVM)



############aafbeeldingen
acc_all_zscore<- data.frame(avgacc_LSVM, avgacc_BR_0, avgacc_BR_05, avgacc_BR_1, avgacc_RF)

#####################################################################
###ACC met feature selection ###############

#functions
#model 1:random forest with LOOV  (two subtypes --> ROCR does not support non binary classes)
#function for RF with LOOV
exec.RFLOOV.top20<- function(samples, labels, soort="prob"){
  library(randomForest)
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    correlations=sapply(1:(dim(samples[-i,])[2]), function(x)(cor(samples[-i,x],as.numeric(factor(labels[-i])), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
    correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
    features.top20=correlations.sorted$ix[1:40] #geeft index van de top 20 meest gecor features
    
    #Now make new version of the data keeping only the top20 features
    samples_top20=samples[,features.top20]
    
    
    set.seed(2020)
    rfmodel=randomForest(samples_top20[-i,], factor(labels[-i]),mtry=5, ntree=500) #labels moeten factor zijn voor RF
    if (soort=="prob"){
      rf_prediction_prob=predict(rfmodel, samples_top20[i,], type="prob") #probabilities van alle drie de typen
      predictions_all <-rbind(predictions_all, rf_prediction_prob)#probabilities bijhouden
    }
    else {
      rf_prediction=predict(rfmodel, samples_top20[i,], type="class") #class van alle drie de typen
      predictions_all <-rbind(predictions_all, as.character(rf_prediction))#classes bijhouden
    }
  }
  return(predictions_all)
}
#model 2: BInomial (ridge) regression met LOOV
exec.BRLOOV.top20<- function(samples, labels, a, soort="prob"){
  library(glmnet)
  #model
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    correlations=sapply(1:(dim(samples[-i,])[2]), function(x)(cor(samples[-i,x],as.numeric(factor(labels[-i])), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
    correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
    features.top20=correlations.sorted$ix[1:20] #geeft index van de top 20 meest gecor features
    
    #Now make new version of the data keeping only the top20 features
    samples_top20=samples[,features.top20]
    
    set.seed(2020)
    cvfit = cv.glmnet(samples_top20[-i,], (factor(labels[-i])), family = "binomial", type.measure = "class",  alpha = a) #labels moeten factor zijn voor RF
    if (soort=="prob"){
      cv_prediction_prob = predict(cvfit, matrix(samples_top20[i,], nrow=1), s = "lambda.min", type="response")
      predictions_all <-rbind(predictions_all, cv_prediction_prob)
    }
    else {
      cv_prediction= predict(cvfit, matrix(samples_top20[i,], nrow=1), s = "lambda.min", type="class")
      predictions_all <-rbind(predictions_all, cv_prediction)
    }
  }
  return(predictions_all)
}

#model 3: LSVM met LOOV
exec.LSVMLOOV.top20<- function(samples, labels, soort="prob"){
  library(e1071)
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    correlations=sapply(1:(dim(samples[-i,])[2]), function(x)(cor(samples[-i,x],as.numeric(factor(labels[-i])), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
    correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
    features.top20=correlations.sorted$ix[1:0] #geeft index van de top 20 meest gecor features
    
    #Now make new version of the data keeping only the top20 features
    samples_top20=samples[,features.top20]
    
    
    set.seed(2020)
    svmmodel=svm(samples_top20[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.1) 
    if (soort=="prob"){
      svm_prediction_prob=predict(svmmodel, matrix(samples_top20[i,], nrow=1), probability=TRUE)
      predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
    }
    else {
      svm_prediction=predict(svmmodel, matrix(samples_top20[i,], nrow=1))
      predictions <- c(predictions, as.character(svm_prediction)) 
    }
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}

#############
#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV.top20(lcase_bins_full, labelscase_LUAD, "acc")
RF_accuracy_LUAD=sum(prob==labelscase_LUAD)/length(prob)
print(RF_accuracy_LUAD)

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV.top20(lcase_bins_full, labelscase_LUSC, "acc")
RF_accuracy_LUSC=sum(prob==labelscase_LUSC)/length(prob)
print(RF_accuracy_LUSC)

#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV.top20(lcase_bins_full, labelscase_SCLC, "acc")
RF_accuracy_SCLC=sum(prob==labelscase_SCLC)/length(prob)
print(RF_accuracy_SCLC)

#calculate average acc
avgacc_RF =(RF_accuracy_LUAD+RF_accuracy_LUSC+RF_accuracy_SCLC)/3
print(avgacc_RF)

######################################
#alpha =1
#execute BR LOOV for LUAD
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUAD, 1, "acc")
BR_accuracy_LUAD_1=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_1)
x<-matrix(c(labelscase_LUAD, prob), nrow = 44, ncol = 2)

#execute BR LOOV for LUSC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUSC, 1, "acc")
BR_accuracy_LUSC_1=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_1)


#execute BR LOOV for SCLC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_SCLC, 1, "acc")
BR_accuracy_SCLC_1=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_1)

#calculate average acc
avgacc_BR_1 =(BR_accuracy_LUAD_1+BR_accuracy_LUSC_1+BR_accuracy_SCLC_1)/3
print(avgacc_BR_1)

#alpha =0
#execute BR LOOV for LUAD
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUAD, 0, "acc")
BR_accuracy_LUAD_0=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_0)

#execute BR LOOV for LUSC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUSC, 0, "acc")
BR_accuracy_LUSC_0=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_0)


#execute BR LOOV for SCLC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_SCLC, 0, "acc")
BR_accuracy_SCLC_0=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_0)

#calculate average acc
avgacc_BR_0 =(BR_accuracy_LUAD_0+BR_accuracy_LUSC_0+BR_accuracy_SCLC_0)/3
print(avgacc_BR_0)

#alpha =0.5
#execute BR LOOV for LUAD
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUAD, 0.5, "acc")
BR_accuracy_LUAD_05=sum(prob==labelscase_LUAD)/length(prob)
print(BR_accuracy_LUAD_05)

#execute BR LOOV for LUSC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_LUSC, 0.5, "acc")
BR_accuracy_LUSC_05=sum(prob==labelscase_LUSC)/length(prob)
print(BR_accuracy_LUSC_05)


#execute BR LOOV for SCLC
prob=exec.BRLOOV.top20(lcase_bins_full, labelscase_SCLC, 0.5, "acc")
BR_accuracy_SCLC_05=sum(prob==labelscase_SCLC)/length(prob)
print(BR_accuracy_SCLC_05)

#calculate average acc
avgacc_BR_05=(BR_accuracy_LUAD_05+BR_accuracy_LUSC_05+BR_accuracy_SCLC_05)/3
print(avgacc_BR_05)
##########################################
set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(lcase_bins_full,factor(labelscase), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear")
print(tuned_parameters) #3fill in in model!!

#execute LSVM LOOV for LUAD + execute ROC curve and AUC (functions Lennart)

prob=exec.LSVMLOOV.top20(lcase_bins_full, labelscase_LUAD, "acc")
LSVM_accuracy_LUAD=sum(prob==labelscase_LUAD)/length(prob)
print(LSVM_accuracy_LUAD)

#execute LSVM LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.LSVMLOOV.top20(lcase_bins_full, labelscase_LUSC, "acc")
LSVM_accuracy_LUSC=sum(prob==labelscase_LUSC)/length(prob)
print(LSVM_accuracy_LUSC)

#execute LSVM LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.LSVMLOOV.top20(lcase_bins_full, labelscase_SCLC, "acc")
LSVM_accuracy_SCLC=sum(prob==labelscase_SCLC)/length(prob)
print(LSVM_accuracy_SCLC)

#calculate average acc
avgacc_LSVM =(LSVM_accuracy_LUAD+LSVM_accuracy_LUSC+LSVM_accuracy_SCLC)/3
print(avgacc_LSVM)



############aafbeeldingen
acc_all_zscore<- data.frame(avgacc_LSVM, avgacc_BR_0, avgacc_BR_05, avgacc_BR_1, avgacc_RF)



##########################
#TWO WAY
##########################

#get annotation of sample subtype
library("readxl")
caseinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=2)
#remove samples without annotation (1 tot 44, want enkel liquid stalen)
lcase_bins_full <- lcase_bins[match(caseinfo$CNV.ID[1:44],rownames(lcase_bins)),]

#get annotation
labelscase <- c(caseinfo$type[1:44])

#get binary annotation
labelscase_celltype <- c(caseinfo$type[1:44])
#transforming labels to binary division
for (i in 1:length(labelscase_celltype)){
  if (labelscase_celltype[i]!="SCLC"){
    labelscase_celltype[i]="NSCLC"
  }
}

#execute RF LOOV 
prob=exec.RFLOOV(lcase_bins_full, labelscase_celltype, "acc")
RF_accuracy=sum(prob==labelscase_celltype)/length(prob)
print(RF_accuracy)

#alpha =1
#execute BR LOOV 
prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 1, "acc")
BR_accuracy_1=sum(prob==labelscase_celltype)/length(prob)
print(BR_accuracy_1)

#alpha =0
#execute BR LOOV 
prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 0, "acc")
BR_accuracy_0=sum(prob==labelscase_celltype)/length(prob)
print(BR_accuracy_0)


#alpha =0.5
#execute BR LOOV 
prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 0.5, "acc")
BR_accuracy_05=sum(prob==labelscase_celltype)/length(prob)
print(BR_accuracy_05)


set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(lcase_bins_full,factor(labelscase), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear")
print(tuned_parameters) #3fill in in model!!

#execute LSVM LOOV for LUAD + execute ROC curve and AUC (functions Lennart)

prob=exec.LSVMLOOV(lcase_bins_full, labelscase_celltype, "acc")
LSVM_accuracy=sum(prob==labelscase_celltype)/length(prob)
print(LSVM_accuracy)

Acc_all_two_way_zscore<- data.frame(LSVM_accuracy, BR_accuracy_0, BR_accuracy_05, BR_accuracy_1, RF_accuracy)


################
#acc final model on solid case






