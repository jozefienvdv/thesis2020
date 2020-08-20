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
exec.RFLOOV<- function(samples, labels){
  library(randomForest)
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    rfmodel=randomForest(samples[-i,], factor(labels[-i]),mtry=5, ntree=500) #labels moeten factor zijn voor RF
    rf_prediction_prob=predict(rfmodel, samples[i,], type="prob") #probabilities van alle drie de typen
    predictions_all <-rbind(predictions_all, rf_prediction_prob)#probabilities bijhouden
  }
  return(predictions_all)
}
#model 2: BInomial (ridge) regression met LOOV
exec.BRLOOV<- function(samples, labels, a){
  library(glmnet)
  #model
  predictions <- c()
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    cvfit = cv.glmnet(samples[-i,], (factor(labels[-i])), family = "binomial", type.measure = "class",  alpha = a) #labels moeten factor zijn voor RF
    cv_prediction_prob = predict(cvfit, matrix(samples[i,], nrow=1), s = "lambda.min", type="response")
    predictions_all <-rbind(predictions_all, cv_prediction_prob)
  }
  return(predictions_all)
}

#model 3: LSVM met LOOV
exec.LSVMLOOV<- function(samples, labels){
  library(e1071)
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    svmmodel=svm(samples[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.001) #labels moeten factor zijn voor RF
    svm_prediction_prob=predict(svmmodel, matrix(samples[i,], nrow=1), probability=TRUE)
    predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}

################################
#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_LUAD)
result11=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob[,2]), t1="other")
auc.result_LUAD=simple.auc(result11)
print(auc.result_LUAD)
plot.ROC(result11)

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_LUSC)
result21=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob[,2]), t1="other")
auc.result_LUSC=simple.auc(result21)
print(auc.result_LUSC)
plot.ROC(result21)

#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_SCLC)
result31=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob[,2]), t1="SCLC")
auc.result_SCLC=simple.auc(result31)
print(auc.result_SCLC)
plot.ROC(result31)


#calculate average AUC
avgres_RF = exec.combineROC(result11, result21, result31)
auc.result.avg.RF=simple.auc(avgres_RF)
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)


######################################
#execute BR LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUAD, 1)
result12=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob), t1="other")
auc.result_LUAD=simple.auc(result12)
print(auc.result_LUAD)
plot.ROC(result12)


#execute BR LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(lcase_bins_full, labelscase_LUSC, 1)
result22=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob), t1="LUSC")
auc.result_LUSC=simple.auc(result22)
print(auc.result_LUSC)
plot.ROC(result22)

#execute BR LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(lcase_bins_full, labelscase_SCLC, 1)
result32=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob), t1="SCLC")
auc.result_SCLC=simple.auc(result32)
print(auc.result_SCLC)
plot.ROC(result32)

#calculate average AUC
avgres_MR_0 = exec.combineROC(result12, result22, result32)
plot.ROC(avgres_MR_0)
auc.result.avg.MR_0=simple.auc(avgres_MR_0)
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)

avgres_MR_1 = exec.combineROC(result12, result22, result32)
plot.ROC(avgres_MR_1)
auc.result.avg.MR_1=simple.auc(avgres_MR_1)
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)

avgres_MR_05 = exec.combineROC(result12, result22, result32)
plot.ROC(avgres_MR_05)
auc.result.avg.MR_05=simple.auc(avgres_MR_05)
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)
##########################################
set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(lcase_bins_full,factor(labelscase), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear", probability= TRUE)
print(tuned_parameters) #3fill in in model!!

#execute LSVM LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
prob1=exec.LSVMLOOV(lcase_bins_full, labelscase_LUAD)
result1=exec.ROC(
  truth=labelscase_LUAD, vals1=c(prob1[,2]), t1="other")
auc.result_LUAD=simple.auc(result1)
print(auc.result_LUAD)
plot.ROC(result1)

#execute LSVM LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
prob2=exec.LSVMLOOV(lcase_bins_full, labelscase_LUSC)
result2=exec.ROC(
  truth=labelscase_LUSC, vals1=c(prob2[,1]), t1="LUSC")
auc.result_LUSC=simple.auc(result2)
print(auc.result_LUSC)
plot.ROC(result2)

#execute LSVM LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)
prob3=exec.LSVMLOOV(lcase_bins_full, labelscase_SCLC)
result3=exec.ROC(
  truth=labelscase_SCLC, vals1=c(prob3[,1]), t1="other")
auc.result_SCLC=simple.auc(result3)
print(auc.result_SCLC)
plot.ROC(result3)

#calculate average AUC
avgres_LSVM = exec.combineROC(result1, result2, result3)
plot.ROC(avgres_LSVM)
auc.result.avg.LSVM=simple.auc(avgres_LSVM)
avgAUC=(auc.result_LUAD + auc.result_LUSC + auc.result_SCLC)/3
print(avgAUC)


############aafbeeldingen

#!!!! NOG FOUTEN IN, mis opgeslaan, juiste versie zie feature sel script!!!!!!!!!!!!!


AUC_all_zscore <- data.frame(auc.result.avg.LSVM, auc.result.avg.MR_0, auc.result.avg.MR_05, auc.result.avg.MR_1, auc.result.avg.RF)

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\modelling\\")
png("./ROC 5 ML together Z-score segments 26 juli.png",    # create PNG for the heat map        
    width = 8*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot.ROC.combined(avgres_LSVM, avgres_MR_1, avgres_MR_0,avgres_MR_05,avgres_RF)
legend("right",      # location of the legend on the heatmap plot
       legend = c("Linear support vector machine", "LASSO regression", "Ridge regression", "Elastic net regression",  "Random forest"), # category labels
       col = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),  # color key
       lty= 1,             # line style
       lwd = 5,            # line width
       cex = 1.2
)
dev.off()

#create overview image for best method --> BR
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\modelling\\")
png("./ROC RF z-scores.png",    # create PNG for the heat map        
    width = 8*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot.ROC.combined(result11, result21, result31, title=expression(bold("ROC curves - Random forest - Z-scores")))
legend("right",      # location of the legend on the heatmap plot
       legend = c("LUAD vs. others", "LUSC vs. others", "SCLC vs. others"), # category labels
       col = c("#66c2a5", "#fc8d62", "#8da0cb"),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex = 2
)

dev.off()

