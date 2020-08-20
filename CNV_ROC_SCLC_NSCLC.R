############################################################
#modelling case 2 subtypes by ratio/z-score + TWO way ROC curves 03/07
###########################################################

###PREPROCESSING###
##################
#choose RATIO or Z-score
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
###############################################################

###ROC CURVES###
################

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
    svmmodel=svm(samples[-i,], factor(labels[-i]), kernel="linear", gamma=1e-05, probability = TRUE, cost=0.01) #labels moeten factor zijn voor RF
    svm_prediction_prob=predict(svmmodel, matrix(samples[i,], nrow=1), probability=TRUE)
    predictions <- c(predictions, attr(svm_prediction_prob, "probabilities")) 
  }
  predictions_all = t(matrix(predictions, nrow=2))
  return(predictions_all)
}
######################################
#execute RF LOOV + execute ROC curve and AUC (functions Lennart)
prob=exec.RFLOOV(lcase_bins_full, labelscase_celltype)
result11=exec.ROC(
  truth=labelscase_celltype, vals1=c(prob[,2]), t1="SCLC")
auc.result_RF=simple.auc(result11)
print(auc.result_RF)
plot.ROC(result11)

#execute BR LOOV + execute ROC curve and AUC (functions Lennart)
prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 0)
result121=exec.ROC(
  truth=labelscase_celltype, vals1=c(prob), t1="SCLC")
auc.result_BR_0=simple.auc(result121)
print(auc.result_BR_0)
plot.ROC(result121)

prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 1)
result122=exec.ROC(
  truth=labelscase_celltype, vals1=c(prob), t1="SCLC")
auc.result_BR_1=simple.auc(result122)
print(auc.result_BR_1)
plot.ROC(result122)

prob=exec.BRLOOV(lcase_bins_full, labelscase_celltype, 0.5)
result123=exec.ROC(
  truth=labelscase_celltype, vals1=c(prob), t1="SCLC")
auc.result_BR_5=simple.auc(result123)
print(auc.result_BR_5)
plot.ROC(result123)

set.seed(2020)
library(e1071)
tuned_parameters <- tune.svm(lcase_bins_full,factor(labelscase_celltype), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear", probability= TRUE)
print(tuned_parameters) #3fill in in model!!

prob1=exec.LSVMLOOV(lcase_bins_full, labelscase_celltype)
result1=exec.ROC(
  truth=labelscase_celltype, vals1=c(prob1[,2]), t1="SCLC")
auc.result_LSVM=simple.auc(result1)
print(auc.result_LSVM)
plot.ROC(result1)


#################afbeelidingen
AUC_all_two_way_zscore<- data.frame(auc.result_LSVM, auc.result_BR_0, auc.result_BR_5, auc.result_BR_1, auc.result_RF)

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\modelling\\")
png("./ROC 5 ML together z-score two way 26 juli.png",    # create PNG for the heat map        
    width = 8*200,        # 5 x 300 pixels
    height = 5*200,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot.ROC.combined(result1, result122, result121, result123, result11)
legend(0.6, 0.35,      # location of the legend on the heatmap plot
       legend = c("Linear support vector machine", "LASSO regression", "Ridge regression", "Elastic net regression",  "Random forest"), # category labels
       col = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),  # color key
       lty= 1,             # line style
       lwd = 3,            # line width
       cex = 1
)
dev.off()




rf=randomForest(lcase_bins_full, factor(labelscase_celltype),mtry=5, ntree=500, importance=TRUE)
rf_weights=importance(rf)[,4] #geeft importance waarden terug voor alle featurs
rf_weights.sortedweights=sort(rf_weights,decreasing=TRUE,index.return=TRUE)
rf_top20=rf_weights.sortedweights$ix[1:20]


correlations=sapply(1:(dim(lcase_bins_full)[2]), function(x)(cor(lcase_bins_full[,x],as.numeric(factor(labelscase_celltype)), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
features.top20=correlations.sorted$ix[1:20] #geeft index van de top 20 meest gecor features

intersect(rf_top20,features.top20)
