##############################
#ROC CURVES cfRRBS CASE SUBTYPES
##############################

#PREPROCESSING
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")

cluster.light <- readRDS("cluster.light.Rds") 


#binning samples together --> need smaller bins
#btw binnen van de niet CpG eiland geclusterde files lukt niet --> te weinig RAM
library(prospectr)
meth_bins1 <- binning(t(cluster.light), 1000)
meth_bins2 <- binning((cluster.light.subset), 1000)
meth_bins <- binning((cluster.light.subsetcancer), 1000)

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
prob=exec.RFLOOV(meth_bins_case, annot, 5)
result1=exec.ROC(
  truth=annot, vals1=c(prob[,2]), t1="SCLC")
auc.result1=simple.auc(result1)
print(auc.result1)
plot.ROC(result1)



###model 2: BInomial (ridge) regression met LOOV
exec.BRLOOV <- function(samples, labels, a){
  library(glmnet)
  #model
  predictions <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    cvfit = cv.glmnet(samples[-i,], (factor(labels[-i])), family = "binomial",type = "class", alpha = a) #labels moeten factor zijn voor RF
    cv_prediction = predict(cvfit, matrix(samples[i,], nrow=1), s = "lambda.min",type = "response")
    predictions<-rbind(predictions, as.character(cv_prediction))
  }
  return(predictions)
}

prob=exec.BRLOOV(meth_bins_case, annot, 1)
result21=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result21=simple.auc(result21)
print(auc.result21)
plot.ROC(result21)

prob=exec.BRLOOV(meth_bins_case, annot, 0)
result20=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result20=simple.auc(result20)
print(auc.result20)
plot.ROC(result20)

prob=exec.BRLOOV(meth_bins_case, annot, 0.5)
result25=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result25=simple.auc(result25)
print(auc.result25)
plot.ROC(result25)

#model 3: LSVM met LOOV
exec.LSVMLOOV <- function(samples, labels){
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

prob1=exec.LSVMLOOV(meth_bins_case, annot)
result3=exec.ROC(
  truth=annot, vals1=c(prob1[,2]), t1="SCLC")
auc.result3=simple.auc(result3)
print(auc.result3)
plot.ROC(result3)

set.seed(2020)
tuned_parameters <- tune.svm(meth_bins_case,factor(annot), gamma = 10^(-5:-1), cost = 10^(-3:1), kernel="linear", probability= TRUE)
print(tuned_parameters)
####################################

AUC_all_two_way_zscore<- data.frame(auc.result3, auc.result20, auc.result25, auc.result21, auc.result1)
#create overview image
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\cfRRBS\\modelling\\")
png("./ROC 3 ML together two-wayRatios on cancer genes subset .png",    # create PNG for the heat map        
    width = 8*200,        # 5 x 300 pixels
    height = 5*200,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
plot.ROC.combined(result3, result21, result20 ,result25, result1)
legend(0.6, 0.35,      # location of the legend on the heatmap plot
       legend = c("Linear support vector machine", "LASSO regression", "Ridge regression", "Elastic net regression",  "Random forest"), # category labels
       col = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"),  # color key
       lty= 1,             # line style
       lwd = 3,            # line width
       cex = 1
)

dev.off()

################
#feature selection
prob=exec.RFLOOV.top20(meth_bins_case, annot)
result1=exec.ROC(
  truth=annot, vals1=c(prob[,2]), t1="SCLC")
auc.result1=simple.auc(result1)
print(auc.result1)
plot.ROC(result1)

prob=exec.BRLOOV.top20(meth_bins_case, annot, 1)
result21=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result21=simple.auc(result21)
print(auc.result21)
plot.ROC(result21)

prob=exec.BRLOOV.top20(meth_bins_case, annot, 0)
result20=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result20=simple.auc(result20)
print(auc.result20)
plot.ROC(result20)

prob=exec.BRLOOV.top20(meth_bins_case, annot, 0.5)
result25=exec.ROC(
  truth=annot, vals1=c(prob), t1="SCLC")
auc.result25=simple.auc(result25)
print(auc.result25)
plot.ROC(result25)

prob1=exec.LSVMLOOV.top20(meth_bins_case, annot)
result3=exec.ROC(
  truth=annot, vals1=c(prob1[,2]), t1="SCLC")
auc.result3=simple.auc(result3)
print(auc.result3)
plot.ROC(result3)

#create overview image
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\cfRRBS\\modelling\\")
png("./ROC 3 ML together two-wayRatios top 20.png",    # create PNG for the heat map        
    width = 8*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size
par(cex.main=2.5)
plot.ROC.combined(result1, result25, result3)
legend("right",      # location of the legend on the heatmap plot
       legend = c("RF", "BR (a=0.5)", "LSVM"), # category labels
       col = c("dark blue", "dark green", "dark red"),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex = 2
)

dev.off()
