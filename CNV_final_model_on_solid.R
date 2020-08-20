#################
#finaal model op solid stalen
###################


solidcase_x<-na.omit(solidcase_seg) #!!!!
solidcase_x <- sapply(solidcase_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
solidcase_tot <- data.frame(solidcase_x)
#binning samples together --> need smaller bins
library(prospectr)
scase_bins <- binning(t(solidcase_tot), 1000)
#get annotation of sample subtype
library("readxl")
caseinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=2)
#remove samples without annotation 
scase_bins_full <- scase_bins[match(caseinfo$CNV.ID[45:77],rownames(scase_bins)),]



#get binary annotation
labelscase_celltype <- c(caseinfo$type[45:77])
#transforming labels to binary division
for (i in 1:length(labelscase_celltype)){
  if (labelscase_celltype[i]!="SCLC"){
    labelscase_celltype[i]="NSCLC"
  }
}

#get binary annotation
labelscase_celltype_l <- c(caseinfo$type[1:44])
#transforming labels to binary division
for (i in 1:length(labelscase_celltype_l)){
  if (labelscase_celltype_l[i]!="SCLC"){
    labelscase_celltype_l[i]="NSCLC"
  }
}
#########################################
#create model based on LB, make ROC with solid
samples=lcase_bins_full
labels = labelscase_celltype_l
set.seed(2020)
rfmodel=randomForest(samples, factor(labels),mtry=5, ntree=500) #labels moeten factor zijn voor RF
rf_prediction_prob=predict(rfmodel, scase_bins_full, type="prob") #probabilities van alle drie de typen
result11=exec.ROC(
  truth=labelscase_celltype, vals1=c(rf_prediction_prob[,2]), t1="SCLC")
auc.result_RF=simple.auc(result11)
print(auc.result_RF)
plot.ROC(result11)


 
 
 #create overview image for best method --> BR
 setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\modelling\\")
 png("./ROC solid case RF .png",    # create PNG for the heat map        
     width = 8*200,        # 5 x 300 pixels
     height = 5*200,
     res = 400,            # 300 pixels per inch
     pointsize = 5)        # smaller font size
 plot.ROC.final(result11)
 
 dev.off()
 
 
 
 ################
 ##acc
 #get binary annotation
 labelscase_celltype <- c(caseinfo$type[45:77])
 #transforming labels to binary division
 for (i in 1:length(labelscase_celltype)){
   if (labelscase_celltype[i]!="SCLC"){
     labelscase_celltype[i]="NSCLC"
   }
 }
 samples=lcase_bins_full
 labels = labelscase_celltype_l
set.seed(2020)
rfmodel=randomForest(samples, factor(labels),mtry=5, ntree=500) #labels moeten factor zijn voor RF
rf_prediction=predict(rfmodel, scase_bins_full, type="class") #class van alle drie de typen

RF_accuracy=sum(rf_prediction==labelscase_celltype)/length(rf_prediction)
print(RF_accuracy)
 
 