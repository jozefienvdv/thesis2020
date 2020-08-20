#optimize RF

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

liquidcase_tot <- t(liquidcase_tot)
lcase_tot_full <- liquidcase_tot[match(caseinfo$CNV.ID[1:44],rownames(liquidcase_tot)),]
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
exec.RFLOOV_featurs<- function(samples, labels){
  library(randomForest)
  top20featurs <- data.frame(c(rep(0,20)))
  predictions_all <-c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    rf=randomForest(samples[-i,], factor(labels[-i]),mtry=5, ntree=500, importance=T) #labels moeten factor zijn voor RF
    rf_weights=importance(rf)[,4] #geeft importance waarden terug voor alle featurs
    rf_weights.sortedweights=sort(rf_weights,decreasing=TRUE,index.return=TRUE)
    rf_top20=rf_weights.sortedweights$ix[1:20]
    top20featurs <- cbind(top20featurs, rf_top20)
  }
  return(top20featurs)
}


#execute RF LOOV for LUAD + execute ROC curve and AUC (functions Lennart)
top_LUAD=exec.RFLOOV_featurs(lcase_bins_full, labelscase_LUAD)
y_LUAD<-unlist(top_LUAD[2:45])
z_LUAD <- sort(table(y_LUAD))
df_LUAD <-data.frame(z_LUAD)
df1 <-data.frame(LUAD = as.numeric(levels(df_LUAD$y_LUAD))[df_LUAD$y_LUAD], freq = df_LUAD$Freq)

plot(df, type="h", xlim=c(0,1000))



png(paste0("LUAD top regions.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 3*400,
    res = 400)        # smaller font size
# Basic boxplot
par(new=TRUE)
qplot(as.numeric(levels(df_LUAD$y_LUAD))[df_LUAD$y_LUAD], df_LUAD$Freq, geom=c("line"),xlab="", ylab="frequency")
kp <- plotKaryotype(genome="hg38", plot.type = 3, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp, cex =3)
kpAddChromosomeNames(kp, srt=45, cex=2)
kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
dev.off()


top20_LUAD <- df_LUAD[c((dim(df_LUAD)[1]-20):dim(df_LUAD)[1]),]

#execute RF LOOV for LUSC + execute ROC curve and AUC (functions Lennart)
top_LUSC=exec.RFLOOV_featurs(lcase_bins_full, labelscase_LUSC)
y_LUSC<-unlist(top_LUSC[2:45])
z_LUSC<- sort(table(y_LUSC))
df_LUSC <-data.frame(z_LUSC)
df2 <-data.frame(LUSC= as.numeric(levels(df_LUSC$y_LUSC))[df_LUSC$y_LUSC], freq = df_LUSC$Freq)

plot(df, type="h", xlim=c(0,1000))
top20_LUSC <- df_LUSC[c((dim(df_LUSC)[1]-100):dim(df_LUSC)[1]),]

#execute RF LOOV for SCLC+ execute ROC curve and AUC (functions Lennart)top_LUSC=exec.RFLOOV_featurs(lcase_bins_full, labelscase_LUSC)
top_SCLC=exec.RFLOOV_featurs(lcase_bins_full, labelscase_SCLC)
y_SCLC<-unlist(top_SCLC[2:45])
z_SCLC<- sort(table(y_SCLC))
df_SCLC<-data.frame(z_SCLC)

df3 <-data.frame(SCLC = as.numeric(levels(df_SCLC$y_SCLC))[df_SCLC$y_SCLC], freq = df_SCLC$Freq)

plot(df, type="h", xlim=c(0,1000))

top20_SCLC <- df_SCLC[c((dim(df_SCLC)[1]-100):dim(df_SCLC)[1]),]



overview <- data.frame(top20_LUAD,top20_LUSC,top20_SCLC)
newdf <-data.frame(overview$y_LUAD,overview$y_LUSC, overview$y_SCLC)
sort(table(unlist(newdf)))


#########################################################
###compare RF feature to Pearson cor features

exec.RFLOOV.top20 <- function(samples, labels){
  library(randomForest)
  predictions <- c()
  predictions_all <-c()
  top20featurs <- data.frame(c(rep(0,20)))
  for (i in 1:dim(samples)[1]){
    correlations=sapply(1:(dim(samples[-i,])[2]), function(x)(cor(samples[-i,x],as.numeric(factor(labels[-i])), use = "pairwise.complete.obs", method = "pearson")))#as.numeric: klasselabels omzettten naar numerieke waarden (1 en 2)
    correlations.sorted=sort(abs(correlations),decreasing=T,index.return=T)#abs: gebruik absolute waarden, omdat genen ook antigecorreleerd kunnen zijn
    features.top20=correlations.sorted$ix[1:20] #geeft index van de top 20 meest gecor features
    
    #Now make new version of the data keeping only the top20 features
    samples_top20=samples[,features.top20]
    top20featurs <- cbind(top20featurs, features.top20)
  } 
  return(top20featurs)
}
top_LUAD=exec.RFLOOV.top20(lcase_bins_full, labelscase_LUAD)
y_LUAD<-unlist(top_LUAD[2:45])
z_LUAD <- sort(table(y_LUAD))
df_LUAD <-data.frame(z_LUAD)
df4 <-data.frame(LUAD = as.numeric(levels(df_LUAD$y_LUAD))[df_LUAD$y_LUAD], freq = df_LUAD$Freq)

top_LUSC=exec.RFLOOV.top20(lcase_bins_full, labelscase_LUSC)
y_LUSC<-unlist(top_LUSC[2:45])
z_LUSC<- sort(table(y_LUSC))
df_LUSC <-data.frame(z_LUSC)
df5 <-data.frame(LUSC= as.numeric(levels(df_LUSC$y_LUSC))[df_LUSC$y_LUSC], freq = df_LUSC$Freq)

top_SCLC=exec.RFLOOV.top20(lcase_bins_full, labelscase_SCLC)
y_SCLC<-unlist(top_SCLC[2:45])
z_SCLC<- sort(table(y_SCLC))
df_SCLC<-data.frame(z_SCLC)
df6 <-data.frame(SCLC = as.numeric(levels(df_SCLC$y_SCLC))[df_SCLC$y_SCLC], freq = df_SCLC$Freq)




par(mfrow=c(2,3))
plot(df1, type="h", xlim=c(0,1000), xlab="")
plot(df2, type="h", xlim=c(0,1000), xlab="")
plot(df3, type="h", xlim=c(0,1000), xlab="")
plot(df4, type="h", xlim=c(0,1000), xlab="")
plot(df5, type="h", xlim=c(0,1000), xlab="")
plot(df6, type="h", xlim=c(0,1000), xlab="")


length(intersect(df$LUAD, df2$LUAD))
length(df$LUAD)
length(df2$LUAD)

sect <- intersect(df$LUAD, df2$LUAD)
samples <- lcase_bins_full
samples_int=samples[,sect]


