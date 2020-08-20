################################
#sensitivity analysis cfRRBS
#################################

###PREPROCESSING###
##################

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")

cluster.light <- readRDS("cluster.light.Rds") 
#verder werken met deze data = light gesmooth (missing data invullen) + geclusterd per CpG eilandje (minder datapunten)
#raw data

raw_na<-na.omit(raw) #!!!!
raw_na<-na.omit(cluster.raw) #!!!!
 #binning samples together --> need smaller bins
#btw binnen van de niet CpG eiland geclusterde files lukt niet --> te weinig RAM
library(prospectr)
meth_bins <- binning(t(cluster.light), 5000)

meth_bins <- binning(cluster.light.subset, 5000)
meth_bins <- binning(t(raw_na), 5000)


#annotation case control
library("readxl")
names_case<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=3)
names_control<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=4)
info_all <-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=2)
annot <- c(info_all$annotatie)


#annotation subtypes
annot_types <- c(info_all$type)

#match color to annotation
annotcol_type <-c()
for (i in 1:length(annot_types)){
  if (annot_types[i]=="LUAD"){
    annotcol_type <-c(annotcol_type, "#66c2a5")
  }
  if (annot_types[i]=="LUSC"){
    annotcol_type <-c(annotcol_type, "#fc8d62")
  }
  if (annot_types[i]=="SCLC"){
    annotcol_type <-c(annotcol_type, "#8da0cb")
  }
}

#select only case samples
meth_bins_case_all <- meth_bins[match(c(names_case$case),c(rownames(meth_bins))),]
meth_bins_case <- meth_bins_case_all[-25,]

meth_bins_case <-meth_bins[match(c(names_case$case),c(rownames(meth_bins))),]
#select only control samples
meth_bins_control <- meth_bins[match(c(names_control$control),c(rownames(meth_bins))),]


############################


exec.score<- function(samples){
  scores <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    sum <- sum(abs(samples[i,]))/5000
    scores <- c(scores, sum)
  }
  return(scores)
}

scores_case = exec.score(meth_bins_case)
print(scores_case)

scores_control = exec.score(meth_bins_control)
print(scores_control)

avg_control=mean(scores_control)
sd_control=sd(scores_control)
FDR_control_low=avg_control-1.96*sd_control 
FDR_control_high=avg_control+1.96*sd_control

FDR_control_l2 = avg_control-(1.96*sd_control)/2
FDR_control_h2 = avg_control+(1.96*sd_control)/2

par(mfrow=c(1,1))
boxplot( scores_case,scores_control)
# Create data
names <- c(rep("LB case", length(scores_case)) , rep("LB control", length(scores_control)))
value <- c(scores_case, scores_control)
data <- data.frame(names,value)
mylevels <- (levels(data$names))
levelProportions <- summary(data$names)/nrow(data)

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\sensitivity")
png(paste0("sensitivity_analysis_noFDR.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
boxplot(data$value ~ data$names ,
        at = c(1,2),
        cex=1,
        cex.axis = 0.8,
        las = 1,
        col = c("#8da0cb", "#e78ac3"),
        border = "black",
        horizontal = T,
        xlab = "sensitivity score",
        ylab=""
)
# Add data points
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- data[data$names==thislevel, "value"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/4)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.8) 
}
# add threshold line
abline(v = c(FDR_control_low, FDR_control_high), col=c("#e78ac3"), lwd=2)
abline(v = c(FDR_control_l2, FDR_control_h2), col=c("#e78ac3"), lwd=2)
dev.off()

FN_l = 0
for (i in 1:length(scores_case)){
  if (scores_case[i]<FDR_control_low | scores_case[i]>FDR_control_high){
    FN_l <- FN_l+1
  }
}
FN_l <- FN_l/length(scores_case)
print(FN_l*40)
print(FN_l)
40-17
17/40

FN_l = 0
for (i in 1:length(scores_case)){
  if (scores_case[i]<FDR_control_l2| scores_case[i]>FDR_control_h2){
    FN_l <- FN_l+1
  }
}
FN_l <- FN_l/length(scores_case)
print(FN_l*40)
print(FN_l)
40-26
26/40
###################################


#annotation
labelscase <- names_case$type

LUAD_indices <- which(labelscase %in% c("LUAD"))
LC_LUAD <- meth_bins_case[LUAD_indices,]

LUSC_indices <- which(labelscase %in% c("LUSC"))
LC_LUSC <- meth_bins_case[LUSC_indices,]

SCLC_indices <- which(labelscase %in% c("SCLC"))
LC_SCLC <- meth_bins_case[SCLC_indices,]

scores_LUAD = exec.score(LC_LUAD)
print(scores_LUAD)

scores_LUSC = exec.score(LC_LUSC)
print(scores_LUSC)
scores_SCLC = exec.score(LC_SCLC)
print(scores_SCLC)


# Create data
names <- c(rep("LUAD", length(scores_LUAD)) , rep("LUSC", length(scores_LUSC)) , rep("SCLC", length(scores_SCLC)))
value <- c(scores_LUAD, scores_LUSC, scores_SCLC)
data <- data.frame(names,value)
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\sensitivity")
png(paste0("sensitivity_analysis subtypes no FDR.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
# Basic boxplot
boxplot(data$value ~ data$names ,
        
        at = c(1,2,3),
        names = c("LUAD", "LUSC", "SCLC"),
        cex=1,
        cex.axis = 0.7,
        las = 1,
        col = c("#a6d854", "#66c2a5", "#fc8d62"),
        border = "black",
        horizontal = T,
        xlab = "sensitivity score",
        ylab=""
)
# Add data points
mylevels <- (levels(data$names))
levelProportions <- summary(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- data[data$names==thislevel, "value"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/4)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.8) 
}
abline(v = c(FDR_control_low, FDR_control_high), col=c("#e78ac3"), lwd=2)
dev.off()

#######################
#repeat on abs values
exec.score.abs<- function(samples){
  scores <- c()
  for (i in 1:length(samples)){
    set.seed(2020)
    s<-abs(samples[i]-avg_control)
    
    scores <- c(scores, s)
  }
  return(scores)
}

abs_case = exec.score.abs(scores_case)
print(abs_case)

abs_control = exec.score.abs(scores_control)
print(abs_control)

avg_control=mean(abs_control)
sd_control=sd(abs_control)
FDR_control_low=avg_control-1.96*sd_control 
FDR_control_high=avg_control+1.96*sd_control


par(mfrow=c(1,1))
boxplot( abs_case,abs_control)
# Create data
names <- c(rep("LB case", length(abs_case)) , rep("LB control", length(abs_control)))
value <- c(abs_case, abs_control)
data <- data.frame(names,value)
mylevels <- (levels(data$names))
levelProportions <- summary(data$names)/nrow(data)

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\sensitivity")
png(paste0("sensitivity_analysis_abs.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
boxplot(data$value ~ data$names ,
        at = c(1,2),
        cex=1,
        cex.axis = 0.8,
        las = 1,
        col = c("#8da0cb", "#e78ac3"),
        border = "black",
        horizontal = T,
        xlab = "sensitivity score",
        ylab=""
)
# Add data points
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- data[data$names==thislevel, "value"]
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/4)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.8) 
}
# add threshold line
abline(v = c(FDR_control_high), col=c("#e78ac3"), lwd=2)
abline(v = c(FDR_control_l2, FDR_control_h2), col=c("#e78ac3"), lwd=2)
dev.off()



FN_l = 0
for (i in 1:length(scores_case)){
  if (scores_case[i]<FDR_control_l2| scores_case[i]>FDR_control_h2){
    FN_l <- FN_l+1
  }
}
FN_l <- FN_l/length(scores_case)
print(FN_l*40)
print(FN_l)


