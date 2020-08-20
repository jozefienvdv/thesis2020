################################
#sensitivity analysis
#################################

#Z scores gebruiken!!
###PREPROCESSING###
##################

#1. liquid case
#removing bins that don't have Z-score for every entry + convert factors to numeric
liquidcase_x<-na.omit(liquidcase) #!!!!
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

#2. solid case
solidcase_x<-na.omit(solidcase) #!!!!
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

#3. liquid control
#removing bins that don't have Z-score for every entry + convert factors to numeric
liquidcontrol_x<-na.omit(liquidcontrol) #!!!!
liquidcontrol_x <- sapply(liquidcontrol_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
liquidcontrol_tot <- data.frame(liquidcontrol_x)
#binning samples together --> need smaller bins
library(prospectr)
lcontrol_bins <- binning(t(liquidcontrol_tot), 1000)
#get annotation of sample subtype
library("readxl")
controlinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=1)
#remove samples without annotation 
lcontrol_bins_full <- lcontrol_bins[match(controlinfo$CNV.ID[1:60],rownames(lcontrol_bins)),]

#4. solid control
#removing bins that don't have Z-score for every entry + convert factors to numeric
solidcontrol_x<-na.omit(solidcontrol) #!!!!
solidcontrol_x <- sapply(solidcontrol_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
solidcontrol_tot <- data.frame(solidcontrol_x)
#binning samples together --> need smaller bins
library(prospectr)
scontrol_bins <- binning(t(solidcontrol_tot), 1000)
#get annotation of sample subtype
library("readxl")
controlinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=1)
#remove samples without annotation
scontrol_bins_full <- scontrol_bins[match(controlinfo$CNV.ID[61:69],rownames(scontrol_bins)),]
#########################################

exec.score<- function(samples){
  scores <- c()
  for (i in 1:dim(samples)[1]){
    set.seed(2020)
    sum <- sum(abs(samples[i,]))/1000
    scores <- c(scores, sum)
  }
  return(scores)
}

scores_lcase = exec.score(lcase_bins_full)
print(scores_lcase)
x<-matrix(scores_lcase)

scores_scase = exec.score(scase_bins_full)
print(scores_scase)

scores_lcontrol = exec.score(lcontrol_bins_full)
print(scores_lcontrol)

scores_scontrol = exec.score(scontrol_bins_full)
print(scores_scontrol)

x <-data.frame(scores_lcontrol)




#gewogen gemiddelde voor groep grootte
avg.control= (mean(scores_lcontrol)+mean(scores_scontrol))/2
#ongewogen gemiddelde voor groep grootte --> veel lager, want veel meer LQB control stalen
avg.control2 = mean(c(scores_lcontrol, scores_scontrol))

SD= sd(c(mean(scores_lcontrol), mean(scores_scontrol)))

avg_scontrol=mean(scores_scontrol)
sd_scontrol=sd(scores_scontrol)
FDR_scontrol=avg_scontrol+1.96*sd_scontrol 

avg_lcontrol=mean(scores_lcontrol)
sd_lcontrol=sd(scores_lcontrol)
FDR_lcontrol=avg_lcontrol+1.96*sd_lcontrol 

#plots score distributrie samen met FDR
par(mfrow=c(1,1))
boxplot( scores_lcase, scores_scase, scores_lcontrol, scores_scontrol)

# Create data
names <- c(rep(" 1 LB case", length(scores_lcase)) , rep(" 2 SB case", length(scores_scase)) , rep("3 LB control", length(scores_lcontrol)), rep("4 SB control", length(scores_scontrol)))
value <- c(scores_lcase, scores_scase, scores_lcontrol, scores_scontrol)
data <- data.frame(names,value)
mylevels <- (levels(data$names))
levelProportions <- summary(data$names)/nrow(data)


setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\sensitivity analysis")
png(paste0("sensitivity_analysis.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
# Basic boxplot
boxplot(data$value ~ data$names ,
        at = c(1,2,3,4),
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
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.8) 
}
# add threshold line
abline(v = c(FDR_lcontrol, FDR_scontrol), col=c("#8da0cb", "#e78ac3"), lwd=2)
dev.off()

#False negative rate = procent case stalen onder threshold
FN_l = 0
for (i in 1:length(scores_lcase)){
  if (scores_lcase[i]<FDR_lcontrol){
    FN_l <- FN_l+1
  }
}
FN_l <- FN_l/length(scores_lcase)
print(FN_l*44)

FN_s = 0
for (i in 1:length(scores_scase)){
  if (scores_scase[i]<FDR_scontrol){
    FN_s <- FN_s+1
  }
}
FN_s <- FN_s/length(scores_scase)
print(FN_s*33)

#################################################################################
#repeat for liquid case
#get annotation of sample subtype
library("readxl")
caseinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=2)
#remove samples without annotation (1 tot 44, want enkel liquid stalen)
lcase_bins_full <- lcase_bins[match(caseinfo$CNV.ID[1:44],rownames(lcase_bins)),]

#annotation
labelscase <- c(caseinfo$type[1:44])

LUAD_indices <- which(labelscase %in% c("LUAD"))
LC_LUAD <- lcase_bins_full[LUAD_indices,]


LUSC_indices <- which(labelscase %in% c("LUSC"))
LC_LUSC <- lcase_bins_full[LUSC_indices,]

SCLC_indices <- which(labelscase %in% c("SCLC"))
LC_SCLC <- lcase_bins_full[SCLC_indices,]

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
# Basic boxplot
boxplot(data$value ~ data$names ,
        main = "Sensitivity analysis liquid case",
        at = c(1,2,3),
        names = c("LUAD", "LUSC", "SCLC"),
        cex=1,
        cex.axis = 0.7,
        las = 1,
        col = c("#8da0cb", "#66c2a5", "#fc8d62"),
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
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.8) 
}


labelscase_solid <- c(caseinfo$type[45:77])

LUAD_indices_solid <- which(labelscase_solid %in% c("LUAD"))
LC_LUAD_solid <- scase_bins_full[LUAD_indices_solid,]

LUSC_indices_solid <- which(labelscase_solid %in% c("LUSC"))
LC_LUSC_solid <- scase_bins_full[LUSC_indices_solid,]

SCLC_indices_solid <- which(labelscase_solid %in% c("SCLC"))
LC_SCLC_solid <- scase_bins_full[SCLC_indices_solid,]

scores_LUAD_solid = exec.score(LC_LUAD_solid)
print(scores_LUAD_solid)
scores_LUSC_solid = exec.score(LC_LUSC_solid)
print(scores_LUSC)
scores_SCLC_solid = exec.score(LC_SCLC_solid)
print(scores_SCLC_solid)


boxplot(scores_LUAD, scores_LUSC, scores_SCLC, scores_LUAD_solid, scores_LUSC_solid, scores_SCLC_solid,
        main = "Sensitivity analysis case subtypes",
        at = c(1,2,3, 8,9,10),
        names = c("LUAD", "LUSC", "SCLC", "LUAD", "LUSC", "SCLC"),
        cex=1,
        cex.axis = 0.7,
        las = 1,
        col = c("#a6d854", "#66c2a5", "#fc8d62"),
        border = "black",
        horizontal = TRUE,
        notch = F
)
abline(v = c(FDR_scontrol, FDR_lcontrol), col=c("black", "grey"))

# Create data
names <- c(rep("LB LUAD", length(scores_LUAD)) , rep("LB LUSC", length(scores_LUSC)) , rep("LB SCLC", length(scores_SCLC)),
           rep("SB LUAD", length(scores_LUAD_solid)) , rep("SB LUSC", length(scores_LUSC_solid)) , rep("SB SCLC", length(scores_SCLC_solid)))
value <- c(scores_LUAD, scores_LUSC, scores_SCLC, scores_LUAD_solid, scores_LUSC_solid, scores_SCLC_solid)
data <- data.frame(names,value)


setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\sensitivity analysis")
png(paste0("sensitivity_analysis_case_subtypes.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
# Basic boxplot
boxplot(data$value ~ data$names ,
        at = c(1,2,3, 4, 5, 6),
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
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(thisvalues, myjitter, pch=20, col="black", cex = 0.6) 
}
abline(v = c(FDR_scontrol, FDR_lcontrol), col=c( "#e78ac3", "#8da0cb"), lwd=2)
dev.off()

#False negative rate = procent case stalen onder threshold
FN_l_LUAD = 0
for (i in 1:length(scores_LUAD)){
  if (scores_LUAD[i]<FDR_lcontrol){
    FN_l_LUAD <- FN_l_LUAD+1
  }
}
FN_l_LUAD <- FN_l_LUAD/length(scores_LUAD)
print(FN_l_LUAD*16)
print(FN_l_LUAD)

FN_l_LUSC = 0
for (i in 1:length(scores_LUSC)){
  if (scores_LUSC[i]<FDR_lcontrol){
    FN_l_LUSC <- FN_l_LUSC+1
  }
}
FN_l_LUSC <- FN_l_LUSC/length(scores_LUSC)
print(FN_l_LUSC*13)
print(FN_l_LUSC)

FN_l_LUAD_solid = 0
for (i in 1:length(scores_LUAD_solid)){
  if (scores_LUAD_solid[i]<FDR_scontrol){
    FN_l_LUAD_solid <- FN_l_LUAD_solid+1
  }
}
FN_l_LUAD_solid <- FN_l_LUAD_solid/length(scores_LUAD_solid)
print(FN_l_LUAD_solid*13)
print(FN_l_LUAD_solid)


FN_l_LUSC_solid = 0
for (i in 1:length(scores_LUSC_solid)){
  if (scores_LUSC_solid[i]<FDR_scontrol){
    FN_l_LUSC_solid <- FN_l_LUSC_solid+1
  }
}
FN_l_LUSC_solid <- FN_l_LUSC_solid/length(scores_LUSC_solid)
print(FN_l_LUSC_solid*13)
print(FN_l_LUSC_solid)
