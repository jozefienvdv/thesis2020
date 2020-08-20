library(dplyr)
library(factoextra)

#############################
#PREPROCESSING
#loading ratios + naming rows
controls_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
controls_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\liquid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)

cases_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
cases_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)


#liquid control
liquidcontrol <- data.frame(matrix(nrow=28761, ncol=0)) #maar tot 28761 --> skip X/Y chrom --> creeert te veel bias
names<- c()
for(path in controls_l){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 5)
  liquidcontrol <-cbind(liquidcontrol,ratio[2:28762,]) 
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcontrol)<-names

#solid control
solidcontrol <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in controls_s){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 5)
  solidcontrol <-cbind(solidcontrol,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcontrol)<-names

#liquid case
liquidcase <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in cases_l){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 5)
  liquidcase <-cbind(liquidcase,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcase)<-names

#solid case
solidcase <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in cases_s){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 5)
  solidcase <-cbind(solidcase,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcase)<-names

###########################################
# load separate liquid case samples


