library(dplyr)
library(factoextra)

#############################
#PREPROCESSING
#loading ZSCORES + naming rows
controls_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
controls_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\liquid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)

cases_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
cases_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)


#liquid control
liquidcontrol_z <- data.frame(matrix(nrow=28761, ncol=0)) #maar tot 28761 --> skip X/Y chrom --> creeert te veel bias
names<- c()
for(path in controls_l){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 6)
  liquidcontrol_z <-cbind(liquidcontrol_z,ratio[2:28762,]) 
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcontrol_z)<-names

#solid control
solidcontrol_z <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in controls_s){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 6)
  solidcontrol_z <-cbind(solidcontrol_z,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcontrol_z)<-names

#liquid case
liquidcase_z <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in cases_l){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 6)
  liquidcase_z <-cbind(liquidcase_z,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcase_z)<-names

#solid case
solidcase_z <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(path in cases_s){
  inhoud <- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN") )
  ratio <- select(inhoud, 6)
  solidcase_z <-cbind(solidcase_z,ratio[2:28762,])
  
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcase_z)<-names

###########################################
# load separate liquid case samples


