#####################
#load Z scores volgens segmenten
############################


library(dplyr)
library(factoextra)

#############################
#PREPROCESSING
#loading ZSCORES in bin files
controls_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
controls_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\liquid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)

cases_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
cases_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
#loading Zscores from SEGMENT files
controls_s_segments <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\", pattern="*_segments.bed", full.names=TRUE, recursive=FALSE)
controls_l_segments <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\liquid control\\", pattern="*_segments.bed", full.names=TRUE, recursive=FALSE)

cases_s_segments <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\", pattern="*_segments.bed", full.names=TRUE, recursive=FALSE)
cases_l_segments <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\", pattern="*_segments.bed", full.names=TRUE, recursive=FALSE)


#liquid control
liquidcontrol_seg <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(i in 1:length(controls_l)){
  inhoud <- read.table(controls_l[i], stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  seg <- read.table(controls_l_segments[i], as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  inhoud$zscore.seg <- NA
  for (j in 1:nrow(seg)){
    s <- which(inhoud$chr == seg$chr[j] & inhoud$start == seg$start[j])
    e <- which(inhoud$chr == seg$chr[j] & inhoud$end == seg$end[j])
    inhoud$zscore.seg[s:e] <- seg$zscore[j]
  }
  zscore_seg <- select(inhoud, 7)
  liquidcontrol_seg <-cbind(liquidcontrol_seg,zscore_seg[1:28761,])
  
  name <- unlist(strsplit(controls_l[i], "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcontrol_seg)<-names

#solid control
solidcontrol_seg <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(i in 1:length(controls_s)){
  inhoud <- read.table(controls_s[i], stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  seg <- read.table(controls_s_segments[i], as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  inhoud$zscore.seg <- NA
  for (j in 1:nrow(seg)){
    s <- which(inhoud$chr == seg$chr[j] & inhoud$start == seg$start[j])
    e <- which(inhoud$chr == seg$chr[j] & inhoud$end == seg$end[j])
    inhoud$zscore.seg[s:e] <- seg$zscore[j]
  }
  zscore_seg <- select(inhoud, 7)
  solidcontrol_seg <-cbind(solidcontrol_seg,zscore_seg[1:28761,])
  
  name <- unlist(strsplit(controls_s[i], "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcontrol_seg)<-names

#liquid case
liquidcase_seg <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(i in 1:length(cases_l)){
  inhoud <- read.table(cases_l[i], stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  seg <- read.table(cases_l_segments[i], as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  inhoud$zscore.seg <- NA
  for (j in 1:nrow(seg)){
    s <- which(inhoud$chr == seg$chr[j] & inhoud$start == seg$start[j])
    e <- which(inhoud$chr == seg$chr[j] & inhoud$end == seg$end[j])
    inhoud$zscore.seg[s:e] <- seg$zscore[j]
  }
  zscore_seg <- select(inhoud, 7)
  liquidcase_seg <-cbind(liquidcase_seg,zscore_seg[1:28761,])
  
  name <- unlist(strsplit(cases_l[i], "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(liquidcase_seg)<-names

#solid case
solidcase_seg <- data.frame(matrix(nrow=28761, ncol=0))
names<- c()
for(i in 1:length(cases_s)){
  inhoud <- read.table(cases_s[i], stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  seg <- read.table(cases_s_segments[i], as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  inhoud$zscore.seg <- NA
  for (j in 1:nrow(seg)){
    s <- which(inhoud$chr == seg$chr[j] & inhoud$start == seg$start[j])
    e <- which(inhoud$chr == seg$chr[j] & inhoud$end == seg$end[j])
    inhoud$zscore.seg[s:e] <- seg$zscore[j]
  }
  zscore_seg <- select(inhoud, 7)
  solidcase_seg <-cbind(solidcase_seg,zscore_seg[1:28761,])
  
  name <- unlist(strsplit(cases_s[i], "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  names <- c(names, name)
}
colnames(solidcase_seg)<-names
###########################################

