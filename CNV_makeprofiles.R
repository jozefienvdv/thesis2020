###################################
#make profiles
###################################


##create echte data met genomic Range

controls_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
controls_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\liquid control\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)

cases_s <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)
cases_l <- list.files(path="C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\", pattern="*_bins.bed", full.names=TRUE, recursive=FALSE)

path <- "C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\CFD1700135_bins.bed"

exec.makeGrange <- function(path, soort){
  suppressMessages(library("GenomicRanges"))
  library(regioneR)
  library(karyoploteR)
  library(CopyNumberPlots)
  bedfile<- read.table(path, as.is=TRUE, stringsAsFactors=F, na.strings=c(NaN,"NAN"," NAN"), header=T )
  #omit NA
  bedfile_na  <-na.omit(bedfile) #!!!!
  bedfile_na <- sapply(bedfile_na, function(x) if(is.factor(x)) {
    as.numeric(as.character(x))
  } else {
    x
  })
  #functie om factors naar numeric JUIST om te zetten
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  #create Grange
  lijst1 <- data.frame(bedfile_na)
  if (soort=="zscore"){
    gr <- GRanges(
      seqnames = Rle(c(paste0("chr", as.character(lijst1$chr)))),
      ranges = IRanges(as.numeric.factor(lijst1$start), end = as.numeric.factor(lijst1$end)),
      lrr = as.numeric.factor(lijst1$zscore),
      cn = as.numeric.factor(lijst1$zscore))
  } else{
    gr <- GRanges(
      seqnames = Rle(c(paste0("chr", as.character(lijst1$chr)))),
      ranges = IRanges(as.numeric.factor(lijst1$start), end = as.numeric.factor(lijst1$end)),
      lrr = as.numeric.factor(lijst1$ratio),
      cn = as.numeric.factor(lijst1$ratio),
      y = as.numeric.factor(lijst1$ratio)
    )
  }
  
  return(gr)
}

#uitvoering robuust
exec.CNVplot.robuust <- function(gr){
  kp2 <- plotKaryotype("hg19", plot.type = 4, main="plot robuust", cex=0.5)
  plotLRR(kp2, snps=gr, ymin=-1, ymax=1, add.axis = FALSE, labels=NA)
  plotCopyNumberCallsAsLines(kp2, cn.calls = gr, style = "segments", ymin=-1, ymax=1, labels="z-score")
}
#uitvoering mooi
exec.CNVplot.mooi <- function(gr, name, col, soort){
  if (soort=="zscore"){
    kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main=paste0("zscore ", name), cex=3)
    kpAddChromosomeNames(kp, srt=45, cex=2)
    plotLRR(kp, gr, ymin=-5, ymax=5, labels = NA,  points.col = col, line.at.0 = FALSE, points.cex = 1, add.axis = FALSE)
    #plotCopyNumberCallsAsLines(kp, cn.calls = gr, ymin=-5, ymax=5, lwd=0.3, add.axis=FALSE, labels = NA)
    kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
    kpAxis(kp, ymin = -5, ymax=5, tick.pos = -4:4, cex=2)
    kpAddLabels(kp, labels = "z-score", cex=3, srt=90, pos=3, label.margin = 0.025)
  } else{
    kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main=name, cex=3)
    kpAddChromosomeNames(kp, srt=45, cex=2)
    plotLRR(kp, gr, ymin=-2, ymax=2, labels = NA,  points.col = col, line.at.0 = FALSE, points.cex = 1, add.axis = FALSE)
    #plotCopyNumberCallsAsLines(kp, cn.calls = gr, ymin=-5, ymax=5, lwd=0.3, add.axis=FALSE, labels = NA)
    kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
    kpAxis(kp, ymin = -2, ymax=2, tick.pos = -4:4, cex=2)
    kpAddLabels(kp, labels = "log2(ratio)", cex=3, srt=90, pos=3, label.margin = 0.025)
  }
}


#uitvoering functions
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\profiles\\")
gr=exec.makeGrange("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\control_data\\solid control\\D1907413_bins.bed", "ratio")
par(mfrow=c(1,1))
exec.CNVplot.robuust(gr)
exec.CNVplot.mooi(gr, "D1907413", "#CCE5FF", "ratio")
png("./D1907413_profile.png",    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400,            # 300 pixels per inch
    pointsize = 5)        # smaller font size

exec.CNVplot.mooi(gr, "D1907413", "#AAAAAA")
dev.off()

create.profile.image<-function(path, name, col, soort="ratio"){
  gr=exec.makeGrange(path, soort)
  exec.CNVplot.mooi(gr, name, col, soort)
  png(paste0("./",soort, "_", name,"_profile.png"),    # create PNG for the heat map        
      width = 12*400,        # 5 x 300 pixels
      height = 5*400,
      res = 400,            # 300 pixels per inch
      pointsize = 5)        # smaller font size
  
  exec.CNVplot.mooi(gr, name, col, soort)
  dev.off()
}
############################################################
#finale figuren maken voor elk staal
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\profiles\\round_2")
for(path in controls_s){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,"#66c2a5") #groen
}
for(path in controls_l){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#fc8d62") #oranje
}
for(path in cases_s){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#e78ac3") #rozepaars
}
for(path in cases_l){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#8da0cb") #blauwpaars
}


#herhalen voor zscores
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\profiles\\profiles zscores")
for(path in controls_s){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,"#66c2a5", "zscore") #paarsblauw pastel
}
for(path in controls_l){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#CCE5FF", "zscore") #blauw pastel
}
for(path in cases_s){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#CCFF99", "zscore") #groen pastel
}
for(path in cases_l){
  name <- unlist(strsplit(path, "[\\]"))[10]
  name <- unlist(strsplit(name, "[_]"))[1]
  create.profile.image(path, name,  "#FFFF99", "zscore") #geel pastel
}

#combine profiles of patients
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\profiles\\")
gr_LB=exec.makeGrange("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\CFD1702661_bins.bed", "ratio")
par(mfrow=c(2,1))
exec.CNVplot.mooi(gr, "CFD1702661", "#CCE5FF", "ratio")
exec.CNVplot.mooi(gr, "CFD1702661", "#CCE5FF", "ratio")

exec.CNVplot.dubbel <-function(gr1, gr2){
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 0
  pp$bottommargin <- 20
  kp <- plotKaryotype(genome="hg38", plot.type = 3, ideogram.plotter = NULL,
                      labels.plotter = NULL, plot.params = pp, cex =3)
  kpAddChromosomeNames(kp, srt=45, cex=2)
  kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
  kpAxis(kp, ymin = -1, ymax=1, cex=2, data.panel = 1)
  kpAddLabels(kp, labels = "log2(ratio)", cex=2, srt=90, pos=3, label.margin = 0.025)
  kpAxis(kp, ymin = 1, ymax=-1, cex=2, data.panel = 2)
  kpAddLabels(kp, labels = "log2(ratio)", cex=2, srt=90, pos=3, label.margin = 0.025, data.panel = 2)
  kpPoints(kp, gr1, pch=".", data.panel = 1, ymin = -1, ymax=1, col="#8da0cb", cex=2 ) #staal 1
  kpPoints(kp, gr2, pch=".", data.panel = 2, ymin = 1, ymax=-1, col="#e78ac3",  cex=2) #staal 2 
}


gr1 <- exec.makeGrange("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\CFD1702661_bins.bed", "ratio")
gr2 <- exec.makeGrange("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\CFD1702661_bins.bed", "ratio")
exec.CNVplot.dubbel(gr1, gr2)

create.profile.dubbel<-function(path1, path2, soort="ratio"){
  gr1=exec.makeGrange(path1, soort)
  gr2=exec.makeGrange(path2, soort)
  exec.CNVplot.dubbel(gr1, gr2)
  png(paste0("./",soort, "_patient21_LUAD","_dubbelprofile.png"),    # create PNG for the heat map        
      width = 12*400,        # 5 x 300 pixels
      height = 5*400,
      res = 400,            # 300 pixels per inch
      pointsize = 5)        # smaller font size
  exec.CNVplot.dubbel(gr1, gr2)
  dev.off()
}

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\profiles\\")

create.profile.dubbel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\liquid case\\CFD1702697_bins.bed", "C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\case_data\\solid case\\D1809461_bins.bed")
