#############################
#CLUSTERING sWGS 2/07 --> finale figuren
#########################


###PREPROCESSING###
##################

#1. liquid case
#removing bins that don't have RATIO for every entry + convert factors to numeric
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


##############################################
#CLUSTERING
library(gplots)
#1. case vs control --> Euclidian
my_palette <- colorRampPalette(c("black", "white", "red"))(n = 99)
all_liquid <-rbind(lcase_bins_full, lcontrol_bins_full)
annot <- c(pairinfo$pt.ID,  rownames(lcontrol_bins_full))
#standard heatmap --> uses default Eucl dist
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\clustering\\")
png(paste0("clustering_case_control_eucl_3.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)        # smaller font size

heatmap.2(all_liquid,
          keysize = 0.8,
          key.title = "Log2(ratio)",
          key.xlab = "Log2(ratio)",
          density.info="none", 
          col=my_palette, 
          labRow = annot,
          dendrogram="row",Colv="NA",  trace="none",  
          RowSideColors = c(    # grouping row-variables into different categories
            rep("#8da0cb", 44),   # roze= control, paars=case
            rep("#e78ac3", 60))
          )
par(xpd=TRUE)
legend(0.9,1.05,      # location of the legend on the heatmap plot
       legend = c("control", "case"), # category labels
       col = c("#e78ac3", "#8da0cb"),  # color key
       lty= 1,             # line style
       lwd = 8            # line width
)
dev.off()

#########################################
#2. 3 case subtypes --> Pearson

#get annotation
annot <- c(caseinfo$type[1:44])

#match color to annotation
annotcol <-c()
for (i in 1:length(annot)){
  if (annot[i]=="LUAD"){
    annotcol <-c(annotcol,  "#a6d854") #lichtgroen
  }
  if (annot[i]=="LUSC"){
    annotcol <-c(annotcol, "#66c2a5") #groen
  }
  if (annot[i]=="SCLC"){
    annotcol <-c(annotcol, "#fc8d62") #oranje
  }
}

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\clustering\\")
library(gplots)
#heatmap with Pearson correlation 
rows.cor <- cor(t(lcase_bins_full), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
my_palette <- colorRampPalette(c("black", "white", "red"))(n =99)
# creates a 5 x 5 inch image
png("./clustering_3subtypes_pearson_2.png",    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)       # smaller font size
heatmap.2(lcase_bins_full, 
          col =my_palette, 
          density.info = "none",
          Rowv = as.dendrogram(hclust.row), 
          dendrogram="row",Colv="NA",  trace="none",   
          labRow = annot,
          RowSideColors = annotcol, 
          keysize=0.8,
          key.xlab = "Log2(ratio)",
          key.title = "Log2(ratio)")
par(xpd=TRUE)
legend(0.9,1.05,      # location of the legend on the heatmap plot
       legend = c("LUAD", "LUSC", "SCLC"), # category labels
       col = c("#a6d854", "#66c2a5", "#fc8d62"),  # color key
       lty= 1,             # line style
       lwd = 8            # line width
)

dev.off()

##########################################
#3. solid liquid patient pairs

####remove liquid samples that don't have solid counterpart
#find liquid pairs
library("readxl")
pairinfo <- read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\sWGS\\overzicht stalen.xlsx",sheet=3)
#remove patients without solid counterpart
pairinfo_complete<-na.omit(pairinfo)
#remove liquid samples that are not in pairinfo_complete (+samples without annotation --> some random samples are not mentioned anywhere in annotation file)
liquidcase_onlypairs <- liquidcase[,match(pairinfo_complete$cfDNA,colnames(liquidcase))]

####removing bins that don't have ratio for every entry + convert factors to numeric
liquidcase_onlypairs_x<-na.omit(liquidcase_onlypairs)
liquidcase_onlypairs_x <- sapply(liquidcase_onlypairs_x, function(x) if(is.factor(x)) {
  as.numeric(as.character(x))
} else {
  x
})
liquidcase_onlypairs_tot <- data.frame(liquidcase_onlypairs_x)
###binning samples together --> need smaller bins
library(prospectr)
lcase_onlypairs_bins <- binning(t(liquidcase_onlypairs_tot), 1000)


### combine all bins in one big table for clustering
all_case_pairs <-rbind(lcase_onlypairs_bins, scase_bins_full)

###get patient annotation for figure
patients <-c(pairinfo_complete$pt.ID)
pairinfo_sorted_for_solid <- pairinfo_complete[match(colnames(solidcase_tot), pairinfo_complete$tDNA),]
patients1 <-c(patients, pairinfo_sorted_for_solid$pt.ID)
patients_juist <-gsub("_", " ", patients1)

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\sWGS\\clustering\\")
library(gplots)
my_palette <- colorRampPalette(c("black", "white", "red"))(n = 99)
rows.cor <- cor(t(all_case_pairs), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
# creates a 5 x 5 inch image
png("./clustering_case_pairs_pearson_2.png",    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)       # smaller font size
heatmap.2(all_case_pairs, col =my_palette, 
           density.info = "none",
          Rowv = as.dendrogram(hclust.row), 
          dendrogram="row",Colv="NA",  trace="none",
        
          RowSideColors = c(    # grouping row-variables into different categories
            rep("#8da0cb", 33),   # oranje= conder
            rep("#e78ac3", 33)),
          labRow = patients_juist, 
          keysize=0.8,
          key.xlab = "Log2(ratio)",
          key.title = "Log2(ratio)")
par(xpd=TRUE)
legend(0.12,1.06,      # location of the legend on the heatmap plot
       legend = c("liquid", "solid"), # category labels
       col = c("#8da0cb", "#e78ac3"),  # color key
       lty= 1,             # line style
       lwd = 8            # line width
)
dev.off()
