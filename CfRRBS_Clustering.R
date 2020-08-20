#############################
#CLUSTERING cfRRBS 27/07 --> finale figuren
#########################


###PREPROCESSING###
##################

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")

cluster.light <- readRDS("cluster.light.Rds") 
#verder werken met deze data = light gesmooth (missing data invullen) + geclusterd per CpG eilandje (minder datapunten)

#binning samples together --> need smaller bins
#btw binnen van de niet CpG eiland geclusterde files lukt niet --> te weinig RAM
library(prospectr)
meth_bins <- binning(t(cluster.light), 1000)

#select only case samples
meth_bins_case <- meth_bins[match(c(names_case$case),c(rownames(meth_bins))),]

#annotation case control
library("readxl")
names_case<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=3)
names_control<-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=4)
info_all <-read_excel("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\overzicht stalen cfRRBS.xlsx",sheet=2)
annot <- c(info_all$annotatie)
#match color to annotation
annotcol <-c()
for (i in 1:length(annot)){
  if (annot[i]=="case"){
    annotcol <-c(annotcol, "#8da0cb")
  }
  if (annot[i]=="control"){
    annotcol <-c(annotcol, "#e78ac3")
  }
}

#annotation subtypes
annot_types <- c(info_all$type)


#match color to annotation
annotcol_type <-c()
for (i in 1:length(annot_types)){
  if (annot_types[i]=="LUAD"){
    annotcol_type <-c(annotcol_type, "#a6d854")
  }
  if (annot_types[i]=="LUSC"){
    annotcol_type <-c(annotcol_type, "#66c2a5")
  }
  if (annot_types[i]=="SCLC"){
    annotcol_type <-c(annotcol_type, "#fc8d62")
  }
}

###staal CFD1902728 = row 25 (belangrijk als je analyse zonder dit staal wilt doen)
###Dit staal clustert steeds apart, was ook al outlier in multiQC analyse

##############################################
#CLUSTERING
library(gplots)
#1. case vs control --> Euclidian

my_palette <- colorRampPalette(c("white", "white", "pink", "red"))(n = 99)
all_liquid <-meth_bins[-25,]

#standard heatmap --> uses default Eucl dist
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\cfRRBS\\Clustering\\")
png(paste0("clustering_case_control_eucl.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)        # smaller font size

heatmap.2(all_liquid,

          keysize = 0.8,
          key.title = "Percent methylated",
          key.xlab = "Percent methylated",
          density.info="none", 
          col=my_palette, 
          dendrogram="both",  
          trace="none",
          RowSideColors = annotcol[-25]
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

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\images final\\cfRRBS\\Clustering\\")
library(gplots)
#heatmap with Pearson correlation 
rows.cor <- cor(t(meth_bins_case), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
my_palette <- colorRampPalette(c("white", "white", "pink", "red"))(n = 99)
# creates a 5 x 5 inch image
png("./clustering_3subtypes_pearson.png",    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)       # smaller font size
heatmap.2(meth_bins_case, 
          col =my_palette, 
          trace = "none", 
          density.info = "none",
          Rowv = as.dendrogram(hclust.row), 
          dendrogram="both", 

          labRow = annot_types,
          RowSideColors = annotcol_type, 
          keysize=0.8,
          key.title = "Percent methylated",
          key.xlab = "Percent methylated")
par(xpd=TRUE)
legend(0.9,1.05,      # location of the legend on the heatmap plot
       legend = c("LUAD", "LUSC", "SCLC"), # category labels
       col = c("#a6d854", "#66c2a5", "#fc8d62"),  # color key
       lty= 1,             # line style
       lwd = 8            # line width
)

dev.off()

meth_bins_case_no25 <- meth_bins_case[-25,]
rows.cor <- cor(t(meth_bins_case_no25), use = "pairwise.complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
my_palette <- colorRampPalette(c("white", "white", "pink", "red"))(n = 99)
# creates a 5 x 5 inch image
png("./clustering_3subtypes_pearson no 25.png",    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 10*400,
    res = 400)       # smaller font size
heatmap.2(meth_bins_case_no25, 
          col =my_palette, 
          trace = "none", 
          density.info = "none",
          Rowv = as.dendrogram(hclust.row), 
          dendrogram="both", 
          main = "Case subtypes - Pearson",
          labRow = annot_types[-25],
          RowSideColors = annotcol_type[-25], 
          keysize=0.8,
          key.title = "Percent methylated",
          key.xlab = "Percent methylated")
par(xpd=TRUE)
legend(0.9,1.05,      # location of the legend on the heatmap plot
       legend = c("LUAD", "LUSC", "SCLC"), # category labels
       col = c("#a6d854", "#66c2a5", "#fc8d62"),  # color key
       lty= 1,             # line style
       lwd = 8            # line width
)

dev.off()
##########################################

heatmap.2(meth_bins_case_no25,
          
          keysize = 0.8,
          key.title = "Percent methylated",
          key.xlab = "Percent methylated",
          density.info="none", 
          col=my_palette, 
          dendrogram="both",  
          trace="none",
          RowSideColors = annotcol_type[-25])
x
dev.off()
print(plot(1)) # Basi
