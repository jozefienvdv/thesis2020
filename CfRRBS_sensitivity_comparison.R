#pairwise comp CNV en cfRRBS sens
avg <- mean(scores_case)
new <- c()
for (i in  scores_case){
  print(i)
  j = abs(avg-i)
  new <- c(new, j)
  print(length(new))
}
cfRRBS.sens <- new

CNV.sens <- scores_lcase[-c(15, 28,31,41)]

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

setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\sensitivity")
png(paste0("sensitivity_analysis vgln.png"),    # create PNG for the heat map        
    width = 12*400,        # 5 x 300 pixels
    height = 5*400,
    res = 400)        # smaller font size
plot( CNV.sens, cfRRBS.sens , pch=21, cex=1.2, bg=annotcol_type, ylab="methylation abs(mean-score)", xlab="CNV sensitivity score")
legend("topright",      # location of the legend on the heatmap plot
       legend = c("LUAD", "LUSC", "SCLC"), # category labels
       col = c("#a6d854", "#66c2a5", "#fc8d62"),  # color key
       lty= 1,             # line style
       lwd = 8,            # line width
       cex = 1.2
)
dev.off()


cor.test(cfRRBS.sens, CNV.sens, "greater", "pearson")
cor(cfRRBS.sens, CNV.sens,  "pairwise.complete.obs")
plot(cor(cfRRBS.sens, CNV.sens,  "pairwise.complete.obs"))
