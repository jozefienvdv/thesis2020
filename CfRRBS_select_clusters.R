
#PREPROCESSING
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\select clusters")


lijst1 <- read.table('clusters_test.txt', sep="\t", header=FALSE)
lijst2 <- read.table('promoter_regions_test.txt', sep="\t", header=FALSE)

selected_clusters <- c()
for (i in 1:dim(lijst1)[1]){
  for (j in 1:dim(lijst2)[1]){
   if (lijst1[i,1]==lijst2[j,1]){
     if (lijst1[i,2]>=lijst2[j,2]) {
       if(lijst1[i,3]<=lijst2[j,3]) {
         nieuwe_cluster <- paste(lijst1[i,1], ":", lijst1[i,2], "-", lijst1[i,3], sep="")
         if (!(nieuwe_cluster %in% selected_clusters)){
           selected_clusters <- c(selected_clusters, nieuwe_cluster)
         }
         
       }
     }
   }
  }
}
print(selected_clusters)

write.table(selected_clusters, file = "selected_clusters.txt", sep = "\t",
            row.names = FALSE)


BiocManager::install("GenomicRanges")
library(GenomicRanges) 

gr <- GRanges(
  seqnames = Rle(c(lijst1$V1)),
  ranges = IRanges(lijst1$V2, end = lijst1$V3)
  )

gr2 <- GRanges(
  seqnames = Rle(c(lijst2$V1)),
  ranges = IRanges(lijst2$V2, end = lijst2$V3),
  strand = Rle(lijst2$V5)
  )

x <- findOverlaps(gr, gr2)
y <- queryHits(x)
z <- unique(y)
nieuwe_clusters <- lijst1[z,]

write.table(nieuwe_clusters, file = "selected_clusters.txt", sep = "\t",
            row.names = FALSE)


gr <- GRanges(
  seqnames = Rle(c(lijst1$V1)),
  ranges = IRanges(lijst1$V2, end = lijst1$V3)
  )
print(seqnames(gr))

gr2 <- GRanges(
  seqnames = Rle(c(lijst2$V1)),
  ranges = IRanges(lijst2$V2, end = lijst2$V3),
  strand = Rle(lijst2$V5)
  )
print(seqnames(gr2))
print("objecten ok")
      
x <- findOverlaps(gr, gr2)
print("find overlap ok")
y <- queryHits(x)
z <- unique(y)
nieuwe_clusters <- lijst1[z,]
print("processing ok")
write.table(nieuwe_clusters, file = "selected_clusters3.txt", sep = "\t",
                  row.names = FALSE)
print("write table ok")

