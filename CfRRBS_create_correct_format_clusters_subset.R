#####################################
#create table with selected clusters
######################################


#PREPROCESSING
#select all clusters
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\Lennart R Pipeline\\output pipeline")
cluster.light <- readRDS("cluster.light.Rds") 

#select subset
setwd("C:\\Users\\Jozefien\\Documents\\2019-2020\\thesis\\cfRRBS\\select clusters\\output - selected clusters")
subset_clusters <-read.delim("selected_clusters5.txt")
#transform subset to right format
subset_juist <-c()
for (i in 1:dim(subset_clusters)[1]){
  entry <- paste(subset_clusters[i,1], ":", subset_clusters[i,2], "-", subset_clusters[i,3], sep="")
  subset_juist <- c(subset_juist, entry)
}
#create table with subset
cluster.light.subset <- cluster.light[match(c(subset_juist),c(rownames(cluster.light))),]
#transform table
cluster.light.subset <- data.frame(t(cluster.light.subset))


subset_clusters <-read.delim("C:\\Users\\Jozefien\\Downloads\\selected_clusters_cancergenes.txt")
#transform subset to right format
subset_juist <-c()
for (i in 1:dim(subset_clusters)[1]){
  entry <- paste(subset_clusters[i,1], ":", subset_clusters[i,2], "-", subset_clusters[i,3], sep="")
  subset_juist <- c(subset_juist, entry)
}
#create table with subset
cluster.light.subsetcancer <- cluster.light[match(c(subset_juist),c(rownames(cluster.light))),]
#transform table
cluster.light.subsetcancer <- data.frame(t(cluster.light.subsetcancer))

