#select clusters in cancr genes
BiocManager::install("biomaRt")

setwd("C:\\Users\\Jozefien\\Downloads")

genes500 <-read.table("gene_set_Custom_Gene_Selection.2020-08-14.tsv", header=T)
genes1000 <-read.table("gene_set_Custom_Gene_Selection_2.2020-08-14.tsv", header=T)

library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
TOP1000_LCgenes<-getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c( 'ensembl_gene_id'),
      values=genes1000,
      mart=ensembl)

lijst1 <- read.table('cluster_pos.txt', sep="\t", header=FALSE)
lijst2 <- data.frame(TOP1000_LCgenes)

suppressMessages(library("GenomicRanges"))
print("library ok")

gr <- GRanges(
  seqnames = Rle(c(lijst1$V1)),
  ranges = IRanges(lijst1$V2, end = lijst1$V3)
)
print(seqnames(gr))
gr2 <- GRanges(
  seqnames = Rle(c(lijst2$chromosome_name)),
  ranges = IRanges(lijst2$start_position, end = lijst2$end_position),
  strand = Rle(lijst2$strand)
)
print(seqnames(gr2))
print("objecten ok")

x <- findOverlaps(gr, gr2)
print("find overlap ok")
y <- queryHits(x)
z <- unique(y)
nieuwe_clusters <- lijst1[z,]
print("processing ok")
write.table(nieuwe_clusters, file = "selected_clusters_cancergenes.txt", sep = "\t",
            row.names = FALSE)
print("write table ok")

