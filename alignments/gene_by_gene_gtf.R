#packages

library(ggtree)
library(rtracklayer)
library(GenomicRanges)
library(seqinr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(ape) 
library(magrittr)
getwd()
library(utils)

g_p <- read.delim(paste0(myPath, "/tables/g_p.txt") 


#can delete

gtf <- read.delim(paste0(myPath, "/tables/XY_gtf_w_introns_012924.txt")) 

gtf <- gtf %>% distinct()
write.table(x = gtf, file = paste0(myPath, "/tables/gtf.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


#make gff files 
 
list_of_GENES <- c("DDX3X", "DDX3Y", "EIF1AX", "EIF1AY", "KDM5C", "KDM5D", "KDM6A", "RPS4X", "RPS4Y1", "USP9X", "USP9Y", "UTY", "ZFX", "ZFY")
list_of_start_coords <- c(41329348, 12900108, 20145838, 20571776, 53229207, 19748939, 44869175, 
                          72281248, 2837602, 
                          41081445, 12658368, 13484673, 24145173, 2931281) #these are the starts in the DNA files

GENE = "DDX3X"
ix = 1
for (ix in 1:length(list_of_GENES)){

GENE = list_of_GENES[ix]
start_coord = list_of_start_coords[ix] 

g_p <- read.delim(paste0(myPath, "/tables/g_p.txt"))

#getting the promoter as defined by coordinates   
promoters <- read.delim(paste0(myPath, "/tables/fantom_promoters_1102.txt")) 

promoters <- promoters %>% filter(gene == GENE)
# XY_gtf_w_introns_100923.txt
gtf <- read.delim(paste0(myPath, "/tables/XY_gtf_w_introns_012924.txt")) %>% filter(gene_name == GENE) #%>% filter(feature == "exon")
colnames(gtf) <- c("seqname", "feature","start","end", "gene_name")

for (row in range(1:nrow(promoters))){
  small_df <- promoters[row,]
  gtf <- rbind(gtf, c(small_df$chr, small_df$X, small_df$promoter_start, small_df$promoter_end, GENE))
  }

#get 'new promoter'
gtf <- gtf[,1:5]

g_p_gene <- g_p %>% filter(gene.x == GENE)

gtf$col2 <- '.'
gtf$start <- as.numeric(gtf$start)
gtf$end <- as.numeric(gtf$end)

if (g_p_gene$direction == "plus"){
gtf$start <- gtf$start - start_coord
gtf$end <- gtf$end - start_coord
}
if (g_p_gene$direction == "minus"){
  gtf$end1 <- start_coord - gtf$start 
  gtf$start1 <- start_coord - gtf$end 
  gtf$end <- gtf$end1
  gtf$start <- gtf$start1
}
gtf$col6 <- "+" 
gtf$seqname <- "human" 
gtf1 <- data.frame(a = gtf$seqname, b = gtf$col2, c=gtf$feature, d = gtf$start, e = gtf$end, f = gtf$col2, g = gtf$col6, h = gtf$col2, i = gtf$gene_name) 
gtf1$diff <- gtf1$e - gtf1$d
coords2 <- GRanges(seqnames = gtf1$c, ranges = IRanges(start = gtf1$d, end = gtf1$e), metadata = gtf1$i)

coords2 <- as.data.frame(reduce(coords2)) #GETS THE OVERLAPPING REGIONS 
  
gtf1 <- data.frame(a = "human", b = '.' , c=as.character(coords2$seqnames), d = coords2$start, e = coords2$end, f = '.', g = '+', h = '.', i = GENE) 
class(gtf1$c)

gtf1$c <- ifelse(gtf1$c == "CDS", "exon", gtf1$c)
write.table(x = gtf1, file = paste0(myPath, "/tables/", GENE, "_gtf_all103023.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}






