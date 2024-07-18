
#notes: 
#need to run on R4.1.0 

#load files 

library(GenomicRanges)
library(dplyr)
library(tidyr)

gtf_file <- read.delim(paste0(myPath, "/tables/XY_gtf_w_introns_012924.txt")

gtf_file_promoters <- read.delim(paste0(myPath, '/tables/fantom.bed'), header = FALSE)

#bed file
genes <- c("DDX3X", "DDX3Y", "EIF1AX", "EIF1AY", "KDM5C", "KDM5D", "KDM6A", "RPS4X", "RPS4Y1", "USP9X", "USP9Y", "UTY", "ZFX", "ZFY")
types <- c("intron", 'CDS')

for (gene in genes){
  for (type in types){

    bed1 <- unique(gtf_file %>% filter(gene_name == gene & feature ==type) %>%  select('seqname', 'start', 'end'))
    coords2 <- GRanges(seqnames = bed1$seqname, ranges = IRanges(start = bed1$start, end = bed1$end))
    reduce(coords2)
    bed1 <- as.data.frame(reduce(coords2))[,1:3]
    
    write.table(x = bed1, file = paste0(myPath, "/bedfiles/", gene, "_", type, "_0328.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  }
}

for (gene in genes){

  bed2 <- gtf_file_promoters %>% filter(V4 == gene) %>% select("V1", "V2", "V3")
  write.table(x = bed2, file = paste0(myPath, "/bedfiles/", gene, "_", "promoter", "_0328.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}



