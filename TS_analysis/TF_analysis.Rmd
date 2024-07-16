---
title: "R Notebook"
output: html_notebook
---
#load packages
```{r} 
library(GenomicRanges)
library(seqinr)
library(dplyr)
library(ggpubr)
library(stringr)
library(ape)
library(BSgenome.Hsapiens.UCSC.hg38)

```
#load files
```{r}
fimo_output <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fimo-9.tsv") 
fimo_output$target <- sapply(strsplit(fimo_output$motif_id, "_"), "[[", 1)
fimo_output <- fimo_output %>% filter(q.value <= 0.05) %>% select("target", "sequence_name", "start", "stop") 
colnames(fimo_output) <- c("target", "chrom", "start", "end")

p_bed <- read.delim("/lab/solexa_page/hannah/220516_mpra/helen_start_sites/fantom.bed", header = FALSE)
```

#overlap promoters w TFs
```{r}
p_bed <- GRanges(seqnames = p_bed$V1, ranges = IRanges(start = p_bed$V2, end = p_bed$V3), metadata = p_bed$V4)

fimo_output <- GRanges(seqnames = fimo_output$chrom, ranges = IRanges(start = fimo_output$start, end = fimo_output$end), metadata = fimo_output$target)

get_overlaps <- GenomicRanges::findOverlaps(p_bed, fimo_output)

overlapping_coordinates <- fimo_output[subjectHits(get_overlaps)]
overlapping_metadata1 <- mcols(overlapping_coordinates)
overlapping_metadata2 <- mcols(p_bed[queryHits(get_overlaps)])

full_merge <- cbind(data.frame(overlapping_metadata2), data.frame(overlapping_metadata1))

full_merge <- cbind(full_merge, overlapping_coordinates)
colnames(full_merge) <- c("genes", "tfs", "seqnames", "start", "end", "width", "strand", "copy_tf")

```

#how many binding sites are  on the promoter 
```{r}

promoters <- read.delim('/lab/solexa_page/hannah/220516_mpra/helen_start_sites/fantom.bed', header = FALSE)
colnames(promoters) <- c('chr', 'start', 'end', 'gene')

promoters_by_twenty <- data.frame("chr" = character(), "start" = numeric(), "end" = numeric(), "gene" = character())
for (rowix in (1:nrow(promoters))){
  row_p <- promoters[rowix,]
  n = 20 #number of segments
  start = row_p$start

  for (num in 1:25){
    end = start + 20 
    new_df = data.frame("chr" = row_p$chr, "start" = start, "end" = end, "gene" = row_p$gene)
    promoters_by_twenty <- rbind(promoters_by_twenty, new_df)
    start = end
}}

fimo_output <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fimo-9.tsv")
fimo_output$target <- sapply(strsplit(fimo_output$motif_id, "_"), "[[", 1)
fimo_output <- fimo_output %>% filter(q.value <= 0.05) %>% select("target", "sequence_name", "start", "stop") #only sig
colnames(fimo_output) <- c("target", "chrom", "start", "end")
fimo_output <- as.data.frame(GRanges(seqnames = fimo_output$chrom, ranges = IRanges(start = fimo_output$start, end = fimo_output$end), metadata = fimo_output$target))

genes <- c()
binders <- c() 

for (rowix in (1:nrow(promoters_by_twenty))){
  
  row_p <- promoters_by_twenty[rowix,]
  start = row_p$start
  end = row_p$end
  gene = row_p$gene
  row_p <- GRanges(seqnames = row_p$chr, ranges = IRanges(start = row_p$start, end = row_p$end), metadata = row_p$gene)
 
  
  fimo_output <- GRanges(seqnames = fimo_output$seqnames, ranges = IRanges(start = fimo_output$start, end = fimo_output$end), metadata = fimo_output$metadata)
  get_overlaps <- GenomicRanges::findOverlaps(row_p, fimo_output)
  
  overlapping_coordinates <- data.frame(fimo_output[subjectHits(get_overlaps)])
  if (nrow(overlapping_coordinates) == 0){
    fimo_output <- as.data.frame(fimo_output)

    next
  }
    
  fimo_output <- as.data.frame(fimo_output)
  
  yes_or_no = 'no'
  for (row_btw in (1:nrow(overlapping_coordinates))){
    print('new_row')
    overlapping_mark <- overlapping_coordinates[row_btw,]
    width_mark <- overlapping_mark$width
      
    if ((overlapping_mark$start >= start) & (overlapping_mark$end <= end)){
      print('first if')
      fimo_output <- anti_join(fimo_output, overlapping_mark)  
      yes_or_no = 'yes'
      next
    }
    if (overlapping_mark$start >= start){ #then it 
      print('second if')
      num_bases_in_region <- end - overlapping_mark$start
      if (num_bases_in_region > (width_mark/3)){ #(width_mark/2)
        print('in here')
        fimo_output <- anti_join(fimo_output, overlapping_mark)  
        yes_or_no = 'yes'
        next 
      }
      next
      
    }
    if (overlapping_mark$end >= start){ #then it 
      print('third if')
      num_bases_in_region <-  overlapping_mark$end - start
      if (num_bases_in_region > (width_mark/3)){
        fimo_output <- anti_join(fimo_output, overlapping_mark)  
        yes_or_no = 'yes'
      }
      next
    }
  }
  
  genes <- c(genes, gene)
  binders <- c(binders, yes_or_no)
  
  
} 


```

#graph binding
```{r}
df_nous <- data.frame('genes' = genes, 'binders' = binders) 

df_nous %<>% group_by(genes) %>% mutate(num_binders = n()) %>% select("genes", "num_binders") %>% unique()

g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/g_p.txt") %>% rename('genes' = gene.x)

m_a_g_p <- merge(df_nous, g_p, by = "genes")

m_a_g_p_wide <- pivot_wider(m_a_g_p, id_cols = 'pair', names_from = 'type_of_gene', values_from = 'num_binders')

stat.test <- wilcox.test(m_a_g_p_wide$X_gene, m_a_g_p_wide$Y_gene, paired = TRUE)
  stat.test$p.value
  
pdf(file=(paste0("/lab/solexa_page/hannah/supp_info/figures/5D_tfs.pdf")), width=2.3, height=2.3, colormodel = "rgb")

m_a_g_p %>% ggplot(aes(x = type_of_gene, y = num_binders)) + 
      geom_boxplot(width = 0.2) + 
      geom_point(stroke = 0) + 
      geom_line(aes(group = pair)) + 
      ylab('number of predicted binding sites') + 
      xlab("") + 
      theme_pubr() + 
     # ggrepel::geom_text_repel(aes(label = genes)) + 
      annotate("text", x='X_gene', y = 28, label = paste("p =", format(stat.test$p.value, digits = 2))) 

dev.off()
```

#motif overlap 
```{r}
g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/g_p.txt") %>% dplyr::rename('genes' = gene.x) %>% filter(pair != 'NLGN4X_NLGN4Y' &
                                                                                                                  pair != 'PRKX_PRKY' & 
                                                                                                                  pair != "TMSB4X_TMSB4Y" &
                                                                                                                  pair != "TXLNG_TXLNGY")
                                                                                                                 
x_genes <- (g_p %>% filter(type_of_gene == "X_gene"))$genes
y_genes <-   (g_p %>% filter(type_of_gene == "Y_gene"))$genes

full_merge_ep <- unique(merge(full_merge, g_p, by = "genes") %>% select("genes", "tfs")) 

df_overlaps <- data.frame('pair' = character(), 'num_overlap' = numeric(), 'total_factors' = numeric(), 'is_homolog' = character()) 

for (xg in x_genes){
  
  for (yg in y_genes){
  
    num_overlap <- (full_merge_ep %>% filter(genes == xg | genes == yg) %>% summarize(num_overlap = sum(tfs[which(genes == xg)] %in% tfs[which(genes == yg)])))$num_overlap
    total_factors <- (full_merge_ep %>% filter(genes == xg | genes == yg) %>% summarize(total_factors = length(unique(tfs))))$total_factors

    pair1 = paste0(xg, "_", yg)
    is_homolog = pair1 %in% g_p$pair
    new_df <- data.frame('pair' = pair1, 'num_overlap' = num_overlap, 'total_factors' = total_factors, 'is_homolog' = is_homolog) 
    df_overlaps <- rbind(df_overlaps, new_df)
    
    for (yg1 in y_genes){
        if (yg == yg1){
          next
        }
          num_overlap <- (full_merge_ep %>% filter(genes == yg1 | genes == yg) %>% summarize(num_overlap = sum(tfs[which(genes == yg1)] %in% tfs[which(genes == yg)])))$num_overlap
          total_factors <- (full_merge_ep %>% filter(genes == yg1 | genes == yg) %>% summarize(total_factors = length(unique(tfs))))$total_factors

          pair1 = paste0(yg1, "_", yg)
          new_df <- data.frame('pair' = pair1, 'num_overlap' = num_overlap, 'total_factors' = total_factors, 'is_homolog' = "y_to_y") 
          df_overlaps <- rbind(df_overlaps, new_df)
    }
}
  
    for (xg1 in x_genes){
      print(xg1)
      if (xg1 == xg){
        next
      }
  
      num_overlap <- (full_merge_ep %>% filter(genes == xg | genes == xg1) %>% summarize(num_overlap = sum(tfs[which(genes == xg)] %in% tfs[which(genes == xg1)])))$num_overlap
      total_factors <- (full_merge_ep %>% filter(genes == xg | genes == xg1) %>% summarize(total_factors = length(unique(tfs))))$total_factors

      pair1 = paste0(xg, "_", xg1)
      #is_homolog = pair1 %in% g_p$pair
      new_df <- data.frame('pair' = pair1, 'num_overlap' = num_overlap, 'total_factors' = total_factors, 'is_homolog' = "x_to_x") 
      df_overlaps <- rbind(df_overlaps, new_df)
  }
  
}

df_overlaps$percent <- df_overlaps$num_overlap / df_overlaps$total_factors *100

my_comparisons <- list( c("FALSE", "TRUE"), c("x_to_x", "FALSE"), c("x_to_x", "TRUE"), c("x_to_x", "y_to_y") )
df_overlaps <- unique(df_overlaps %>% select("num_overlap", "total_factors", "is_homolog", "percent"))

```
#graph motif overlap
```{r}
#Outliers assumption
df_overlaps %>%
  group_by(is_homolog) %>%
  identify_outliers(percent)
# Build the linear model
model  <- lm(percent ~ is_homolog,
             data = df_overlaps)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
shapiro_test(residuals(model))
df_overlaps %>%
  group_by(is_homolog) %>%
  shapiro_test(percent)
ggqqplot(df_overlaps, "percent", ggtheme = theme_bw()) +
  facet_grid( ~ is_homolog)
#Homogeneity of variance assumption
df_overlaps %>% levene_test(percent ~ is_homolog)

res.kruskal <- df_overlaps %>% kruskal_test(percent ~ is_homolog)
pwc <- df_overlaps %>% filter(percent > 0) %>%
  dunn_test(percent ~ is_homolog, p.adjust.method = "holm") 
pwc$p.adj.sci <- format(pwc$p.adj, scientific = FALSE, digits = 3)

pwc <- pwc %>%  add_xy_position(x = "is_homolog",dodge = 1,step.increase = .1)

df_overlaps$is_homolog = factor(df_overlaps$is_homolog, levels=c('TRUE', 'FALSE', 'x_to_x', 'y_to_y')) 


pdf(file=(paste0("/lab/solexa_page/hannah/supp_info/figures/S14_motif_overlaps.pdf")), width=5, height=5, colormodel = "rgb")
ggplot(df_overlaps, aes(x=is_homolog, y=percent)) +
  geom_violin() +
  geom_point() +
#ggrepel::geom_text_repel(aes(label = pair)) +
   stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = c(95, 99, 103, 106)) +
  theme_pubr()
dev.off()

```

#get GC of each 50bp region tested (gc_vs_activity.Rmd)
```{r}
# Load the hg38 genome
genome <- BSgenome.Hsapiens.UCSC.hg38

fib_mpra <-  read.delim("/lab/solexa_page/hannah/supp_info/tables/50bp_bed_test_counts_fullfib.txt") %>% filter(chr != "chr6" & chr != "chr7" & chr != "cmv")
fib_mpra$coords <- paste0(fib_mpra$chr, ":", fib_mpra$start, "-", fib_mpra$end)
genes_mpra <- read.delim("/lab/solexa_page/hannah/supp_info/tables/full_fib_50_1.txt") %>% select(c("coords", "gene")) %>% filter(gene != "USP9Yb" & gene != "PRKX"
                       & gene != "PRKY"
                       & gene != "NLGN4X"
                       & gene != "NLGN4Y"
                       & gene != "TMSB4X"
                       & gene != "TMSB4Y"
                       & gene != "TXLNG"
                       & gene != "TXLNGY" 
                       &  gene != "TSPY1")
full_fib <- merge(fib_mpra, genes_mpra, by = "coords")
  #"/lab/solexa_page/hannah/220516_mpra/full_fib_distance_cropped.txt")

get_gc <- function(chromosome, start_coord, end_coord) {
  # Create a GenomicRanges object
gr <- GRanges(seqnames = chromosome, ranges = IRanges(start = start_coord, end = end_coord))

# Extract the sequence
sequence <<- getSeq(genome, gr)

gc_content <- ((oligonucleotideFrequency(sequence, width = 1)[1,][2] + oligonucleotideFrequency(sequence, width = 1)[1,][3]) / sum(oligonucleotideFrequency(sequence, width = 1)[1,])) *100

return(unname(gc_content))

}

gcs <- c()
for (row in (1:nrow(full_fib))){
  sm_f <- full_fib[row,]
  gc_percent <- get_gc(sm_f$chr,sm_f$start,sm_f$end)
  gcs <- c(gcs, gc_percent)
}

full_fib$gc <- gcs

```


#correlate num motifs w activity 
```{r}

fib_mpra <-  read.delim("/lab/solexa_page/hannah/supp_info/tables/50bp_bed_test_counts_fullfib.txt") %>% filter(chr != "chr6" & chr != "chr7" & chr != "cmv")
fib_mpra$coords <- paste0(fib_mpra$chr, ":", fib_mpra$start, "-", fib_mpra$end)
genes_mpra <- read.delim("/lab/solexa_page/hannah/supp_info/tables/full_fib_50_1.txt") %>% select(c("coords", "gene")) %>% filter(gene != "USP9Yb" & gene != "PRKX"
                       & gene != "PRKY"
                       & gene != "NLGN4X"
                       & gene != "NLGN4Y"
                       & gene != "TMSB4X"
                       & gene != "TMSB4Y"
                       & gene != "TXLNG"
                       & gene != "TXLNGY" 
                       &  gene != "TSPY1")
fib_mpra <- merge(fib_mpra, genes_mpra, by = "coords")

#fib_mpra <- fib_mpra %>% filter(distance <= 0)
fib_mpra$gene_log2 <- paste(fib_mpra$gene, fib_mpra$log2, sep = "_")

fimo_output <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fimo-10_mpra.tsv") 
fimo_output$target <- sapply(strsplit(fimo_output$motif_id, "_"), "[[", 1)
fimo_output <- fimo_output %>% filter(q.value <= 0.05) %>% select("target", "sequence_name", "start", "stop") 
colnames(fimo_output) <- c("target", "chrom", "start", "end")

p_bed <- GRanges(seqnames = fib_mpra$chr, ranges = IRanges(start = fib_mpra$start, end = fib_mpra$end), metadata = fib_mpra$gene_log2)

fimo_output <- GRanges(seqnames = fimo_output$chrom, ranges = IRanges(start = fimo_output$start, end = fimo_output$end), metadata = fimo_output$target)

get_overlaps <- GenomicRanges::findOverlaps(p_bed, fimo_output)

overlapping_coordinates <- fimo_output[subjectHits(get_overlaps)]
overlapping_metadata1 <- mcols(overlapping_coordinates)
overlapping_metadata2 <- mcols(p_bed[queryHits(get_overlaps)])

full_merge <- cbind(data.frame(overlapping_metadata2), data.frame(overlapping_metadata1))

full_merge <- cbind(full_merge, overlapping_coordinates)
colnames(full_merge) <- c("genes", "tfs", "seqnames", "start", "end", "width", "strand", "copy_tf")

full_merge %<>% group_by(genes) %>% mutate(num_tf = (n_distinct(tfs))) #instead of n()

not_in_merge <- fib_mpra %>% filter(!(gene_log2 %in% full_merge$genes))
not_in_merge$num_tf <- 0
not_in_merge <- not_in_merge %>% select(c("gene", "log2", "gene_log2", "num_tf"))


full_merge_a <- unique(full_merge %>% dplyr::select("genes", "num_tf"))
full_merge_a_1 <- str_split(full_merge_a$genes, pattern = "_")
full_merge_a <- cbind(as.data.frame(do.call(rbind, full_merge_a_1)),full_merge_a) 
colnames(full_merge_a) <- c("gene",'log2', 'gene_log2','num_tf')
full_merge_a <- rbind(full_merge_a, not_in_merge)
g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/g_p.txt") 
colnames(g_p) <- c("gene",'pair', 'type_of_gene')

m_a_g_p <- merge(full_merge_a, g_p, by = "gene")
m_a_g_p$log2 <- as.numeric(m_a_g_p$log2)

full_fib$gene_log2 <- paste(full_fib$gene, full_fib$log2, sep = "_") 

full_fib1 <- full_fib %>% select("gene_log2", "gc") #from /lab/solexa_page/hannah/220516_mpra/gc_vs_activity.Rmd

m_a_g_p <- left_join(m_a_g_p, full_fib1, by = "gene_log2")

pdf(file=(paste0("/lab/solexa_page/hannah/supp_info/figures/S15.pdf")), width=5, height=5, colormodel = "rgb")
m_a_g_p %>% ggplot(aes(x = gc  , y = num_tf , color = log2)) + #color = type_of_gene
  #facet_wrap(~gene.x) + 
  geom_point(alpha = 0.5, stroke = 0) + 
  ylab('number of TFs') + 
  xlab("G+C content") + 
  theme_pubr() + 
  #xlim(0,100) + 
  stat_cor(method = "pearson", label.x = 25, label.y = 30) + 
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
  scale_color_gradient(low = "blue", high = "red") 
dev.off()

pdf(file=(paste0("/lab/solexa_page/hannah/supp_info/figures/S15B.pdf")), width=5, height=5, colormodel = "rgb")

m_a_g_p %>% ggplot(aes(x = log2  , y = gc , color = num_tf)) + #color = type_of_gene
  #facet_wrap(~gene.x) + 
  geom_point(alpha = 0.5, stroke = 0) + 
  ylab('G+C content') + 
  xlab("log2(RNA/DNA)") + 
  theme_pubr() + 
  #xlim(0,100) + 
  stat_cor(method = "pearson", label.x = -2, label.y = 80)  + 
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
  scale_color_gradient(low = "blue", high = "red") 
dev.off()

pdf(file=(paste0("/lab/solexa_page/hannah/supp_info/figures/5C_new_gc.pdf")), width=2.5, height=2.5, colormodel = "rgb")

m_a_g_p %>% ggplot(aes(x = gc  , y = log2 , color = type_of_gene)) + #color = type_of_gene
  #facet_wrap(~gene.x) + 
  geom_point(alpha = 0.5, stroke = 0) + 
  xlab('G+C content') + 
  ylab("log2(RNA/DNA)") + 
  theme_pubr() + 
  #xlim(0,100) + 
  stat_cor(method = "pearson", label.y = -2, label.x= 80)  + 
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
#+ 
scale_colour_manual(values = c("#e46915", "#857bc6"))
dev.off()

# correlation_result <- cor.test(m_a_g_p$log2, m_a_g_p$gc, method = "pearson")
# 
# p_value <- correlation_result$p.value

# m_a_g_p %>%
#   ggplot(aes(x = log2, y = gc, color = type_of_gene)) +
#   geom_point(alpha = 0.5, stroke = 0) +
#   ylab('G+C content') +
#   xlab("log2(RNA/DNA)") +
#   theme_pubr() +
#   geom_smooth(method='lm', formula= y ~ x, se = FALSE) +
#   scale_colour_manual(values = c("#e46915", "#857bc6")) +
#   annotate("text", x = -2, y = 80, label = paste("p-value =", format(p_value, scientific = TRUE)))

```
