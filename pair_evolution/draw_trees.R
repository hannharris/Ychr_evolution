
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



#calculate pairwise distances 

GENES = c("EIF1AX", "EIF1AY", "KDM5D" , "KDM5C","UTY", "KDM6A", "ZFY", "ZFX", "DDX3Y" ,"DDX3X", "USP9Y" , "USP9X", "RPS4Y1", "RPS4X") 

xgenes = c("EIF1AX", "KDM5C", "KDM6A",  "ZFX", "DDX3X", "USP9X",  "RPS4X") 

types = c('exon', 'intron' ,'promoter') 

spec2 = c('chimp', 'gorilla', 'pileatedgibbon', 'orangutan', 'mac', 'marm') 

dist_dna_t <- data.frame(gene = character(), species = character(), species2 = character(), branchlen = numeric(), type = character())

cant_calculate <- data.frame(gene = character(), species = character(), species2 = character(), type = character())

for (gene1 in GENES){
  for (type1 in types){
    for (specy in spec2){
       tryCatch({ 
         alignment = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1, '/humanmasked.', specy, '_', type1, '_msa_filtered.phy'))
         distances <- dist.dna(alignment, model = "F84")
         test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = type1)
         dist_dna_t <- rbind(dist_dna_t, test_df)

      if (type1 == 'intron'){ #& !(gene1 %in% xgenes)){   
          alignment1 = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1,'/', 'humanmasked.', specy, '_intron_msa_GCadjfor_exon.phy')) #since i know it's an intron
          distances <- dist.dna(alignment1, model = "F84")
          test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = 'intron_GC_forexon')
          dist_dna_t <- rbind(dist_dna_t, test_df)
          alignment2 = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1,'/', 'humanmasked.', specy, '_intron_msa_GCadjfor_promoter.phy')) #since i know it's an intron
          print('no_error')
          distances <- dist.dna(alignment2, model = "F84")
          test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = 'intron_GC_forpromoter')
          dist_dna_t <- rbind(dist_dna_t, test_df)
        }

      }, error = function(e){
          new_df <- data.frame(gene = gene1, species = 'human', species2 = specy, type = type1)
          cant_calculate <- rbind(cant_calculate,new_df )
      }, warning = function(w){
          new_df <- data.frame(gene = gene1, species = 'human', species2 = specy, type = type1)
          cant_calculate <- rbind(cant_calculate,new_df )
      })
    
    }
  closeAllConnections()
  }
}


#write.table(dist_dna_t, paste0(myPath, "/tables/dist_dna_t_1030.txt"), quote = FALSE, row.names = FALSE)


#process dist_dna

dist_dna_t <- read.delim(paste0(myPath, "/tables/dist_dna_t_1030.txt"),sep = " ") 
dist_dna_t %<>% filter(type != 'intron') %<>% filter(type != 'promoter_a')
dist_dna_t %<>% filter(gene != 'NLGN4X') %>% filter(gene != 'NLGN4Y') %>% filter(gene != 'TMSB4X') %>% filter(gene != 'TMSB4Y') %>% filter(gene != 'TXLNG') %>% filter(gene != 'TXLNGY')


get_graphs_normalized_gc_exons(dist_dna_t, g_p = g_p, file_suffix = 'normalized_both_fantom_promoters')


#RATIO BETWEEN 
my_distances <- dist_dna_t #TEST 

g_p$gene <- g_p$gene.x
merged_distances <- merge(my_distances, g_p, by = "gene")

merged_distances_1 <- merged_distances %>% filter(type != "intron_GC_forpromoter") %>% group_by(pair, type, species2) %>% filter(n() == 2) %>% mutate(len_ratio_yvx = (branchlen[which(type_of_gene == "Y_gene")] / branchlen[which(type_of_gene == "X_gene")])) #%>% ungroup()

merged_distances_1 <- unique(merged_distances_1 %>% dplyr::select("pair", "type", "species2", "len_ratio_yvx"))

merged_distances_1 <- merged_distances_1 %>% na.omit() %>% group_by(type) %>% mutate(average_ratio = mean(len_ratio_yvx))
merged_distances_1 <- merged_distances_1 %>% na.omit() %>% group_by(type) %>% mutate(std = sd(len_ratio_yvx)/length(len_ratio_yvx)) 

#shapiro test and wilcox test 
shapiro_test(merged_distances_1$len_ratio_yvx)
merged_distances_1 %>% na.omit() %>% group_by(type) %>% wilcox_test(len_ratio_yvx ~ 1, mu = 1) 

pdf(file=(paste0(myPath, "/S8_y_to_x_branch_len.pdf")), width=2, height=4, colormodel = "rgb") #check over this file 
merged_distances_1 %>% ggplot(aes(x = type, y = len_ratio_yvx)) + 
  geom_boxplot() + 
  #geom_violin() +
  geom_point(position = position_dodge(0.8)) + 
  geom_hline(yintercept=1) + 
  ggpubr::theme_pubr() #+ 
 # ylim(0,2.5)
dev.off()


#calculate alpha

dist_dna_t <- read.delim(paste0(myPath, "/tables/dist_dna_t_1030.txt", sep = " ") 
dist_dna_t %<>% filter(gene != 'NLGN4X') %>% filter(gene != 'NLGN4Y') %>% filter(gene != 'TMSB4X') %>% filter(gene != 'TMSB4Y') %>% filter(gene != 'TXLNG') %>% filter(gene != 'TXLNGY')
my_distances <- dist_dna_t #TEST 

g_p$gene <- g_p$gene.x
merged_distances <- merge(my_distances, g_p, by = "gene")

merged_distances_1 <- merged_distances %>% filter(type != "intron_GC_forpromoter" & type != "intron_GC_forexon") %>% group_by(pair, type, species2) %>% filter(n() == 2) %>% mutate(len_ratio_yvx = (branchlen[which(type_of_gene == "Y_gene")] / branchlen[which(type_of_gene == "X_gene")])) #%>% ungroup()

merged_distances_1 <- unique(merged_distances_1 %>% dplyr::select("pair", "type", "species2", "len_ratio_yvx"))

merged_distances_1 <- merged_distances_1 %>% na.omit() %>% group_by(type) %>% mutate(average_ratio = mean(len_ratio_yvx))
merged_distances_1 <- merged_distances_1 %>% na.omit() %>% group_by(type) %>% mutate(std = sd(len_ratio_yvx)/length(len_ratio_yvx)) #

merged_distances_1$alpha_c <- (2*merged_distances_1$len_ratio_yvx) / (3 - merged_distances_1$len_ratio_yvx)


median_value <- merged_distances_1 %>%
  filter(type == "intron") %>%
  summarise(median_alpha = median(alpha_c)) %>%
  pull()


sd_value <- merged_distances_1 %>%
  filter(type == "intron") %>% filter(alpha_c < 10 & alpha_c > -10) %>%
  summarise(median_alpha = sd(alpha_c)) %>%
  pull()
 
merged_distances_1

#plot
pdf(file=(paste0(myPath, "/S8b_y_to_x_branch_len_alpha.pdf"), width=2, height=3, colormodel = "rgb")

ggplot(merged_distances_1 %>%
  filter(type == "intron"), aes(x = type, y= alpha_c)) + 
  geom_boxplot(width = 0.1) + 
  geom_point(stroke = 0) + 
  ylim(0,10) + 
  #xlim(0,3) + 
 # geom_hline(yintercept = 2.3, linetype = "dashed")  + 
  theme_pubr() 
dev.off() 



#get_graphs_normalized_gc_exons


dist_dna_t<- read.delim(paste0(myPath, "/tables/dist_dna_t_1030.txt",sep = " ")) 
dist_dna_t %<>% filter(type != 'intron') %<>% filter(type != 'promoter_a')
dist_dna_t %<>% filter(gene != 'NLGN4X') %>% filter(gene != 'NLGN4Y') %>% filter(gene != 'TMSB4X') %>% filter(gene != 'TMSB4Y') %>% filter(gene != 'TXLNG') %>% filter(gene != 'TXLNGY')


get_graphs_normalized_gc_exons <- function(dist_dna_t, g_p = g_p, file_suffix = 'GC_new'){

two_colors_sexchr <- c("#e48034", "#7772b5")  
three_colors <- c("#9E64A1","#489A52", "#83b9e3", "brown")
two_colors_ep <- c("#a661a3", "#83b9e3")

my_distances <- dist_dna_t  

g_p$gene <- g_p$gene.x
merged_distances <- merge(my_distances, g_p, by = "gene")


merged_distances$species2 = factor(merged_distances$species2, 
    levels=c("chimp", "gorilla","orangutan", "pileatedgibbon", "mac", "marm"))


pdf(file=(paste(myPath, "/2B_branch_lens_0222.pdf"), width=7, height=5.5, colormodel = "rgb")

print(ggplot(merged_distances %>% filter(type != "intron_GC_forpromoter") %>% mutate(branchlen = branchlen * 100), aes(x = type, y = branchlen, fill = type_of_gene)) +
        #geom_violin() +
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75), stroke = 0) + 
  theme_pubr() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
          facet_wrap(~species2, scales = "free") + 
  ylab("branch length (subst./site)") +
  scale_fill_manual(values = two_colors_sexchr) + 
  #theme(axis.title = element_text(size = 14),
    #axis.text = element_text(size = 14),
    #plot.title = element_text(size = 14))+labs(fill = NULL) +
 # scale_fill_manual(values = xy_colors)  +
   stat_compare_means(
                     method = "wilcox.test", 
                     label = "p.format") ) 
dev.off()


#calculating constraint
my_distances2 <- merged_distances %>% group_by(gene,species2) %>% filter(n() == 4) %>% mutate(promoternorm = branchlen[which(type == "promoter")] / branchlen[which(type == "intron_GC_forpromoter")]) %>% mutate(exonnorm = branchlen[which(type == "exon")] / branchlen[which(type == "intron_GC_forexon")]) 


pdf(file=(paste0(myPath, "/3C_correlationplot.pdf")), width=6, height=3, colormodel = "rgb") 

my_distances2 %>% select("pair", "type_of_gene", "exonnorm", "promoternorm") %>% distinct %>% filter(promoternorm < 3.25 & exonnorm < 0.5) %>%
  
  ggplot(aes(x=exonnorm, y=promoternorm)) + 
  geom_point(stroke = 0) + #aes(color = pair
  facet_wrap(~type_of_gene, scales = "free") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  stat_cor(method = "spearman") + #label.x = xmin, label.y = ymax
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
  theme_pubr()
dev.off()

#to filter: 
dist_norm<- unique(my_distances2 %>% dplyr::filter(type != "intron") %>% dplyr::select("gene", "species2", "type_of_gene", "promoternorm", "exonnorm",  "pair"))
dist_norm <- tidyr::pivot_longer(dist_norm, cols = c(4,5), names_to = 'type')

dist_norm_nums <- dist_norm %>% na.omit() %>% group_by(type_of_gene, type) %>% mutate(mean_val = median(value)) %>% ungroup()

#calculate constraint ratio 
dist_norm_nums <- dist_norm_nums %>% na.omit() %>% group_by(species2, pair, type) %>% filter(n() == 2) %>% mutate(constr_ratio = value[which(type_of_gene == "Y_gene")] / value[which(type_of_gene == "X_gene")]) %>% ungroup() %>% group_by(type_of_gene, type) %>% mutate(med_constr_ratio = median(constr_ratio)) 

pdf(file=(paste0(myPath, "/3B_cds_promoter.pdf")), width=2, height=3, colormodel = "rgb")

dist_norm_nums %>% filter(constr_ratio < 5.5) %>% ggplot(aes(x=type, y = constr_ratio)) + 
  geom_boxplot(width = 0.2) + 
  geom_point(stroke = 0) + 
  stat_compare_means(
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = 5) + 
  scale_fill_manual(values = two_colors_ep) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
  theme_pubr() #+ 
  #ylim(0,5.25)
dev.off()

pdf(file=(paste0(myPath, "/3A_ratios.pdf")), width=3.5, height=3, colormodel = "rgb")
print(dist_norm %>% filter(value < 2.1) %>% 
  ggplot(aes(x=type, y = value, fill = type_of_gene)) + #
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.2) + 
  geom_point(position = position_dodge(width = 0.9), stroke = 0) + 

   scale_fill_manual(values = two_colors_sexchr) + 
  stat_compare_means(
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = 2) +
   scale_y_continuous(expand = expansion(mult = c(0, 0))) +  
  theme_pubr()) 
dev.off()

}


#stats on dist_norm


#Outliers assumption
dist_norm %>%
  group_by(type_of_gene, type) %>%
  identify_outliers(value)
# Build the linear model
model  <- lm(value ~ type_of_gene*type,
             data = dist_norm)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
shapiro_test(residuals(model))
dist_norm %>%
  group_by(type_of_gene, type) %>%
  shapiro_test(value)
ggqqplot(dist_norm, "value", ggtheme = theme_bw()) +
  facet_grid(type_of_gene ~ type)
dist_norm$name <- paste(dist_norm$type_of_gene, dist_norm$type)
res.kruskal <- dist_norm %>% group_by('name') %>% kruskal_test(value ~ name)

dist_norm %>% group_by('name') %>% kruskal_effsize(value ~ name)
pwc <- dist_norm %>% group_by('name') %>% dunn_test(value ~ name, p.adjust.method = "holm") 
pwc$p.adj.sci <- format(pwc$p.adj, scientific = FALSE, digits = 3)

pwc <- pwc %>%  add_xy_position(x = "region",dodge = 1,step.increase = .1)


#random introns

GENES = c("USP9Y") 

spec2 = c('chimp', 'gorilla', 'pileatedgibbon', 'orangutan',  'mac') 
lengths = c(100,200,300,500,1000)
dist_dna_t <- data.frame(gene = character(), species = character(), species2 = character(), branchlen = numeric(), number_rep = numeric(), length_r = numeric()) 

for (length in lengths){

  for (specy in spec2){

  num = 100
  while(num > 0){
        print(num)
        print(specy)
        #specy = "chimp"
        alignment = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/random_samples/USP9Y', specy, num, "_", length, '.phy'))

        distances <- dist.dna(alignment, model = "F84")
        test_df <- data.frame(gene = 'USP9Y_intron' , species = 'human', species2 = specy, branchlen = distances[1], number_rep = num, length_r = length)
        dist_dna_t <- rbind(dist_dna_t, test_df)
        num = num - 1 
  }
   closeAllConnections()
  }
}

write.table(x = dist_dna_t, file = paste0(myPath, "/tables/usp9y_introns.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

dist_dna_t$length_r <- as.character(dist_dna_t$length_r)

dist_dna_t1 <- dist_dna_t %>% group_by(species2,length_r) %>% mutate(std = sd(branchlen)) %>% mutate(average = mean(branchlen)) %>% select("gene", "species2", "std", "average", "length_r")
dist_dna_t1 <- unique(dist_dna_t1)

dist_dna_t1$species2 = factor(dist_dna_t1$species2, 
    levels=c("chimp", "gorilla","orangutan", "pileatedgibbon", "mac", "marm"))

dist_dna_t1$length_r = factor(dist_dna_t1$length_r, 
    levels=c("100", "200", "300", "500","1000"))

pdf(file=(paste0(myPath, "/S7_allspecies.pdf")), width=4, height=4, colormodel = "rgb")
ggplot(dist_dna_t1, aes(x = species2, y = average, color = length_r)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = average - std, ymax = average + std), width = 0.2, position = position_dodge(width = 0.5)) +
  theme_pubr() + 
  ylab("subst./site")
dev.off()



#get GC and alignment perc. through alignment


GENES = c("USP9Y") 
spec2 = c('chimp', 'gorilla', 'pileatedgibbon', 'orangutan',  'mac', 'marm') 
gc_dist <- data.frame("GC" = numeric(), "dist" = numeric(), "spec" = character())

for (specy in spec2){
alignment = read.dna(paste0(myPath, '/multiz_7sp/USP9Y/humanmasked.', specy, '_intron_msa_filtered', '.phy'))

intervals <- seq(from = 1, to = ncol(alignment) - 1000, by = 1000)
for (num in intervals){
  start = num 
  end = num + 1000 
  new_a <- alignment[,start:end]
  distances <- dist.dna(new_a, model = "F84")[1]
  gc_cont <- seqinr::GC(as.character(alignment["humanmaske",][start:end]))
  new_df <- data.frame("GC" = gc_cont, "dist" = distances, "spec" = specy)
  gc_dist <- rbind(gc_dist, new_df)
}}

write.table(x = gc_dist, file = paste0(myPath, "/tables/usp9y_GCdist_introns.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

pdf(file=(paste0(myPath, "/S17_gcvssubst.pdf")), width=7, height=6, colormodel = "rgb")

ggplot(gc_dist, aes(x = GC, y = dist)) + 
  facet_wrap(~spec, scales = "free") + 
  geom_point() +
  stat_cor(method = "spearman") + #label.x = 0, label.y = 0.1) + 
  theme_pubr()
dev.off()





