---
title: "R Notebook"
output: html_notebook
---

#packages
```{r}
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
getwd(
)
library(utils)
g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") 
```

#calculate pairwise distances 
```{r} 

GENES = c("EIF1AX", "EIF1AY", "KDM5D" , "KDM5C","UTY", "KDM6A", "ZFY", "ZFX", "DDX3Y" ,"DDX3X", "USP9Y" , "USP9X", "RPS4Y1", "RPS4X") 

xgenes = c( "EIF1AX", "KDM5C", "KDM6A",  "ZFX", "DDX3X", "USP9X",  "RPS4X") 

types = c('exon', 'intron' ,'promoter', 'promoter_a') 

spec2 = c('chimp', 'gorilla', 'pileatedgibbon', 'orangutan', 'mac', 'marm') 

dist_dna_t <- data.frame(gene = character(), species = character(), species2 = character(), branchlen = numeric(), type = character(), align_length = numeric())

cant_calculate <- data.frame(gene = character(), species = character(), species2 = character(), type = character())

for (gene1 in GENES){
  for (type1 in types){
    for (specy in spec2){
      #print(paste0(gene1, type1, specy))
       tryCatch({ 
         alignment = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1, '/humanmasked.', specy, '_', type1, '_msa_filtered.phy'))
         # alignment = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1, '/humanmasked.', specy, '_', type1, '_msa.phy'))
      
        distances <- dist.dna(alignment, model = "N")
        test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = type1, align_length = length(alignment))
        dist_dna_t <- rbind(dist_dna_t, test_df)

      if (type1 == 'intron'){ #& !(gene1 %in% xgenes)){
          print('okok')
          print(gene1)
        #  "all_species/" + key[0] + "_intron_msa_GCadj" + key[1] + ".phy", "phylip")
          alignment1 = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1,'/', 'humanmasked.', specy, '_intron_msa_GCadjfor_exon.phy')) #since i know it's an intron
          distances <- dist.dna(alignment1, model = "N")
          test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = 'intron_GC_forexon', align_length = length(alignment1))
          dist_dna_t <- rbind(dist_dna_t, test_df)
          #gene1 = "DDX3X"
         alignment2 = read.dna(paste0('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_7sp/', gene1,'/', 'humanmasked.', specy, '_intron_msa_GCadjfor_promoter.phy')) #since i know it's an intron
          print('no_error')
          distances <- dist.dna(alignment2, model = "N")
          test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = 'intron_GC_forpromoter', align_length = length(alignment2))
          dist_dna_t <- rbind(dist_dna_t, test_df)
        }

      }, error = function(e){
          print('here')
          new_df <- data.frame(gene = gene1, species = 'human', species2 = specy, type = type1)
          cant_calculate <- rbind(cant_calculate,new_df )
      }, warning = function(w){
          new_df <- data.frame(gene = gene1, species = 'human', species2 = specy, type = type1)
          cant_calculate <- rbind(cant_calculate,new_df )
          print(new_df)
          print('warned')
      })
    
    }
  closeAllConnections()
  }
}

print(dist_dna_t)
#dist_dna__raw <- dist_dna_t
# dist_dna_t$align_length1 <- dist_dna_t$align_length /2
# dist_dna_t$ratio <- dist_dna_t$branchlen / dist_dna_t$align_length

dist_dna_raw

```

```{r}

dist_dna_t_mac <- dist_dna_t  %>% filter(type == "exon" | type == "intron" | type == "promoter") %>% select(c("gene", "species", "species2", "branchlen", "type")) #%>% filter(species2 == 'orangutan')
                                                                                                                                            
dist_dna_t_mac %<>% pivot_wider(values_from = "branchlen", names_from = 'type')
dist_dna_t_mac$promoter_diverg <- dist_dna_t_mac$promoter #promoter
dist_dna_t_mac$intron_diverg <- dist_dna_t_mac$intron



polymorphism <- read.delim('/lab/solexa_page/hannah/1000genomes/polymorphism_gnomad_V4.txt')[,c(1,2,4)]  #polymorphism_gnomad.txt
polymorphism %<>% pivot_wider(values_from = "sum_v", names_from = 'region')
polymorphism$promoter_polym <- polymorphism$promoter
polymorphism$intron_polym <- polymorphism$intron

merged <- merge(dist_dna_t_mac, polymorphism, by = "gene")
#merged <- merged %>% filter(promoter_diverg > 10)
```
#merged
{\displaystyle \alpha =1-{\frac {D_{s}P_{n}}{D_{n}P_{s}}}}

```{r}
#substitutions in promoters and introns - get the actual number ! 
merged$cds_div <- merged$exon/ merged$intron_diverg
merged$cds_poly <- merged$CDS/ merged$intron_polym
merged$CDS_NI <- merged$cds_poly / merged$cds_div 

merged$div <- merged$promoter_diverg/merged$intron_diverg
merged$poly <- merged$promoter_polym/merged$intron_polym

merged$NI <- merged$poly / merged$div 
merged$gene1 <-merged$gene

# DoS = D(n)/(D(n) + D(s)) - P(n)/(P(n) + P(s))
merged$DOS <- (merged$promoter_diverg / (merged$promoter_diverg + merged$intron_diverg) )- (merged$promoter_polym/(merged$promoter_polym + merged$intron_polym))

  
merged$alpha <- 1 - ((merged$intron_diverg * merged$promoter_polym) / (merged$promoter_diverg*merged$intron_polym))
g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") %>% rename( 'gene' = 'gene.x')

merged <- merge(merged, g_p, by = "gene")

ggplot(merged, aes(x=gene1, y=NI, color = pair)) + 
  geom_boxplot() + 
  geom_point() +  
  geom_hline(yintercept = 1, linetype = "dashed") + 
  theme_pubr() + 
  xlab("")

ggplot(merged, aes(x=gene1, y=DOS, color = pair)) + 
  geom_boxplot() + 
  geom_point() +  
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_pubr() + 
  xlab("")

ggplot(merged, aes(x=gene1, y=CDS_NI, color = pair)) + 
  geom_boxplot() + 
  geom_point() +  
  geom_hline(yintercept = 1, linetype = "dashed") + 
  theme_pubr() + 
  xlab("")

```
#statistics 
```{r}

#could add up all the polymorphisms synonymous/nonsynonymous

simplified <- merged %>% filter(species2 == "mac") %>% na.omit() %>% group_by(gene) %>% mutate(sum_promoter_diverg = sum(promoter_diverg)) %>% mutate(sum_intron_diverg = sum(intron_diverg)) %>% select("gene", "sum_promoter_diverg", "sum_intron_diverg", "promoter_polym", "intron_polym") %>% distinct()

#TEST WITH CDS
# simplified <- merged %>% filter(species2 == "mac") %>% na.omit() %>% group_by(gene) %>% mutate(sum_promoter_diverg = sum(exon)) %>% mutate(sum_intron_diverg = sum(intron_diverg)) %>% select("gene", "sum_promoter_diverg", "sum_intron_diverg", "CDS", "intron_polym") %>% distinct() 


#add up substitutions from  many branches
results_df <- data.frame()

for (num_row in 1:nrow(simplified)){
  
  one_line <- simplified[num_row,]

  matrix <- matrix(c(one_line$sum_intron_diverg, one_line$sum_promoter_diverg, one_line$intron_polym, one_line$promoter_polym ), nrow = 2,
 	              dimnames =list(c("promoter", "intron"), c("divergence", "polymorphism")))
  
 #  matrix <- matrix(c(one_line$sum_intron_diverg, one_line$sum_promoter_diverg, one_line$intron_polym, one_line$CDS ), nrow = 2,
 # 	              dimnames =list(c("promoter", "intron"), c("divergence", "polymorphism")))


  #adjusted_p_values <- p.adjust(p_values, method = "BH")
  result <- fisher.test(matrix)
  result$estimate
  odds_ratio <- result$estimate
  p_value <- result$p.value
  conf_int_1 <- result$conf.int[1]
  conf_int_2 <- result$conf.int[2]
  
  
  new_df <- data.frame('odds_ratio' = odds_ratio, 'pval' = p_value, 'conf_low' = conf_int_1, 'conf_high' = conf_int_2)

  new_df1 <- cbind(one_line, new_df)
  results_df <- rbind(results_df, new_df1)
  
}

results_df$adj_pval <- p.adjust(results_df$pval, method = "BH")

merge_res <- merge(results_df, g_p, by = "gene")


```
#plot ORs
```{r}
#df$Category <- factor(df$Category, levels = c("A", "B", "C"))
merge_res$gene <- factor(merge_res$gene, levels = c("DDX3X", "DDX3Y", "EIF1AX", "EIF1AY", "KDM5C", "KDM5D", "KDM6A", "UTY", "RPS4X", "RPS4Y1", "USP9X", "USP9Y", "ZFX", "ZFY")) 

merge_res$log2odds <- log2(merge_res$odds_ratio) 
merge_res$conf_low <- log2(merge_res$conf_low)
merge_res$conf_high <- log2(merge_res$conf_high)
#pdf('/lab/solexa_page/hannah/1000genomes/CDS_NI.pdf', width = 5, height = 5)
pdf('/lab/solexa_page/hannah/1000genomes/NI_log2_V4.pdf', width = 5, height = 5)
merge_res %>% 
ggplot(aes(x=gene, y = (log2odds), color = type_of_gene)) + 
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.1) + 
 #ylim(0,6) + 
  theme_pubr() + 
  geom_text(aes(label = ifelse(adj_pval < 0.05, "*", "")), 
            size = 10, vjust = -1) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = 'grey')  
dev.off()
```






















