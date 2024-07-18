
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
getwd(
)
library(utils)

myPath = #PATH TO GITHUB 

g_p <- read.delim(myPath, "/tables/g_p.txt") 


#calculate pairwise distances 
 

GENES = c("EIF1AX", "EIF1AY", "KDM5D" , "KDM5C","UTY", "KDM6A", "ZFY", "ZFX", "DDX3Y" ,"DDX3X", "USP9Y" , "USP9X", "RPS4Y1", "RPS4X") 

xgenes = c( "EIF1AX", "KDM5C", "KDM6A",  "ZFX", "DDX3X", "USP9X",  "RPS4X") 

types = c('exon', 'intron' ,'promoter') 

spec2 = c('chimp', 'gorilla', 'pileatedgibbon', 'orangutan', 'mac', 'marm') 

dist_dna_t <- data.frame(gene = character(), species = character(), species2 = character(), branchlen = numeric(), type = character(), align_length = numeric())

cant_calculate <- data.frame(gene = character(), species = character(), species2 = character(), type = character())

for (gene1 in GENES){
  for (type1 in types){
    for (specy in spec2){
       tryCatch({ 
         alignment = read.dna(paste0(myPath, '/multiz_7sp/', gene1, '/humanmasked.', specy, '_', type1, '_msa_filtered.phy'))
         distances <- dist.dna(alignment, model = "N")
         test_df <- data.frame(gene = gene1, species = 'human', species2 = specy, branchlen = distances[1], type = type1, align_length = length(alignment))
         dist_dna_t <- rbind(dist_dna_t, test_df)
         
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


dist_dna_t_mac <- dist_dna_t  %>% filter(type == "exon" | type == "intron" | type == "promoter") %>% select(c("gene", "species", "species2", "branchlen", "type")) 
                                                                                                                                            
dist_dna_t_mac %<>% pivot_wider(values_from = "branchlen", names_from = 'type')
dist_dna_t_mac$promoter_diverg <- dist_dna_t_mac$promoter #promoter
dist_dna_t_mac$intron_diverg <- dist_dna_t_mac$intron

polymorphism <- read.delim('/lab/solexa_page/hannah/1000genomes/polymorphism_gnomad_V4.txt')[,c(1,2,4)]  #polymorphism_gnomad.txt
polymorphism %<>% pivot_wider(values_from = "sum_v", names_from = 'region')
polymorphism$promoter_polym <- polymorphism$promoter
polymorphism$intron_polym <- polymorphism$intron

merged <- merge(dist_dna_t_mac, polymorphism, by = "gene")


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
g_p <- read.delim(paste0(myPath, "/tables/g_p.txt")) %>% rename('gene' = 'gene.x')

merged <- merge(merged, g_p, by = "gene")

ggplot(merged, aes(x=gene1, y=NI, color = pair)) + 
  geom_boxplot() + 
  geom_point() +  
  geom_hline(yintercept = 1, linetype = "dashed") + 
  theme_pubr() + 
  xlab("")



#plot ORs

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























