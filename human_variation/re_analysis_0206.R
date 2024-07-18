

library(ggplot2)
library(ggpubr)

myPath = #PATH TO GITHUB

#gnomad V3

g_p <- read.delim(paste0(myPath, "/tables/g_p.txt")) %>% rename(gene = gene.x)
all_variants <- read.delim('variants_gnomad.csv',  sep = ",") #variants_gnomad.csv variants_0328_prom_norm_gc.csv
all_variants1 <- read.delim('variants_gnomad_CDS.csv',  sep = ",") 
all_variants <- rbind(all_variants, all_variants1)
all_variants <- all_variants %>% filter(AC > 0 & qual == "[]" & VARIANT_CLASS == "SNV") #allele count greater than 1 & VARIANT_CLASS == "SNV" & AF < 0.5 & AF > 0.05 & AF < 0.95 

all_var_tolook <- merge(all_variants, g_p, by = "gene")
all_var_tolook$X <- NULL
all_var_tolook <- all_var_tolook %>% select(c("gene", "chrom", "position", "region", "sum_len", "pair", "type_of_gene", "direction")) %>% distinct()

all_var_1 <- all_var_tolook %>% group_by(gene, region) %>% mutate(sum_v = n()) 

all_var_1$ratio <- all_var_1$sum_v / all_var_1$sum_len

all_var_1 <- all_var_1 %>% select("gene", "region", "ratio", "sum_v") %>% distinct() %>% filter(region != "intron_GC") #ned something there 
write.table(x = all_var_1, file = paste0(myPath,'/tables/polymorphism_gnomad.txt'), quote = FALSE, sep = "\t",  row.names = FALSE)


#gnomad V4

g_p <- read.delim(paste0(myPath, "/tables/g_p.txt")) %>% rename(gene = gene.x)
all_variants <- read.delim('variants_gnomad_V4.csv',  sep = ",") #variants_gnomad.csv variants_0328_prom_norm_gc.csv
all_variants <- all_variants %>% filter(AC > 0 & qual == "[]") # & VARIANT_CLASS == "SNV") #allele count greater than 1 & VARIANT_CLASS == "SNV" & AF < 0.5 & AF > 0.05 & AF < 0.95
nrow(all_variants)

all_var_tolook <- merge(all_variants, g_p, by = "gene")
all_var_tolook$X <- NULL
all_var_tolook <- all_var_tolook %>% select(c("gene", "chrom", "position", "region", "sum_len", "pair", "type_of_gene", "direction")) %>% distinct()

all_var_1 <- all_var_tolook %>% group_by(gene, region) %>% mutate(sum_v = n()) 

all_var_1$ratio <- all_var_1$sum_v / all_var_1$sum_len

all_var_1 <- all_var_1 %>% select("gene", "region", "ratio", "sum_v") %>% distinct() %>% filter(region != "intron_GC") #ned something there 
write.table(x = all_var_1, file = paste0(myPath, '/tables/polymorphism_gnomad_V4.txt'), quote = FALSE, sep = "\t",  row.names = FALSE)








