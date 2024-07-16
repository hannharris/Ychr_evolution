
#all code for mpra analysis and figures

#load packages

#install.packages("DescTools")
#library(DescTools)

library(DescTools) 
library(ggplot2)
library(ggbreak)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(purrr)


#colors for plots

scalescolors <- c("#66C2A5", "#FC8D62", "#8DA0CB")
xy_colors <- c('#1B9E77','#7570B3') 
expr_vs_prom_colors <- c('#A6CEE3','#1F78B4') 


#load data 

start_sites <- read.delim("/lab/solexa_page/hannah/supp_info/tables/starts.txt")
start_sites <- rbind(start_sites, data.frame("gene" = c("EEF1A1", "ACTB", "CMV"), "start" = c(73521032, 5530709, 1))) #positive ctls 

fib_counts_1 <- read.delim("/lab/solexa_page/hannah/220516_mpra/Count_Basic/work/32/cc2c02d5c079f4b217c5d91737a57e/Fib_7169_1_counts.tsv")
fib_counts_2 <- read.delim("/lab/solexa_page/hannah/220516_mpra/Count_Basic/work/96/bc3fc6b4ba1178ac8933a8a9611c6b/Fib7276_2_counts.tsv") 
lcl_counts <- read.delim("/lab/solexa_page/hannah/220516_mpra/LCL_1_counts.tsv")
lcl_counts <- lcl_counts[,c(1,5,6)]
colnames(lcl_counts) <- c('name', 'lcl_log2', 'lcl_n_obs_bc')

#gene_pairs <- read.delim("/lab/solexa_page/hannah/220516_mpra/g_p.txt")
gene_pairs <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") 
gene_pairs$gene <- gene_pairs$gene.x 
gene_pairs$gene.x <- NULL
gene_pairs <- gene_pairs %>%  filter(gene != "USP9Yb" & gene != "PRKX"
                       & gene != "PRKY"
                       & gene != "NLGN4X"
                       & gene != "NLGN4Y"
                       & gene != "TMSB4X"
                       & gene != "TMSB4Y"
                       & gene != "TXLNG"
                       & gene != "TXLNGY")

metadata <- read.delim("/lab/solexa_page/hannah/220516_mpra/metadata_0119.txt")
metadata <- metadata %>% filter(gene != "USP9Yb" & gene != "PRKX"
                       & gene != "PRKY"
                       & gene != "NLGN4X"
                       & gene != "NLGN4Y") 

td <- full_join(fib_counts_1, fib_counts_2, by = "name")[,c(1,5,6,10, 11)] %>% filter(name != "no_BC")
colnames(td) <- c("name", "fib_1_log2", "fib_1_n_obs_bc", "fib_2_log2", "fib_2_n_obs_bc")

m_td <- merge(td, metadata, by = "name") 
m_td <- full_join(m_td, lcl_counts, by = "name")
m_td$full_fib_log2 <- (m_td$fib_1_log2 + m_td$fib_2_log2) /2

full_fib <- td 
full_fib$n_obs_bc <- full_fib$fib_1_n_obs_bc + full_fib$fib_2_n_obs_bc
full_fib$log2 <- (full_fib$fib_1_log2 + full_fib$fib_2_log2) /2
full_fib <- na.omit(full_fib)
full_fib <- full_fib[,c(1,3,4,5,7,6)]




#make bed for ucsc genome browser

beds <- make_bedgraph(full_fib, 'file_nameUCSC.txt', 'name')

beds <- beds[[2]] %>% filter(chromosome != "cmv")
beds <- separate(beds, c(3), c("A", "b"), sep = "_")
beds <- beds %>% filter(b != "shuff")
beds$b <- NULL
write.table(x = beds, file = "/lab/solexa_page/hannah/supp_info/file_nameUCSC.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



#fig. S10A library composition

df_n_t <- data.frame('type' = c("negative_ctls", "positive_ctls", 'test_seqs'), 'values' = c(240, 103, 1232))

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S10_lib_content.pdf"), width=4, height=4, colormodel = "rgb") 

ggplot(df_n_t, aes(x=type, y=values)) + 
  geom_col(width = 0.2) + 
  theme_pubclean() + 
  ylab("number of sequences") + 
  xlab("")

dev.off()


#correlation between libs 

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S10_corr_libs.pdf"), width=4, height=4, colormodel = "rgb") 
ggplot(m_td, aes(x= fib_1_log2, y= fib_2_log2)) + 
  stat_cor(method = "pearson") + 
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
  geom_point(aes(color = type), alpha = 0.5, stroke = 0, size = 2) + 
  theme_pubr() + 
  geom_text_repel(aes(label = gene)) +  
  scale_color_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) + 
  stat_regline_equation(label.y = 5) 
dev.off()

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/4C_corr_full_fib_lcl.pdf"), width=4, height=4, colormodel = "rgb")
ggplot(m_td %>% filter(type != "positive ctl" & type != "negative ctl" & gene != "USP9Yb" & gene != "PRKX"
                       & gene != "PRKY"
                       & gene != "NLGN4X"
                       & gene != "NLGN4Y" & gene != "TMSB4Y"
                       & gene != "TMSB4X" & gene != "TXLNG"
                       & gene != "TXLNGY" &lcl_n_obs_bc >= 2 & fib_1_n_obs_bc >=2 & fib_2_n_obs_bc >= 2), aes(x= full_fib_log2, y= lcl_log2)) + 
  stat_cor(method = "pearson") + 
  geom_smooth(method='lm', formula= y~x, se= FALSE) + 
  geom_point(aes(color = chromosome), alpha = 0.8, stroke = 0, size = 2) +
  theme_pubr() + 
  #geom_text_repel(aes(label = gene)) +  
  scale_color_manual(values = c("#e46915", "#857bc6")) + 
  stat_regline_equation(label.y = 2.5) +
  xlim(-2,3) +
  ylim(-2,3)
dev.off()

#histogram of n_bcs 
pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S10_rep1_fib_bcs.pdf"), width=4, height=4, colormodel = "rgb") 
ggplot(m_td, aes(x=fib_1_n_obs_bc)) +
  geom_histogram() +
  stat_bin(bins = 50)+
  theme_pubclean()
dev.off()
pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S10_fib_bcs.pdf"), width=4, height=4, colormodel = "rgb") 
ggplot(m_td, aes(x=fib_2_n_obs_bc)) +
geom_histogram() +
  stat_bin(bins = 50) +
  theme_pubclean()
dev.off()
pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S10_lcl_bcs.pdf"), width=4, height=4, colormodel = "rgb") 
ggplot(m_td %>% filter(name!= 'no_BC'), aes(x=lcl_n_obs_bc)) +
geom_histogram() +
  stat_bin(bins = 50) +
  theme_pubclean()
dev.off()

#median(fib_counts_1$n_obs_bc)



#create dataframes 

lcl_counts <- read.delim("/lab/solexa_page/hannah/220516_mpra/LCL_1_counts.tsv")
lcl_50 <- write_50_df("/lab/solexa_page/hannah/supp_info/tables/lcl_50_1.txt", metadata, lcl_counts)
lcl_200 <- write_200_df("/lab/solexa_page/hannah/supp_info/tables/lcl_200.txt", metadata, lcl_counts)

fib_1_50 <- write_50_df("/lab/solexa_page/hannah/supp_info/tables/fib_1_50_1.txt", metadata, fib_counts_1)
fib_2_50 <- write_50_df("/lab/solexa_page/hannah/supp_info/tables/fib_2_50_1.txt", metadata, fib_counts_2)
full_fib_50 <- write_50_df("/lab/solexa_page/hannah/supp_info/tables/full_fib_50_1.txt", metadata, full_fib)
fib_1_200 <- write_200_df("/lab/solexa_page/hannah/supp_info/tables/fib_1_200_1.txt", metadata, fib_counts_1)
fib_2_200 <- write_200_df("/lab/solexa_page/hannah/supp_info/tables/fib_2_200_1.txt", metadata, fib_counts_2)
full_fib_200 <- write_200_df("/lab/solexa_page/hannah/supp_info/tables/full_fib_200_1.txt", metadata, full_fib)


#load dataframes 

lcl_50 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/lcl_50_1.txt")
lcl_200 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/lcl_200.txt")
fib_1_50 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fib_1_50_1.txt")
fib_2_50 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fib_2_50_1.txt")
full_fib_50 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/full_fib_50_1.txt")
fib_1_200 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fib_1_200_1.txt")
fib_2_200 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/fib_2_200_1.txt")
full_fib_200 <- read.delim("/lab/solexa_page/hannah/supp_info/tables/full_fib_200_1.txt")



#get distances 
  for control plots and for graphing controls based on distance to TSS


lcl_50_dist <- get_dist_TSS(lcl_50, start_sites)
lcl_200_dist <- get_dist_TSS(lcl_200, start_sites)
full_fib_50_dist <- get_dist_TSS(full_fib_50, start_sites)
full_fib_200_dist <- get_dist_TSS(full_fib_200, start_sites) 

make_bedgraph(full_fib_50_dist, file_name, column_to_sep)

#graph controls

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/4B_fib_ctls.pdf"), width=4, height=4, colormodel = "rgb")

full_fib_dist_200_1 <- full_fib_200_dist %>% filter(distance < 500 & distance > -500)
res.kruskal <- full_fib_dist_200_1 %>% kruskal_test(log2 ~ type)
full_fib_dist_200_1 %>% kruskal_effsize(log2 ~ type)
pwc <- full_fib_dist_200_1 %>% 
  dunn_test(log2 ~ type, p.adjust.method = "holm") 
pwc$p.adj.sci <- format(pwc$p.adj, scientific = FALSE, digits = 3)

pwc <- pwc %>%  add_xy_position(x = "type",dodge = 1,step.increase = .1)
full_fib_dist_200_1$type = factor(full_fib_dist_200_1$type, levels=c('negative ctl', 'positive ctl', 'test')) 

ggviolin(full_fib_dist_200_1, x = "type", y = "log2", trim = TRUE) + stat_pvalue_manual(pwc,   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 0)) + stat_boxplot( width = 0.1)
dev.off()

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S11_lcl_ctls.pdf"), width=4, height=4, colormodel = "rgb")

lcl_200_dist_1 <- lcl_200_dist %>% filter(distance < 500 & distance > -500)
res.kruskal <- lcl_200_dist_1 %>% kruskal_test(log2 ~ type)
lcl_200_dist_1 %>% kruskal_effsize(log2 ~ type)
pwc <- lcl_200_dist_1 %>% 
  dunn_test(log2 ~ type, p.adjust.method = "holm") 
pwc$p.adj.sci <- format(pwc$p.adj, scientific = FALSE, digits = 3)

pwc <- pwc %>%  add_xy_position(x = "type",dodge = 1,step.increase = .1)
lcl_200_dist_1$type = factor(lcl_200_dist_1$type, levels=c('negative ctl', 'positive ctl', 'test')) 

ggviolin(lcl_200_dist_1, x = "type", y = "log2", trim = TRUE) + stat_pvalue_manual(pwc,   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 0)) + stat_boxplot( width = 0.1)
dev.off()



#graph promoters

# per_pair_dist_graphs(full_fib_200_dist, 500, 'test100fib', gene_pairs)
# per_pair_dist_graphs(lcl_200_dist, 500, 'test100lcl', gene_pairs)

per_pair_dist_graphs(lcl_50_dist, 500, "S12", gene_pairs)
per_pair_dist_graphs(full_fib_50_dist, 500, "5A", gene_pairs)

per_pair_dist_graphs(full_fib_50_dist, 500, "5A_noerror", gene_pairs)




#get and plot AUC 

aucsm_full_fib <- get_AUC(full_fib_200_dist, 500, gene_pairs)
aucsm_lcl <- get_AUC(lcl_200_dist, 500, gene_pairs)


aucsm_full_fib <- get_AUC(full_fib_50_dist, 500, gene_pairs)
aucsm_lcl <- get_AUC(lcl_50_dist, 500, gene_pairs)

write.table(x = aucsm_full_fib, file = "/lab/solexa_page/hannah/supp_info/tables/aucsm_full_fib.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(x = aucsm_lcl, file = "/lab/solexa_page/hannah/supp_info/tables/aucsm_lcl.txt", quote = FALSE, sep = "\t", row.names = FALSE)

plot_AUC(aucsm_full_fib, '/lab/solexa_page/hannah/supp_info/figures/5C_auc_scatter_fib.pdf', 620, 620)
plot_AUC(aucsm_lcl, '/lab/solexa_page/hannah/supp_info/figures/S13_auc_scatter_lcl.pdf', 1000, 1000)



#get and plot max peak

full_fib_50_dist1 <- full_fib_50_dist %>% filter(distance < 100 & distance > -300)
max_full_fib <- max_peak(full_fib_50_dist1)

lcl_50_dist1 <- lcl_50_dist %>% filter(distance < 100 & distance > -300)

max_lcl <- max_peak(lcl_50_dist1)

write.table(x = max_full_fib, file = "/lab/solexa_page/hannah/supp_info/tables/max_full_fib.txt", quote = FALSE, sep = "\t", row.names = FALSE) #redo this
write.table(x = max_lcl, file = "/lab/solexa_page/hannah/supp_info/tables/max_lcl.txt", quote = FALSE, sep = "\t", row.names = FALSE) #redo this

max_full_fib <- merge(max_full_fib, gene_pairs, by = 'gene')
max_lcl <- merge(max_lcl, gene_pairs, by = 'gene')

plot_maxpeak(max_full_fib, '/lab/solexa_page/hannah/supp_info/figures/5_maxpeak_fib.pdf', 2, 2)
plot_maxpeak(max_lcl, '/lab/solexa_page/hannah/supp_info/figures/S13_maxpeak_lcl.pdf', 2.5, 2.5)






#AUC < or > 0 fraction


full_fib_50_dist_5hun <- full_fib_50_dist %>% filter(distance <= 500 & distance >= -500)
full_fib_50_dist_5hun$distunder0 <- ifelse(full_fib_50_dist_5hun$distance <= 0, "yes", "no" )

lcl_50_dist_5hun <- lcl_50_dist %>% filter(distance <= 500 & distance >= -500)
lcl_50_dist_5hun$distunder0 <- ifelse(lcl_50_dist_5hun$distance <= 0, "yes", "no" )


#switch to AUC
full_fib_50_dist_5hun$log2_int <- ifelse(full_fib_50_dist_5hun$log2 > 0, full_fib_50_dist_5hun$log2,  0) 

full_fib_50_dist_5hun <- full_fib_50_dist_5hun %>% filter(!(gene == "ZFX" & distance > 22))

full_fib_50_dist_5hun %<>%
  group_by(gene, distunder0) %>%
  mutate(area = AUC(distance,log2_int, method="trapezoid")) %>% ungroup()

full_fib_50_dist_5hun_g <- (full_fib_50_dist_5hun %>% select(c("gene", "distunder0", "area")) %>% distinct()) %>% filter(gene != "ACTB" & gene != "CMV" & gene != "EEF1A1") 


full_fib_50_dist_5hun_g %<>% pivot_wider(id_cols = 'gene', names_from = 'distunder0', values_from = 'area')

full_fib_50_dist_5hun_g$fraction <- full_fib_50_dist_5hun_g$yes / (full_fib_50_dist_5hun_g$no + full_fib_50_dist_5hun_g$yes)*100

full_fib_50_dist_5hun_g$celltype <- 'Fib'
full_fib_50_dist_5hun_g$up_frac <- full_fib_50_dist_5hun_g$yes / (full_fib_50_dist_5hun_g$no + full_fib_50_dist_5hun_g$yes)*100
full_fib_50_dist_5hun_g[13,6] <- 100 #make 100 bc was NA 

full_fib_50_dist_5hun_g$down_frac <- full_fib_50_dist_5hun_g$no / (full_fib_50_dist_5hun_g$no + full_fib_50_dist_5hun_g$yes)*100

###
lcl_50_dist_5hun$log2_int <- ifelse(lcl_50_dist_5hun$log2 > 0, lcl_50_dist_5hun$log2,  0) 

lcl_50_dist_5hun %<>%
  group_by(gene, distunder0) %>%
  mutate(area = AUC(distance,log2_int, method="trapezoid")) %>% ungroup()

lcl_50_dist_5hun_g <- (lcl_50_dist_5hun %>% select(c("gene", "distunder0", "area")) %>% distinct()) %>% filter(gene != "ACTB" & gene != "CMV" & gene != "EEF1A1") 


lcl_50_dist_5hun_g %<>% pivot_wider(id_cols = 'gene', names_from = 'distunder0', values_from = 'area')

lcl_50_dist_5hun_g$fraction <- lcl_50_dist_5hun_g$yes / (lcl_50_dist_5hun_g$no + lcl_50_dist_5hun_g$yes)*100

lcl_50_dist_5hun_g$celltype <- 'LCL'

lcl_50_dist_5hun_g$up_frac <- lcl_50_dist_5hun_g$yes / (lcl_50_dist_5hun_g$no + lcl_50_dist_5hun_g$yes)*100
lcl_50_dist_5hun_g[13,6] <- 100 #make 100 bc was NA 

lcl_50_dist_5hun_g$down_frac <- lcl_50_dist_5hun_g$no / (lcl_50_dist_5hun_g$no + lcl_50_dist_5hun_g$yes)*100


#fraction of total activity 500bp upstream of TSS 
bound_tog <- rbind(full_fib_50_dist_5hun_g,lcl_50_dist_5hun_g)

bound_tog <- bound_tog %>% group_by(celltype) %>%
  mutate(gene = reorder(gene, up_frac))

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/frac_AUC_TSS.pdf"), width=6, height=4, colormodel = "rgb")

ggplot(bound_tog, aes(x=gene, y = up_frac, fill = celltype)) +
  geom_col(position = "dodge") + 
  theme_pubr()
dev.off()









