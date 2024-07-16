
#run with R 4.2.1

#load libraries

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(stringr)
library(rstatix)
library(gt)
library(gtExtras)

#load files

g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") 
gc <- read.delim("/lab/solexa_page/hannah/supp_info/tables/GC_1030.csv", sep = ",")[2:5]
gc <- gc %>% filter(pairs != "NLGN4X_NLGN4Y" & pairs != "TXLNG_TXLNGY" & pairs != "TMSB4X_TMSB4Y")
gc_metrics <- gc %>% group_by(pairs, region) %>% mutate(difference_in_gc = GC_perc[which(gene == "X_gene")] - GC_perc[which(gene == "Y_gene")]) 
gc_metrics <- gc_metrics %>% group_by(region) %>% mutate(median_gc = mean(difference_in_gc))
perc <- read.delim("/lab/solexa_page/hannah/supp_info/tables/percent_alignment_1030.csv", sep = ",")[2:4]
perc <- perc %>% filter(pair != "NLGN4X_NLGN4Y" & pair != "TXLNG_TXLNGY" & pair != "TMSB4X_TMSB4Y")

perc_scramble <- read.delim("/lab/solexa_page/hannah/supp_info/tables/percent_alignment_scramble.csv",sep = ",")[2:4]

#GC graphs

stat.test <- gc %>% filter(GC_perc > 0 ) %>% group_by(region) %>%
  wilcox_test(GC_perc ~ gene, paired = TRUE) %>%
  add_significance()


gc1 <- gc %>% filter(GC_perc > 0)  %>% pivot_wider(names_from = "gene", values_from = "GC_perc") 
    
    
pdf(file=("/lab/solexa_page/hannah/supp_info/figures/1D_GC_scatter.pdf"), width=3, height=1.5, colormodel = "rgb")
gc1 %>% ggplot(aes(x=X_gene, y=Y_gene)) + 
  geom_point(stroke = 0) + 
  facet_grid(~region) + 
  theme_pubr() + 
  #geom_text_repel(aes(label = pairs)) + 
  xlim(30,80) + 
  ylim(30,80) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
dev.off()

 

#graph percents

perc$is_pair <- ifelse(perc$pair %in% g_p$pair, paste(perc$region, "yes"), paste(perc$region, "no"))
perc_scramble$is_pair<- "promoter scram"
perc <- rbind(perc, perc_scramble)
perc$is_pair <- factor(perc$is_pair, levels = c("exon yes", "intron yes", "promoter yes", "promoter scram", "promoter no"))

#sample 7 random values for random combos of promoters#
random_values <- perc %>%
  filter(is_pair == "promoter no") %>%
  sample_n(7, replace = FALSE)

perc <- perc %>% filter(percent > 0 ) %>% filter(is_pair != 'promoter no')
perc <- rbind(random_values, perc)


#Outliers assumption
perc %>%
  group_by(is_pair) %>%
  identify_outliers(percent)
# Build the linear model
model  <- lm(percent ~ is_pair,
             data = perc)
# # Create a QQ plot of residuals
ggqqplot(residuals(model))
shapiro_test(residuals(model))
perc %>%
  group_by(is_pair) %>%
  shapiro_test(percent)
ggqqplot(perc, "percent", ggtheme = theme_bw()) +
facet_grid(~is_pair)
# #Homogeneity of variance assumption
perc %>% levene_test(percent ~ is_pair)

res.kruskal <- perc %>% filter(is_pair != "promoter no") %>% kruskal_test(percent ~ is_pair)

library(emmeans)
perc %>%
  anova_test(percent ~ is_pair, error = model)
pwc <- perc %>% 
  emmeans_test(percent ~ is_pair, p.adjust.method = "holm") 
pwc <- pwc %>%  add_xy_position(x = "is_pair",dodge = 1,step.increase = .1)


violin_plot <- ggviolin(perc %>% filter(percent > 0), x = "is_pair", y = "percent",trim = TRUE ) + #add = "point"
   geom_point(
    stroke = 0,  # Set the stroke width
  ) 

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/1C_percentalign.pdf"), width=3, height=4, colormodel = "rgb")
violin_plot + 
  #stat_pvalue_manual(pwc[c(1,2,4,6),],   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 8))  
  stat_pvalue_manual(pwc[c(1,2,5,8,9,10),],   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 0))
dev.off()



#make table 1b

tableb <- read.delim('/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/pair_align/fig1b.txt', check.names = FALSE)
tableb$`approximate divergence time (MYA)` <- NULL
tableb$CDS <- tableb$`% identity CDS` 
tableb$promoter <- tableb$`% identity promoter`
tableb$`% identity CDS` <- NULL
tableb$`% identity promoter` <- NULL
  
table_gt <- tableb %>% gt() %>%
   tab_options(table.font.names = c("Helvetica Neue")) %>% 
  #  tab_header(
  #   title = "X-Y homologs") %>% 
  tab_spanner(
    label = "X-Y percent identity",
    columns = c(CDS, promoter)
  ) %>% 
  cols_align(
  align = c( "center"),
  columns = c(1, 2,3,4) 
) 

table_gt
gt::gtsave(table_gt, "my_table.docx") #can just drag into illustrator



#graph_rolling


perc_rolling <- read.delim("/lab/solexa_page/hannah/supp_info/tables/perc_rolling_1208.csv", sep = ",")

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S3_supprolling.pdf"), width=8, height=6, colormodel = "rgb")

perc_rolling %>% pivot_longer(perc_rolling, cols = c(5,6), names_to = "place", values_to = "coord") %>% 
  ggplot(aes(coord, percent)) + 
  facet_wrap(~pair) + 
  geom_line() + 
  theme_pubr() 

dev.off()



#CpG content


cpg <- read.delim("/lab/solexa_page/hannah/supp_info/tables/CpG_1127.csv", sep = ",")

cpg1 <- cpg %>% select(!"raw_cpGs" & !"X") %>% pivot_wider(names_from = "gene", values_from = "CpG_norm") 
wilcox.test(cpg1$X_gene, cpg1$Y_gene, paired = TRUE)

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S6_CpG.pdf"), width=4, height=4, colormodel = "rgb")
cpg1 %>% ggplot(aes(x=X_gene, y=Y_gene)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  geom_text_repel(aes(label = pairs)) +
  #xlim(0,15) +
  #ylim(0,15) +
  scale_x_continuous(limits = c(0.5, 1.2), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0.5, 1.2), expand = expansion(mult = c(0, 0))) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

dev.off()

cpg1 <- cpg %>% select(!"CpG_norm" & !"X") %>% pivot_wider(names_from = "gene", values_from = "raw_cpGs") 

wilcox.test(cpg1$X_gene, cpg1$Y_gene, paired = TRUE)

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/figS6_scatter.pdf"), width=4, height=4, colormodel = "rgb")

cpg1 %>% ggplot(aes(x=X_gene, y=Y_gene)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  geom_text_repel(aes(label = pairs)) +
  #xlim(0,15) +
  #ylim(0,15) +
  scale_x_continuous(limits = c(0, 15), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 15), expand = expansion(mult = c(0, 0))) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
dev.off()



#human-dog-bull

perc_1 <- perc %>% filter(is_pair == "intron yes" | is_pair == "exon yes" | is_pair == "promoter yes")
perc_1$is_pair <- NULL
perc_1$species <- 'humanX:humanY'
hdb <- read.delim('/lab/solexa_page/hannah/supp_info/tables/percent_alignment_1211_dogbull.csv', sep = ",")[2:5]
hdb <- rbind(perc_1, hdb)

#Outliers assumption
hdb %>%
  group_by(region, species) %>%
  identify_outliers(percent)
# Build the linear model
model  <- lm(percent ~ region*species,
             data = hdb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
shapiro_test(residuals(model))
hdb %>%
  group_by(region, species) %>%
  shapiro_test(percent)
ggqqplot(hdb, "percent", ggtheme = theme_bw()) +
  facet_grid(region ~ species)
#Homogeneity of variance assumption
hdb %>% levene_test(percent ~ region*species)

res.kruskal <- hdb %>% group_by(region) %>% kruskal_test(percent ~ species)
hdb %>% group_by(region) %>% kruskal_effsize(percent ~ species)
pwc <- hdb %>% filter(percent > 0) %>% group_by(region) %>% 
  dunn_test(percent ~ species, p.adjust.method = "holm") 
pwc$p.adj.sci <- format(pwc$p.adj, scientific = FALSE, digits = 3)

pwc <- pwc %>%  add_xy_position(x = "region",dodge = 1,step.increase = .1)

hdb$region = factor(hdb$region, levels=c('exon', 'intron', 'promoter')) 

violin_plot <- ggviolin(hdb %>% filter(percent > 0), x = "region", y = "percent", fill = "species", trim = TRUE)  +
  geom_point(aes(x=region, y=percent, fill = species), position = position_dodge(width = 0.8),
    stroke = 0,  
  ) 

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S4_percentIDdogbull.pdf"), width=4, height=4, colormodel = "rgb")
violin_plot + 
  stat_pvalue_manual(pwc %>% filter(region != "exon"),   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 17)) + 
  stat_pvalue_manual(pwc %>% filter(region == "exon"),   label = "p.adj.signif", tip.length = 0.01, position = position_nudge(y = 5))
dev.off()



#dog-bull GC percents


gc_p <- read.delim('/lab/solexa_page/hannah/supp_info/tables/percent_GC_1211_dogbull.csv', sep = ",")[,2:5]

g_p$gene <- g_p$gene.x

gc_ms <- unique(merge(gc_p, g_p, by = 'gene') %>% select("type_of_gene", "species", "gc", "region", "gene"))
gc_ms1 <- gc_ms %>% pivot_wider(names_from = "species", values_from = "gc") 

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S5_dog_GC.pdf"), width=5, height=2.25, colormodel = "rgb")
gc_ms1 %>% ggplot(aes(x=human, y=dog)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  #geom_text_repel(aes(label = gene)) +
  facet_wrap(~region) + 
  xlim(30,80) +
  ylim(30,80) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
dev.off()

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/S5_bull_GC.pdf"), width=5, height=2.25, colormodel = "rgb")

gc_ms1 %>% ggplot(aes(x=human, y=bull)) +
  geom_point() +
  theme_pubr() +
  geom_text_repel(aes(label = gene)) +
  facet_wrap(~region) + 
  xlim(30,80) +
  ylim(30,80) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

dev.off()


#human/chicken promoter GC

g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") 

gc <- read.delim("/lab/solexa_page/hannah/supp_info/tables/GC_1030.csv", sep = ",")[2:5]
gc <- gc %>% filter(pairs != "NLGN4X_NLGN4Y" & pairs != "TXLNG_TXLNGY" & pairs != "TMSB4X_TMSB4Y")

gc_chicken <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_dog_bull/chicken/chicken_gc.txt", sep = ",")
g_p <- g_p %>% rename('gene' = 'gene.x')

gc_chicken <- merge(g_p, gc_chicken, by = "gene")

gc_X <- gc %>% filter(gene == "X_gene" & region == 'promoter') %>% rename('pair' = 'pairs')
gc_Y <- gc %>% filter(gene == "Y_gene" & region == 'promoter') %>% rename('pair' = 'pairs')

x_genes_chicken_human <- merge(gc_chicken, gc_X, by = "pair")
x_genes_chicken_human$xg <- "x to x"
y_genes_chicken_human <- merge(gc_chicken, gc_Y, by = "pair")
y_genes_chicken_human$xg <- "y to y"
rbound <- rbind(x_genes_chicken_human, y_genes_chicken_human)

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/1_chicken_GC.pdf"), width=2.25, height=2.25, colormodel = "rgb")

rbound %>% ggplot(aes(x=gc, y=GC_perc, color = xg)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  #geom_text_repel(aes(label = gene)) +
 # facet_wrap(~region) + 
  xlim(40,83) +
  ylim(40,83) +
  geom_text_repel(aes(label = gene.x)) +
    theme(legend.position = "none") + 

  ylab("human sex chr promoters") + 
  xlab("chicken promoter G+C%") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

dev.off()

rbound_X <- rbound %>% filter(xg == "x to x")
rbound_Y <- rbound %>% filter(xg == "y to y")

wilcox_test_resultX <- wilcox.test(rbound_X$gc, rbound_X$GC_perc, paired = TRUE)
wilcox_test_resultY <- wilcox.test(rbound_Y$gc, rbound_Y$GC_perc, paired = TRUE)



#human/opossum promoter GC

g_p <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/g_p.txt") 

gc <- read.delim("/lab/solexa_page/hannah/supp_info/tables/GC_1030.csv", sep = ",")[2:5]
gc <- gc %>% filter(pairs != "NLGN4X_NLGN4Y" & pairs != "TXLNG_TXLNGY" & pairs != "TMSB4X_TMSB4Y")

gc_chicken <- read.delim("/lab/solexa_page/hannah/220516_mpra/msa/long_alignments/multiz_dog_bull/opossum/opossum_gc.txt", sep = ",")
g_p <- g_p %>% rename('gene' = 'gene.x')

gc_chicken <- merge(g_p, gc_chicken, by = "gene")

gc_X <- gc %>% filter(gene == "X_gene" & region == 'promoter') %>% rename('pair' = 'pairs')
gc_Y <- gc %>% filter(gene == "Y_gene" & region == 'promoter') %>% rename('pair' = 'pairs')

x_genes_chicken_human <- merge(gc_chicken, gc_X, by = "pair")
x_genes_chicken_human$xg <- "x to x"
y_genes_chicken_human <- merge(gc_chicken, gc_Y, by = "pair")
y_genes_chicken_human$xg <- "y to y"
rbound <- rbind(x_genes_chicken_human, y_genes_chicken_human)

pdf(file=("/lab/solexa_page/hannah/supp_info/figures/1_opossum_GC.pdf"), width=2.25, height=2.25, colormodel = "rgb")
rbound %>% ggplot(aes(x=gc, y=GC_perc, color = xg)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  #geom_text_repel(aes(label = gene)) +
  #facet_wrap(~region) + 
  xlim(40,83) +
  ylim(40,83) +
  geom_text_repel(aes(label = gene.x)) + 
  theme(legend.position = "none") + 
  ylab("human sex chr promoters") + 
  xlab("opossum promoter G+C%") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
dev.off()



pdf(file=("/lab/solexa_page/hannah/supp_info/figures/2_opossum_GC.pdf"), width=2.25, height=2.25, colormodel = "rgb")
rbound1 <- rbound %>% filter(pair != "RPS4X_RPS4Y1" & pair != "KDM5C_KDM5D")
rbound1 %>% ggplot(aes(x=gc, y=GC_perc, color = xg)) +
  geom_point(stroke = 0) +
  theme_pubr() +
  #geom_text_repel(aes(label = gene)) +
  #facet_wrap(~region) + 
  xlim(40,83) +
  ylim(40,83) +
  geom_text_repel(aes(label = gene.x)) + 
  theme(legend.position = "none") + 
  ylab("human sex chr promoters") + 
  xlab("opossum promoter G+C%") + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
dev.off()

# wilcox_test_result <- wilcox.test(gc1_promoter$X_gene, gc1_promoter$Y_gene, paired = TRUE)

rbound_X <- rbound %>% filter(xg == "x to x")
rbound_Y <- rbound %>% filter(xg == "y to y")

wilcox_test_resultX <- wilcox.test(rbound_X$gc, rbound_X$GC_perc, paired = TRUE)
wilcox_test_resultY <- wilcox.test(rbound_Y$gc, rbound_Y$GC_perc, paired = TRUE)




