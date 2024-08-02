# Weakening of Y promoters shaped the evolution of dosage-sensitive gene expression on the human sex chromosomes

## This repository contains code necessary to recreate analyses in: 

Harris HL, Skaletsky H, Keys H, Viswanathan K, Page DC. 2024. Weakening of Y promoters shaped the evolution of dosage-sensitive gene expression on the human sex chromosomes. 
___ 

Each folder contains code and processed data that may be used to recreate one or more figures and tables:  
 
software versions used: 
analyses were conducted using either Python (v3.8.10) or R (v4.1.0 or v4.2.0). Software packages used in Python included Bio (v1.76), pandas (v1.4.2), pyranges (v0.0.79), genomicranges (v0.2.11), and gtfparse. Software packages used in R included GenomicRanges (v1.46.1), seqinr (v4.2-8), dplyr (v1.0.7), ggpubr (v0.4.0), stringr (v1.4.0), ape (v5.6), BSgenome (v1.62.0), tidyr (v1.1.4), ggplot2 (v3.3.5), ggrepel (v0.9.1), rstatix (v0.7.0), gt (v0.10.1), edgeR (v3.36.0), and purrr (v0.3.4).

## Generate annotation files by running:
  1. alignments/get_gtf_git.ipynb
  2. alignments/gene_by_gene_gtf.R

## Figure 1,2 + supplement: X vs Y homolog sequence comparisons, Y promoters have lost G+C content 

 generate alignment, calculate % alignment, calculate GC % content:  
 * alignments/pair_align_multiz_final_git.ipynb

 dog/bull % alignment, GC content: 
 * alignments/align_humandogbull_git.ipynb 

 opossum/chicken GC content: 
 * GC_analysis/chicken_GC_git.ipynb
 * GC_analysis/op_GC_git.ipynb

 make figures:
 * GC_analysis/f1_graphs.R

## Figure 3 + supplement: Multi-species alignments reveal faster rate of Y evolution 

 generate primate alignment: 
 * alignments/multiz_7 git.ipynb

 polish primate alignments:
 * alignments/multiz_7sp_final_git.ipynb

 get GC normalized intron alignments: 
 * GC_analysis/gc_content_normalization_git.ipynb

 make figures:
 * pair_evolution/draw_trees.R 

## Figure 4: Y homologs are under less constraint than X homologs 

  make figure 4A + constraint analysis: 
  * pair_evolution/draw_trees.R 

  make figure 4B, negative selection testing: 
  * human_variation
  * human_variation/pull_vcf_info.py to make VCF files per gene/region 
  * human_variation/make_bedfiles.R to make bedfiles of gene regions
  * human_variation/get_regions.py and human_variation/get_regions.sh to make csv of all variants per gene/region
  * human_variation/re_analysis_0206.R to generate file of variant counts per gene/region
  * human_variation/divergence_v1.R to compare to between species variation and generate figure 


## Figure 5: Quantitative assessment of promoter activity

  make pool of elements
  *mpra_design/makeElementPool_HH_new.py

  associate barcodes with candidate regulatory sequences (CRSs)
  *mpra_design/new_mpra_pipeline.sh

  run count pipeline: 
  *mpra_design/count_pipeline.sh

## Figure 6: mpra analysis

  get distance of each start site to CRS start: 
  *mpra_analysis/distance_calcs_TSS_git.ipynb

  mpra functions: 
  *mpra_analysis/mpra_functions.R

  make mpra figures: 
  *mpra_analysis/mpra_final.R

## Figure 7: 
 get number of genes: 
 *pair_gene_expression/fig7.R

## Supplement TF analysis: 

 figures and analysis: 
 *TF_analysis/TF_analysis.R


## Supplement gene expression: 

 *pair_gene_expression/plot_tpms_v1_git.Rmd 
 
## ZFX promoter evolution: 

 human-opossum alignment: 
 *ZFX_promoter_evolution/align_humanopossum_zfx.py
 
 make figures:
 *ZFX_promoter_evolution/human_opossum_graphs_git.Rmd
















