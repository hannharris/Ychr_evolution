---
title: "mpra functions"
output: html_notebook
---
#write_50_df
```{r}
write_50_df <- function(filename, metadata, dataframe){
merged <- merge(metadata, dataframe, by = "name")
merged <- merged %>% filter(type != "negative ctl")
merged <- merged %>% filter(n_obs_bc > 1) #n_obs > 2 
fib_1_50 <<- calculate_50bp_windows(merged, metadata, "counts") #this also saves a BED file
#colnames(fib_1_50)[colnames(fib_1_50) == "OG_coords"] <- "name"
m_join <<- left_join(fib_1_50, metadata, by = "name" )
m_join <- m_join[,c("coords", "log2", "sd", "lower_val", "upper_val", "gene", "orientation", "type")]
#colnames(m_join) <- c("coords", "log2", "gene", "orientation", "type")
write.table(x = m_join, file = filename, quote = FALSE, sep = "\t",  row.names = FALSE)

return(fib_1_50)
}
```
#write_200_df
```{r}
write_200_df <- function(filename, metadata, dataframe){
merged <- merge(metadata, dataframe, by = "name")
merged <<- merged %>% filter(n_obs_bc > 1) #n_obs > 2 


m_join <- merged[,c("name", "gene", "orientation", "type", "log2")]
m_join$coords <- m_join$name
m_join$name <- NULL
write.table(x = m_join, file = filename, quote = FALSE, sep = "\t",  row.names = FALSE)
return(m_join)
}

```

#make bed
```{r}
make_bedgraph <- function(input_file, file_name, column_to_sep){
  #input needs to have distances (so you can get rid of the shuff seqeucnes)
  #function removes CMV 
  #saves the bed file to a specified name and returns the bed 
  
  make_correct_bed_to_edit <- function(input_file) {
     tryCatch(expr =  {return(input_file %>% dplyr::filter(!is.na(distance)) %>% dplyr::filter(gene.y != "CMV"))
    }, 
    warning = function(w) {
      return(input_file)
    }, error = function(e){
      return(input_file)
  })}
 
  bed <<- make_correct_bed_to_edit(input_file)
 
  bed <- separate(bed, column_to_sep, c("chromosome", 'coords1'), sep = (":"), remove = FALSE)
  bed <- separate(bed, coords1, c("start", 'end'), sep = ("-"), remove = FALSE)

  new_bed <- bed %>% dplyr::select("chromosome", "start", "end", "log2")

  write.table(x = new_bed, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE)
  return(list(bed, new_bed))
}

```
#calculate_50bp_windows
```{r}
calculate_50bp_windows <- function(input_file, metadata, cell_type){ #metadata,
  #takes in a file like lcl_distances, metadata (new_meta_), and the cell type 
  #creates a bed file
  #saves two files - a bedgraph of the log2 for every 50bp and also a file with genes, coords, annotation of direction of gene for future distance calculations 

  bed <- make_bedgraph(input_file, "bed", "name")[[1]]

  list_chromosome <- c()
  list_starts <- c()
  list_ends <- c()
  list_b <- c()
  OG_coords <- c()
  bed$start <- as.numeric(bed$start)
  bed$end <- as.numeric(bed$end)
 
for (row in 1:nrow(bed)){
    #print('here')
    new_df <- bed[row,]
    name <-  as.character(new_df$name)
    coord2 <- new_df$start + 50  
    coord3 <- coord2 + 50
    coord4 <- coord3 + 50
    list_starts <- c(list_starts, new_df$start, coord2, coord3, coord4)
    list_ends <- c(list_ends, coord2, coord3, coord4, new_df$end + 1)
    list_b <- c(list_b, new_df$log2, new_df$log2, new_df$log2, new_df$log2)
    list_chromosome <- c(list_chromosome, new_df$chromosome,new_df$chromosome, new_df$chromosome, new_df$chromosome)
    OG_coords <- c(OG_coords, name, name, name, name)
}

fifty_bases_df <<- data.frame('start' = list_starts, 'end' = list_ends, 'chromosome' = list_chromosome, 'OG_coords' = OG_coords, 'log2' = list_b)
fifty_bases_df$coords <<- paste(fifty_bases_df$chromosome, ":", fifty_bases_df$start, "-", fifty_bases_df$end, sep = "")

fifty_bases_df_to_avg <<- fifty_bases_df %>% dplyr::select("log2", "coords")

coord_options <- unique(fifty_bases_df_to_avg$coords)
coordins <- c()
sds <- c()
for (coord in coord_options){
  uniq <-fifty_bases_df_to_avg %>% filter(coords == coord)
  stddev <- sd(uniq$log2)
  sds <- c(sds, stddev)
  coordins <- c(coordins, coord)
}
sd_df <- data.frame('coords' = coordins, 'sd' = sds)


get_averages <- fifty_bases_df_to_avg %>% pivot_wider(names_from = coords, values_from = log2, values_fn = function(x){mean(x,na.rm = T)}) #check to make sure this actually works correctly 

fifty_bp_avg <<- get_averages %>% pivot_longer(cols = 1:length(get_averages), names_to = "coords", values_to = "log2")

filename_bed = paste("50bp_bed_test_", cell_type, ".txt", sep = "")

new_bed <- make_bedgraph(fifty_bp_avg, filename_bed, "coords")

uni <- fifty_bases_df %>% dplyr::select("OG_coords", 'coords')
uni <- dplyr::distinct(uni, coords, .keep_all = TRUE)

fifty_bp_averages_w_genenames <- merge(fifty_bp_avg, sd_df, by = "coords")
fifty_bp_averages_w_genenames[is.na(fifty_bp_averages_w_genenames)] = 0 #set some NA sd to 0's
fifty_bp_averages_w_genenames$lower_val <- fifty_bp_averages_w_genenames$log2 - fifty_bp_averages_w_genenames$sd
fifty_bp_averages_w_genenames$upper_val <- fifty_bp_averages_w_genenames$log2 +  fifty_bp_averages_w_genenames$sd

fifty_bp_averages_w_genenames <- dplyr::left_join(fifty_bp_averages_w_genenames, uni, by = "coords")
fifty_bp_averages_w_genenames$coordinate <- fifty_bp_averages_w_genenames$OG_coords
fifty_bp_averages_w_genenames$name <- fifty_bp_averages_w_genenames$coordinate
fifty_bp_averages_w_genenames <<- merge(fifty_bp_averages_w_genenames, metadata, by = 'name')

return(fifty_bp_averages_w_genenames)
}

```

#distance from TSS 
```{r}

get_dist_TSS <- function(lcl_50, start_sites){
  
  lcl_1_dist <- merge(lcl_50, start_sites,by = "gene") 
  lcl_1_dist <- separate(lcl_1_dist, 'coords', c("chromosome", 'coords1'), sep = (":"), remove = FALSE)
  lcl_1_dist <- separate(lcl_1_dist, coords1, c("mpra_start", 'mpra_end'), sep = ("-"), remove = FALSE)
  lcl_1_dist$mpra_start <- as.numeric(lcl_1_dist$mpra_start)
  lcl_1_dist$mpra_end <- as.numeric(lcl_1_dist$mpra_end)
  lcl_1_dist$mpra_end <- as.numeric(lcl_1_dist$mpra_start)
  lcl_1_dist$start <- as.numeric(lcl_1_dist$start)

  lcl_1_dist$distance <- ifelse(lcl_1_dist$orientation == 'plus',lcl_1_dist$mpra_start - lcl_1_dist$start ,  lcl_1_dist$start - lcl_1_dist$mpra_end)
  return(lcl_1_dist)

}
```

#graph promoters
```{r}

per_pair_dist_graphs <- function(distance_df, distance_cutoff, cell_type, gene_pairs){

  promoters_dist1 <<- distance_df %>% filter(distance < distance_cutoff & distance > (-1*distance_cutoff))
  meg_dist_g_p <<- merge(promoters_dist1, gene_pairs, by = "gene") 
  
  filename = paste("/lab/solexa_page/hannah/supp_info/figures/", cell_type, ".pdf", sep = "")
  pdf(file=(filename), width=6, height=6, colormodel = "rgb")
  print(ggplot(meg_dist_g_p, aes(x=distance, y=log2, color = type_of_gene)) +  #, shape =TSS_num
    geom_line() + 
    geom_point(stroke = 0) + 
    facet_wrap(vars(pair)) + 
    ggpubr::theme_pubr() + 
   #geom_pointrange(aes(ymin = lower_val, ymax = upper_val)) + 
    scale_colour_manual(values = c("#e46915", "#857bc6")) +
    ylab("log2(RNA/DNA)")) 
  dev.off()

}

```

#get AUC
```{r}
get_AUC <- function(promoters_dist1, distance_cutoff, gene_pairs) {
  #finds AUC values for each X-Y gene pair, returns a df of those values
  lcl_distance_cropped <<- promoters_dist1 %>% filter(distance < distance_cutoff & distance > (-1*distance_cutoff))

  lcl_distance_cropped$log2_int <- ifelse(lcl_distance_cropped$log2 > 0, lcl_distance_cropped$log2,  0) 
  areas <- c()
  genes <- c()
  for (gene_1 in unique(lcl_distance_cropped$gene)){
    
    crop_lcl_distance_cropped <- lcl_distance_cropped %>% filter(gene == gene_1)
    if (gene_1 == "ZFY" | gene_1 == "ZFX"){
      crop_lcl_distance_cropped <- crop_lcl_distance_cropped %>% filter(distance >= -280 & distance <= 50)
    }
    
    area <- AUC(crop_lcl_distance_cropped$distance, crop_lcl_distance_cropped$log2_int, method="trapezoid")
    areas <- c(areas, area)
    genes <- c(genes, gene_1)
  }
  aucs <- data.frame(areas, gene = genes)
  print(aucs)
  aucsm <- merge(aucs, gene_pairs, by = "gene")
  
  return(aucsm)
}
```

#plot AUC
```{r}

plot_AUC <- function(f_auc, filename, xmax, ymax){
  
  f_auc_new <- f_auc %>% select("type_of_gene", "pair", "areas")
  f_auc_new <- pivot_wider(f_auc_new, id_cols = 'pair', names_from = 'type_of_gene', values_from = 'areas')

  f_auc_new$X_gene <- as.numeric(f_auc_new$X_gene)
  f_auc_new$Y_gene <- as.numeric(f_auc_new$Y_gene)
  
  stat.test <- wilcox.test(f_auc_new$X_gene, f_auc_new$Y_gene, paired = TRUE)
  stat.test$p.value
pdf(file=(filename), width=2.3, height=2.3, colormodel = "rgb")
 
  ggplot(f_auc_new, aes(x = X_gene, y = Y_gene)) + 
    geom_point(stroke = 0) + 
    theme_pubr() + 
    geom_text_repel(aes(label = pair)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0))) + 
    
    annotate("text", x = 100, y = 500, 
           label = paste("p =", format(stat.test$p.value, digits = 2)),
           size = 5, color = "black")
#dev.off()
}


```
#get maxpeak
```{r}

max_peak <- function(lcl_sub){
 # lcl_sub <- full_fib_50_dist1
  lcl_sub$gene.x <- lcl_sub$gene
  lcl_sub$gene <- NULL
  lcl_subsetted_peaks_top <- data.frame()
  genes <- unique(lcl_sub$gene.x)
  for (gene in genes) { #finding top 5 around the max peak 
      print(gene)
      lcl_sub_gene <- lcl_sub %>% filter(gene.x == gene)
      single_line <- lcl_sub_gene %>% filter(log2 == max(lcl_sub_gene$log2))
      lcl_subsetted_peaks_top <- rbind(lcl_subsetted_peaks_top, single_line) 
  }
  lcl_subsetted_peaks_top$gene <- lcl_subsetted_peaks_top$gene.x 
  lcl_subsetted_peaks_top$gene.x <- NULL
  return(lcl_subsetted_peaks_top)
}

```
#plot maxpeak
```{r}

plot_maxpeak <- function(f_auc, filename, xmax, ymax){
 # f_auc <- max_full_fib
  f_auc_new <- f_auc %>% select("type_of_gene", "pair", "log2")
  f_auc_new <- pivot_wider(f_auc_new, id_cols = 'pair', names_from = 'type_of_gene', values_from = 'log2')

  f_auc_new$X_gene <- as.numeric(f_auc_new$X_gene)
  f_auc_new$Y_gene <- as.numeric(f_auc_new$Y_gene)
  
  stat.test <- wilcox.test(f_auc_new$X_gene, f_auc_new$Y_gene, paired = TRUE)
  print(stat.test$p.value)
  
  
pdf(file=(filename), width=2.3, height=2.3, colormodel = "rgb")
  ggplot(f_auc_new, aes(x = X_gene, y = Y_gene)) + 
    geom_point(stroke = 0) + 
    theme_pubr() + 
    geom_text_repel(aes(label = pair)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    #xlim(0,xmax) + 
    #ylim(0,ymax) + 
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0))) + 
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0))) + 
    annotate("text", x = 0.25, y = 2, 
           label = paste("p =", format(stat.test$p.value, digits = 2)),
           size = 5, color = "black")
#dev.off()
}

```







