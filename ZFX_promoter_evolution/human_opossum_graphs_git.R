

myPath <- #PATH TO GITHUB FOLDER


h_o <- read.delim(paste0(myPath, "/tables/mutation_types_ZFX_promoter.txt"), sep = ",")
colnames(h_o) <- c("mut_types", "X", "Y")

# new_row <- c("total", sum(h_o$X), sum(h_o$Y ))
# 
# h_o$X_total <- 271 #sum(h_o$X)
# h_o$Y_total <- 271 #sum(h_o$Y)

h_o$A_total <- 67 #sum(h_o$X)
h_o$C_total <- 204 #sum(h_o$Y)


h_ocg <- h_o[1:3,]
h_ocg$perc_X <- (h_ocg$X / h_ocg$C_total)*100
h_ocg$perc_Y <- (h_ocg$Y / h_ocg$C_total)*100


h_oat <- h_o[4:6,]
h_oat$perc_X <- (h_oat$X / h_oat$A_total)*100
h_oat$perc_Y <- (h_oat$Y / h_oat$A_total)*100

h_o <- rbind(h_ocg, h_oat)

h_o <- h_o[,c(1,6,7)]


#graph


h_o1 <- as.data.frame(t(as.matrix(h_o)))

h_o2 <- h_o %>% pivot_longer(cols = c(2,3))


#pdf(file=paste0(myPath, "/figures/1_opossum_ZFX_muttypes.pdf"), width=3, height=3, colormodel = "rgb")

ggplot(h_o2, aes(x=mut_types, y = value, fill = name)) + 
  geom_col(position = "dodge", width = 0.5) +
  theme_pubr() + 
  ylab("substitution type percent")+
    theme(legend.position = "none") 

dev.off()




#introns


h_o <- read.delim(paste0(myPath, "/tables/mutation_types_ZFX_intron.txt"), sep = ",")
colnames(h_o) <- c("mut_types", "X", "Y")



h_o$A_total <- 210 #sum(h_o$X)
h_o$C_total <- 95 #sum(h_o$Y)
#
h_ocg <- h_o[1:3,]
h_ocg$perc_X <- (h_ocg$X / h_ocg$C_total)*100
h_ocg$perc_Y <- (h_ocg$Y / h_ocg$C_total)*100
#
#
h_oat <- h_o[4:6,]
h_oat$perc_X <- (h_oat$X / h_oat$A_total)*100
h_oat$perc_Y <- (h_oat$Y / h_oat$A_total)*100
#
h_o <- rbind(h_ocg, h_oat)

h_o <- h_o[,c(1,6,7)]

h_o1 <- as.data.frame(t(as.matrix(h_o)))

h_o2 <- h_o %>% pivot_longer(cols = c(2,3))

ggplot(h_o2, aes(x=mut_types, y = value, fill = name)) + 
  geom_col(position = "dodge", width = 0.5) +
  theme_pubr() + 
  ylab("substitution type percent") +
  scale_fill_manual(values = c("#DE8034", "#7772AF"))








