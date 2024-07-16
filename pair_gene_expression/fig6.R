---
title: "R Notebook"
output: html_notebook
---
#load files


s2_asr <- read.delim('S2_ASR.txt')
s6_asr <- read.delim('S6_ASR.txt')
s3_naqvi <- read.delim('Table_S3_naqvi.txt')

s26 <- full_join(s2_asr, s6_asr, by = "Gene")


#split into lcl 


s26_lcl <- s26 %>% filter(!is.na(deltaEx_LCL))
nrow(s26_lcl %>% filter(Revisedstatus == 'Escape'))
nrow(s26_lcl %>% filter(Revisedstatus == 'Subject'))
nrow(s26_lcl %>% filter(Revisedstatus == 'No call'))

s26_lcl_XY <- nrow(s26_lcl %>% filter(Revisedstatus == 'Escape' & NPX.NPY_gene == TRUE))

s26_lcl_escape <- nrow(s26_lcl %>% filter(Revisedstatus == 'Escape' & NPX.NPY_gene == FALSE))

s26_lcl_subject <- nrow(s26_lcl %>% filter(Revisedstatus == 'Subject' & NPX.NPY_gene == FALSE))

s26_lcl_subject_XY <- nrow(s26_lcl %>% filter(Revisedstatus == 'Subject' & NPX.NPY_gene == TRUE))

print(paste(s26_lcl_XY, "_", s26_lcl_escape, "_", s26_lcl_subject, "_", s26_lcl_subject_XY))

#split into fib 

s26_lcl <- s26 %>% filter(!is.na(deltaEx_Fib))
nrow(s26_lcl %>% filter(Revisedstatus == 'Escape'))
nrow(s26_lcl %>% filter(Revisedstatus == 'Subject'))
nrow(s26_lcl %>% filter(Revisedstatus == 'No call'))
nrow(s26_lcl %>% filter(is.na(Revisedstatus)))

s26_lcl_XY <- nrow(s26_lcl %>% filter(Revisedstatus == 'Escape' & NPX.NPY_gene == TRUE))

s26_lcl_escape <- nrow(s26_lcl %>% filter(Revisedstatus == 'Escape' & NPX.NPY_gene == FALSE))

s26_lcl_subject <- nrow(s26_lcl %>% filter(Revisedstatus == 'Subject' & NPX.NPY_gene == FALSE))

s26_lcl_subject_XY <- nrow(s26_lcl %>% filter(Revisedstatus == 'Subject' & NPX.NPY_gene == TRUE))

print(paste(s26_lcl_XY, "_", s26_lcl_escape, "_", s26_lcl_subject, "_", s26_lcl_subject_XY))




