
library(stringr)
a=read.csv('../phenotypes/Ames/All_Crosses_list.csv')
b=read.csv('../phenotypes/Ames/BLUEs.csv')
fq=read.table('../../282_load/282_fqnames.txt', header=F)
fq$name=str_split_fixed(fq$V1,'_',3)[,1]
a$name=str_split_fixed(a$Hybrid_pheno_name, paste0('x', a$Tester), 2)[,1]

sum(a$name %in% fq$name) ## still only 235 not 282, but a start!
sum(toupper(unique(a$name)) %in% toupper(fq$name)) ## up to 244

## 283 individuals
length(unique(fq$name))
sum(toupper(unique(fq$name)) %in% toupper(a$name)) ## 138 oh duh because there are two testers...

## what am i missing??
unique(fq$name)[!toupper(unique(fq$name)) %in% toupper(a$name)]


