library(rtracklayer)
library(vcfR)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)
library(data.table)
library(dplyr)
library(stringr)



a=fread('../segregation_distortion/allNAM_hapids.TEbp.sup.2022-03-17.txt')

rs= a %>% group_by(refrange) %>% summarize(hapbprange=max(hapwidth, na.rm=T)-min(hapwidth,na.rm=T), tebprange=max(TEbp, na.rm=T)-min(TEbp, na.rm=T))
rs$filter500kb=rs$hapbprange>500000
rs$filter1Mb=rs$hapbprange>1e6

rs$KEEPlowrange=!rs$filter1Mb

### filter based on b73 genotyping
b73correct=read.table('refranges_B73correctlygenotyped.2022-03-22.txt', header=T)

tmp=merge(rs, b73correct, keep='all')

rr=read.csv('maize1_reference_ranges.csv')
rr$genic=rr$name=='FocusRegion' 

##  number of refranges we keep
sum(tmp$KEEPlowrange & tmp$KEEPb73correct)

## number of refranges we remove
sum(!tmp$KEEPlowrange | !tmp$KEEPb73correct)

## amount of genome we are considering
sum(rr[rr$ref_range_id %in% tmp$refrange[tmp$KEEPb73correct & tmp$KEEPlowrange],]$range_end+1-rr[rr$ref_range_id %in% tmp$refrange[tmp$KEEPb73correct & tmp$KEEPlowrange],]$range_start)



tmp$KEEPfinalFilter=tmp$KEEPlowrange&tmp$KEEPb73correct


write.table(tmp[,c('refrange', 'KEEPfinalFilter')], paste0('refranges_B73correctlygenotypedAND1Mbrangeremoved.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
