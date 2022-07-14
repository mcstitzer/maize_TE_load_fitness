## change java memory before starting java

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx48g", "-Xms48g"))
gc()
#options(java.parameters = c("-Xmx48g", "-Xms48g"))



library(rtracklayer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(dplyr)
library(stringr)



library(rPHG)

rPHG::startLogger(fullPath = NULL, fileName = NULL)

###### if brapi server is working, do this
myCon <- rPHG::BrapiCon("cbsudc01.biohpc.cornell.edu")
myCon

myCon %>% rPHG::availablePHGMethods() ## need to pick a method

myMethod <- "NonMergedReadMapping_AllNamParents_Haploid" ## this is unideal, but diploid mappings against just the 2 parents aren't available yet

samples=myCon %>% 
  rPHG::PHGMethod(myMethod) %>%
  rPHG::readSamples()

refranges=myCon %>% 
  rPHG::PHGMethod(myMethod) %>%
  rPHG::readRefRanges()

## eh let's just try to get them all because why not...
all.haps=myCon %>% 
  rPHG::PHGMethod(myMethod) %>%
  rPHG::readTable(index = F) 
all.haps[all.haps==-1]=NA

#################################
### if brapi server isnt' working, do this on a cbsu machine or cornell vpn

configPath='phg_config.txt' ## has a password so not putting on github...
pathMethod=rPHG::pathsForMethod(configFile=configPath, pathMethod=myMethod)
dim(pathMethod)
#write.table(pathMethod, 'phg_hapids_namrils_dbout.txt', row.names=T, col.names=T, sep='\t', quote=F)
#all.haps=read.table('phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps=pathMethod
all.haps[all.haps==-1]=NA



## this is the values for each haplotype of how many TEs there are
atorrm=read.table('allNAM_hapids.TEbpUpdate.sup.2022-07-12.txt', header=T, comment.char='')

gsmat=sapply(1:ncol(all.haps), function(x) atorrm$hapwidth[match(all.haps[,x], atorrm$hapid)])
temat=sapply(1:ncol(all.haps), function(x) atorrm$TEbp[match(all.haps[,x], atorrm$hapid)])
nontemat=sapply(1:ncol(all.haps), function(x) atorrm$nonTEbp[match(all.haps[,x], atorrm$hapid)])
dhhmat=sapply(1:ncol(all.haps), function(x) atorrm$DHHbp[match(all.haps[,x], atorrm$hapid)])
dtamat=sapply(1:ncol(all.haps), function(x) atorrm$DTAbp[match(all.haps[,x], atorrm$hapid)])
dtcmat=sapply(1:ncol(all.haps), function(x) atorrm$DTCbp[match(all.haps[,x], atorrm$hapid)])
dthmat=sapply(1:ncol(all.haps), function(x) atorrm$DTHbp[match(all.haps[,x], atorrm$hapid)])
dtmmat=sapply(1:ncol(all.haps), function(x) atorrm$DTMbp[match(all.haps[,x], atorrm$hapid)])
dttmat=sapply(1:ncol(all.haps), function(x) atorrm$DTTbp[match(all.haps[,x], atorrm$hapid)])
rilmat=sapply(1:ncol(all.haps), function(x) atorrm$RILbp[match(all.haps[,x], atorrm$hapid)])
ritmat=sapply(1:ncol(all.haps), function(x) atorrm$RITbp[match(all.haps[,x], atorrm$hapid)])
rlcmat=sapply(1:ncol(all.haps), function(x) atorrm$RLCbp[match(all.haps[,x], atorrm$hapid)])
rlgmat=sapply(1:ncol(all.haps), function(x) atorrm$RLGbp[match(all.haps[,x], atorrm$hapid)])
rlxmat=sapply(1:ncol(all.haps), function(x) atorrm$RLXbp[match(all.haps[,x], atorrm$hapid)])
knob180mat=sapply(1:ncol(all.haps), function(x) atorrm$knob180bp[match(all.haps[,x], atorrm$hapid)])
tr1mat=sapply(1:ncol(all.haps), function(x) atorrm$tr1bp[match(all.haps[,x], atorrm$hapid)])
centromeremat=sapply(1:ncol(all.haps), function(x) atorrm$centromerebp[match(all.haps[,x], atorrm$hapid)])
telomeremat=sapply(1:ncol(all.haps), function(x) atorrm$telomerebp[match(all.haps[,x], atorrm$hapid)])
ribosomalmat=sapply(1:ncol(all.haps), function(x) atorrm$ribosomalbp[match(all.haps[,x], atorrm$hapid)])




rownames(gsmat)=rownames(all.haps)
rownames(temat)=rownames(all.haps)
rownames(dhhmat)=rownames(all.haps)
rownames(dtamat)=rownames(all.haps)
rownames(dtcmat)=rownames(all.haps)
rownames(dthmat)=rownames(all.haps)
rownames(dtmmat)=rownames(all.haps)
rownames(dttmat)=rownames(all.haps)
rownames(rilmat)=rownames(all.haps)
rownames(ritmat)=rownames(all.haps)
rownames(rlcmat)=rownames(all.haps)
rownames(rlgmat)=rownames(all.haps)
rownames(rlxmat)=rownames(all.haps)
rownames(knobmat)=rownames(all.haps)
rownames(centromeremat)=rownames(all.haps)
rownames(telomeremat)=rownames(all.haps)
rownames(ribosomalmat)=rownames(all.haps)

colnames(gsmat)=colnames(all.haps)
colnames(temat)=colnames(all.haps)
colnames(dhhmat)=colnames(all.haps)
colnames(dtamat)=colnames(all.haps)
colnames(dtcmat)=colnames(all.haps)
colnames(dthmat)=colnames(all.haps)
colnames(dtmmat)=colnames(all.haps)
colnames(dttmat)=colnames(all.haps)
colnames(rilmat)=colnames(all.haps)
colnames(ritmat)=colnames(all.haps)
colnames(rlcmat)=colnames(all.haps)
colnames(rlgmat)=colnames(all.haps)
colnames(rlxmat)=colnames(all.haps)
colnames(knobmat)=colnames(all.haps)
colnames(centromeremat)=colnames(all.haps)
colnames(telomeremat)=colnames(all.haps)
colnames(ribosomalmat)=colnames(all.haps)

                    
matrixList=list(gsmat, temat, nontemat, dhhmat, dtamat, dtcmat, dthmat, dtmmat, dttmat, rilmat, ritmat, rlcmat, rlgmat, rlxmat, knob180mat, tr1mat, centromeremat, telomeremat, ribosomalmat)
names(matrixList)=c('gsmat', 'temat', 'nontemat', 'dhhmat', 'dtamat', 'dtcmat', 'dthmat', 'dtmmat', 'dttmat', 'rilmat', 'ritmat', 'rlcmat', 'rlgmat', 'rlxmat', 'knob180mat', 'tr1mat', 'centromeremat', 'telomeremat', 'ribosomalmat')

                    
saveRDS(matrixList, 'matrixList_bpEachRefRange.RDS')

   
                    
## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt', header=T)

                    
rilsum=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat, na.rm=T), tebp=rowSums(temat, na.rm=T))

rilsumKeep=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat[,colnames(gsmat)%in% paste0('X',b$refrange[b$KEEPfinalFilter])], na.rm=T), tebp=rowSums(temat[,colnames(gsmat)%in% paste0('X',b$refrange[b$KEEPfinalFilter])], na.rm=T))
            
write.table(rilsumKeep, paste0('ril_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
                    
                    
                    
