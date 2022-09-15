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

parentsp=rPHG::graphBuilder(configFile=configPath, methods='mummer4_NAM_released_assemblies')
parents=t(SummarizedExperiment::assays(parentsp)$hapID)
parents[parents==-1]=NA
#write.table(parents, 'phg_hapids_parents_dbout.txt', row.names=T, col.names=T, sep='\t', quote=F)
#all.haps=read.table('phg_hapids_parents_dbout.txt', header=T)
### okay, i know this is horrible, but renaming parents to all.haps here because I don't want to rewrite all fo this
## all.haps=parents


## this is the values for each haplotype of how many TEs there are
atorrm=read.table('allNAM_hapids.TEbpUpdate.sup.2022-08-30.txt', header=T, comment.char='')

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
allknobmat=sapply(1:ncol(all.haps), function(x) atorrm$knobbp[match(all.haps[,x], atorrm$hapid)])
knob180mat=sapply(1:ncol(all.haps), function(x) atorrm$knob180bp[match(all.haps[,x], atorrm$hapid)])
tr1mat=sapply(1:ncol(all.haps), function(x) atorrm$tr1bp[match(all.haps[,x], atorrm$hapid)])
centromeremat=sapply(1:ncol(all.haps), function(x) atorrm$centromerebp[match(all.haps[,x], atorrm$hapid)])
telomeremat=sapply(1:ncol(all.haps), function(x) atorrm$telomerebp[match(all.haps[,x], atorrm$hapid)])
ribosomalmat=sapply(1:ncol(all.haps), function(x) atorrm$ribosomalbp[match(all.haps[,x], atorrm$hapid)])

b73mat=sapply(1:ncol(all.haps), function(x) all.haps[,x] %in% atorrm$hapid[atorrm$genotype=='B73'])
## can get this one from subsetting gsmat by b73 mat, I HOPE
#b73bpmat=sapply(1:ncol(all.haps), function(x) all.haps[,x] %in% atorrm$hapid[atorrm$genotype=='B73')
                                                                           



rownames(gsmat)=rownames(all.haps)
rownames(temat)=rownames(all.haps)
rownames(nontemat)=rownames(all.haps)
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
rownames(allknobmat)=rownames(all.haps)
rownames(knob180mat)=rownames(all.haps)
rownames(tr1mat)=rownames(all.haps)
rownames(centromeremat)=rownames(all.haps)
rownames(telomeremat)=rownames(all.haps)
rownames(ribosomalmat)=rownames(all.haps)
rownames(b73mat)=rownames(all.haps)

colnames(gsmat)=colnames(all.haps)
colnames(temat)=colnames(all.haps)
colnames(nontemat)=colnames(all.haps)
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
colnames(allknobmat)=colnames(all.haps)
colnames(knob180mat)=colnames(all.haps)
colnames(tr1mat)=colnames(all.haps)
colnames(centromeremat)=colnames(all.haps)
colnames(telomeremat)=colnames(all.haps)
colnames(ribosomalmat)=colnames(all.haps)
colnames(b73mat)=colnames(all.haps)

b73bpmat=gsmat
b73bpmat[!b73mat]=0
                    
matrixList=list(gsmat, temat, nontemat, dhhmat, dtamat, dtcmat, dthmat, dtmmat, dttmat, rilmat, ritmat, rlcmat, rlgmat, rlxmat, allknobmat, knob180mat, tr1mat, centromeremat, telomeremat, ribosomalmat, b73mat, b73bpmat)
names(matrixList)=c('gsmat', 'temat', 'nontemat', 'dhhmat', 'dtamat', 'dtcmat', 'dthmat', 'dtmmat', 'dttmat', 'rilmat', 'ritmat', 'rlcmat', 'rlgmat', 'rlxmat','allknobmat', 'knob180mat', 'tr1mat', 'centromeremat', 'telomeremat', 'ribosomalmat', 'b73mat', 'b73bpmat')

                    
saveRDS(matrixList, paste0('matrixList_bpEachRefRange.', Sys.Date(), '.RDS'))

   
                    
## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotyped.2022-03-22.txt', header=T)

                    
rilsum=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat, na.rm=T), 
                  tebp=rowSums(temat, na.rm=T),
                 nontebp=rowSums(nontemat, na.rm=T),
                 dhhbp=rowSums(dhhmat, na.rm=T),
                 dtabp=rowSums(dtamat, na.rm=T),
                 dtcbp=rowSums(dtcmat, na.rm=T),
                 dthbp=rowSums(dthmat, na.rm=T),
                 dtmbp=rowSums(dtmmat, na.rm=T),
                 dttbp=rowSums(dttmat, na.rm=T),
                 rilbp=rowSums(rilmat, na.rm=T),
                 ritbp=rowSums(ritmat, na.rm=T),
                 rlcbp=rowSums(rlcmat, na.rm=T),
                 rlgbp=rowSums(rlgmat, na.rm=T),
                 rlxbp=rowSums(rlxmat, na.rm=T),
                 allknobbp=rowSums(allknobmat, na.rm=T),
                 knob180bp=rowSums(knob180mat, na.rm=T),
                 tr1bp=rowSums(tr1mat, na.rm=T),
                 centromerebp=rowSums(centromeremat, na.rm=T),
                 telomerebp=rowSums(telomeremat, na.rm=T),
                 ribosomalbp=rowSums(ribosomalmat, na.rm=T))#,
 #                b73bp=rowSums(b73mat, na.rm=T))

                    
for(minB73correct in c(0,5,10,20,22,37)){

KEEPb73=b$refrange[b$nB73correct>minB73correct]
rilsumKeep=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T), 
                      tebp=rowSums(temat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     nontebp=rowSums(nontemat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dhhbp=rowSums(dhhmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dtabp=rowSums(dtamat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dtcbp=rowSums(dtcmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dthbp=rowSums(dthmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dtmbp=rowSums(dtmmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     dttbp=rowSums(dttmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     rilbp=rowSums(rilmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     ritbp=rowSums(ritmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     rlcbp=rowSums(rlcmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     rlgbp=rowSums(rlgmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     rlxbp=rowSums(rlxmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     allknobbp=rowSums(allknobmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     knob180bp=rowSums(knob180mat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     tr1bp=rowSums(tr1mat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     centromerebp=rowSums(centromeremat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     telomerebp=rowSums(telomeremat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     ribosomalbp=rowSums(ribosomalmat[,colnames(gsmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     b73bp=rowSums(b73bpmat[,colnames(b73bpmat)%in% paste0('X',KEEPb73)], na.rm=T),
                     b73rr=rowSums(b73mat[,colnames(b73mat) %in% paste0('X', KEEPb73)],na.rm=T))
            print(paste0('so now done with ', minB73correct))
            head(rilsumKeep)
write.table(rilsumKeep, paste0('ril_bp_repeats.greaterthan', minB73correct, 'B73correct.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
}                    
                    
     parentsumKeep=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T), 
                      tebp=rowSums(temat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     nontebp=rowSums(nontemat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dhhbp=rowSums(dhhmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dtabp=rowSums(dtamat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dtcbp=rowSums(dtcmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dthbp=rowSums(dthmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dtmbp=rowSums(dtmmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     dttbp=rowSums(dttmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     rilbp=rowSums(rilmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     ritbp=rowSums(ritmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     rlcbp=rowSums(rlcmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     rlgbp=rowSums(rlgmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     rlxbp=rowSums(rlxmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     allknobbp=rowSums(allknobmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     knob180bp=rowSums(knob180mat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     tr1bp=rowSums(tr1mat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     centromerebp=rowSums(centromeremat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     telomerebp=rowSums(telomeremat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     ribosomalbp=rowSums(ribosomalmat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T),
                     b73bp=rowSums(b73mat[,colnames(gsmat)%in% paste0('R',b$refrange[b$KEEPfinalFilter])], na.rm=T))
               
write.table(parentsumKeep, paste0('parent_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
 
