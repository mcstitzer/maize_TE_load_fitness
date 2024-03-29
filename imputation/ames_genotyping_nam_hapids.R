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




## use Ames genotyping to pull out B73 replicates
amesMethod <- "Ames_NonMergedReadMapping_AllLines_Haploid" ## there are several B73's - this will hopefully ask if they show differences from refgen
nam=read.table('../figures/nam_fams.txt')

myCon <- rPHG::BrapiCon("cbsudc01.biohpc.cornell.edu")
myCon

myCon %>% rPHG::availablePHGMethods() ## need to pick a method


## pull out all NAM parents
samples=myCon %>% 
    rPHG::PHGMethod(amesMethod) %>%
    rPHG::readSamples() %>%
    filter(str_split_fixed(toupper(sampleName), '-', 2)[,1] %in% toupper(c('B73',nam$V2[-17]))) ## don't want Mo17

################################
### if brapi server isnt' working, do this on a cbsu machine or cornell vpn

## this gets all AMES, then subsample to my nam parent taxa

configPath='phg_config.txt' ## has a password so not putting on github...
pathMethod=rPHG::pathsForMethod(configFile=configPath, pathMethod=amesMethod)
dim(pathMethod)
#write.table(pathMethod, 'phg_hapids_namrils_dbout.txt', row.names=T, col.names=T, sep='\t', quote=F)
#all.haps=read.table('phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps=pathMethod
all.haps[all.haps==-1]=NA

ames.haps=all.haps[rownames(all.haps) %in% samples$sampleName,]

## eh, easier than renaming the next steps :(
all.haps=ames.haps


## this is the values for each haplotype of how many TEs there are
atorrm=read.table('allNAM_hapids.TEbpUpdate.sup.2022-08-30.txt', header=T, comment.char='')

all.haps.namname=str_split_fixed(rownames(all.haps), '-', 2)[,1]
## count per sample how many "self" refranges are called
selfmat=sapply(1:ncol(all.haps), function(x) toupper(atorrm$genotype[match(all.haps[,x], atorrm$hapid)])==toupper(all.haps.namname))

               
colnames(selfmat)=colnames(all.haps)
rownames(selfmat)=rownames(all.haps)
selfrr=rowSums(selfmat[,colnames(selfmat) %in% paste0(KEEPb73)],na.rm=T)

sr=data.frame(selfrr)
dim(sr)
head(sr)
sr$nam=all.haps.namname

               
rr=read.csv('../imputation/maize1_reference_ranges.csv')

               
               
rr %>% group_by(chrom) %>% summarize(max(ref_range_id))
cm=rr %>% group_by(chrom) %>% summarize(max(ref_range_id))
pdf('~/transfer/ames_refrange_count.pdf', 24,38)
ggplot(sr, aes(x=nam, group=nam, y=selfrr)) + geom_boxplot()
ggplot(sr, aes(x=nam, label=nam, y=selfrr)) + geom_text()
ggplot(melt(selfmat), aes(x=Var2, color=value, fill=value, y=Var1)) + geom_tile() + scale_fill_brewer(palette='Set
1', na.value='white') + scale_color_brewer(palette='Set1', na.value='white') + geom_vline(
head(cm)
colnames(cm)[2]='maxrr'
ggplot(melt(selfmat), aes(x=Var2, color=value, fill=value, y=Var1)) + geom_tile() + scale_fill_brewer(palette='Set1', na.value='white') + scale_color_brewer(palette='Set1', na.value='white') + geom_vline(x_intercept=cm$maxrr)
ggplot(melt(selfmat), aes(x=Var2, color=value, fill=value, y=Var1)) + geom_tile() + scale_fill_brewer(palette='Set1', na.value='white') + scale_color_brewer(palette='Set1', na.value='white') + geom_vline(x=intercept=cm$maxrr)
ggplot(melt(selfmat), aes(x=Var2, color=value, fill=value, y=Var1)) + geom_tile() + scale_fill_brewer(palette='Set1', na.value='white') + scale_color_brewer(palette='Set1', na.value='white') + geom_vline(xintercept=cm$maxrr)
dev.off()
               
               
               
rc=read.table('~/transfer/amesGBS_linecount_sodivthisby4.txt', header=F)
rc$reads=rc$V1/4
rc$name=gsub('amesDemultiplex_', '', gsub('.fastq', '', rc$V2))
temp=str_split_fixed(rc$name, ':', 2)
rc$nname=paste0(temp[,1], '-', temp[,2])
rc$nnname=str_split_fixed(rc$nname, '_', 2)[,1]
rc$nnnname=gsub(':', '_', rc$nnname)
sr$readcount=rc$reads[match(rownames(sr), rc$nnnname)]
        
pdf('~/transfer/ames_readcount_coverage.pdf', 14,8)
ggplot(sr, aes(x=readcount, y=selfrr)) + geom_text(aes(label=nam))
ggplot(sr, aes(x=readcount, y=selfrr)) + geom_text(aes(label=nam))
ggplot(sr, aes(x=readcount, y=selfrr, color=nam=='B73')) + geom_text(aes(label=nam, alpha=0.8)) + scale_color_manual(values=c('black', 'red')) + xlab('GBS read count') + ylab('Refrange called with own haplotype')
dev.off()
 ## this history is saved in ohnothisisbad.R
               
               

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
                    
rownames(selfmat)=rownames(all.haps)

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

colnames(selfmat)=colnames(all.haps)
             
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


## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotyped.2022-03-22.txt', header=T)
minB73correct=22
KEEPb73=b$refrange[b$nB73correct>minB73correct]


rilsumKeep=data.frame(id=rownames(all.haps), genomesize=rowSums(gsmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T), 
                      tebp=rowSums(temat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     nontebp=rowSums(nontemat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dhhbp=rowSums(dhhmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dtabp=rowSums(dtamat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dtcbp=rowSums(dtcmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dthbp=rowSums(dthmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dtmbp=rowSums(dtmmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     dttbp=rowSums(dttmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     rilbp=rowSums(rilmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     ritbp=rowSums(ritmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     rlcbp=rowSums(rlcmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     rlgbp=rowSums(rlgmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     rlxbp=rowSums(rlxmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     allknobbp=rowSums(allknobmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     knob180bp=rowSums(knob180mat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     tr1bp=rowSums(tr1mat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     centromerebp=rowSums(centromeremat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     telomerebp=rowSums(telomeremat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     ribosomalbp=rowSums(ribosomalmat[,colnames(gsmat)%in% paste0(KEEPb73)], na.rm=T),
                     b73bp=rowSums(b73bpmat[,colnames(b73bpmat)%in% paste0(KEEPb73)], na.rm=T),
                     b73rr=rowSums(b73mat[,colnames(b73mat) %in% paste0(KEEPb73)],na.rm=T),
                     selfrr=rowSums(selfmat[,colnames(selfmat) %in% paste0(KEEPb73)],na.rm=T))
            print(paste0('so now done with ', minB73correct))
            head(rilsumKeep)
              rilsumKeep$nam=str_split_fixed(rilsumKeep$id, '-', 2)[,1]
write.table(rilsumKeep, paste0('ames_bp_repeats.greaterthan', minB73correct, 'B73correct.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
