library(rtracklayer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(dplyr)
library(stringr)
library(reshape2)


all.haps=read.table('phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps[all.haps==-1]=NA


parents=read.table('phg_hapids_parents_dbout.txt', header=T)


## this is the values for each haplotype of how many TEs there are

## this is the values for each haplotype of how many TEs there are
atorrm=read.table('coregene_distance_bins.txt', header=T, comment.char='')
## reshape to wid eformat so we only count each hapid once
atorrm=dcast(atorrm, hapid~genecat)
atorrm[is.na(atorrm)]=0 ## ncie and this mkes sense because ranges alternate between genic intergenic and 0 and value ingene

## add columns for umr
a2=read.table('umr_bins.txt', header=T, comment.char='')
a2=dcast(a2, hapid~umrCountnonzero)
a2[is.na(a2)]=0 ## ncie and this mkes sense because ranges alternate between genic intergenic and 0 and value ingene
colnames(a2)[2:3]=c('noUMR', 'hasUMR')

## add columns for recentinsertions
a3=read.table('recentinsertion_bins.txt', header=T, comment.char='')
a3=dcast(a3, hapid~recentinsertion)
a3[is.na(a3)]=0 ## ncie and this mkes sense because ranges alternate between genic intergenic and 0 and value ingene
colnames(a3)[2:3]=c('olderInsertion', 'recentInsertion')


## this part takes a long long time (2 hours??) - so it's good i've pulled it out from the below lapply :~)
matchMatrix=data.frame(lapply(1:ncol(all.haps), function(refrange) match(all.haps[,refrange], atorrm$hapid)))
colnames(matchMatrix)=colnames(all.haps)

## okay leaving for a run

atorrm=merge(atorrm, a2)
atorrm=merge(atorrm, a3)
                              
## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt', header=T)

keepColnames=colnames(all.haps)%in% paste0('X',b$refrange[b$KEEPfinalFilter])
          
          
                            
tefams=data.frame(RIL=rownames(all.haps))

for(tefam in colnames(atorrm)[2:ncol(atorrm)]){
 tempMat=sapply(1:ncol(all.haps), function(refrange) atorrm[matchMatrix[,refrange],tefam])
 rownames(tempMat)=rownames(all.haps)
 colnames(tempMat)=colnames(all.haps)
 tefams[,tefam]=rowSums(tempMat[,keepColnames], na.rm=T)
}

                    

write.table(tefams, paste0('ril_tecategories_bp.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
