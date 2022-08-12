library(rtracklayer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(dplyr)
library(stringr)



all.haps=read.table('phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps[all.haps==-1]=NA


parents=read.table('phg_hapids_parents_dbout.txt', header=T)


## this is the values for each haplotype of how many TEs there are
atorrm=read.table('allNAM_hapids.FamiliesUpdate.sup.2022-08-09.txt', header=T, comment.char='')

## this part takes a long long time (2 hours??) - so it's good i've pulled it out from the below lapply :~)
matchMatrix=data.frame(lapply(1:ncol(all.haps), function(refrange) match(all.haps[,refrange], atorrm$hapid)))
colnames(matchMatrix)=colnames(all.haps)

### make a matrix and save separately for each!
dir.create('TEfamily_matrices/')
                              
tefams=data.frame(RIL=rownames(all.haps))
                              
## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt', header=T)

keepColnames=colnames(all.haps)%in% paste0('X',b$refrange[b$KEEPfinalFilter])
                              
for(tefam in colnames(atorrm)[13:ncol(atorrm)]){
 tempMat=sapply(1:ncol(all.haps), function(refrange) atorrm[matchMatrix[,refrange],tefam])
 rownames(tempMat)=rownames(all.haps)
 colnames(tempMat)=colnames(all.haps)
 tefams[,tefam]=rowSums(tempMat[,keepColnames], na.rm=T)
  ##writing output takes forever let's just get this done!!!
# write.table(tempMat, paste0('TEfamily_matrices/',tefam, '.bp.txt'), row.names=T, col.names=T, quote=F, sep='\t')
}
#names(famMatrices)=colnames(atorrm)[13:ncol(atorrm)]
#### wowowowow this might actually work!!!!
 ## lol saving in an array took 500 Gb so I changed to save each to disk separately :)


                    

write.table(tefams, paste0('ril_TEFam_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')

### okay, i know this is horrible, but renaming parents to all.haps here because I don't want to rewrite all fo this
all.haps=parents
matchMatrix=data.frame(lapply(1:ncol(all.haps), function(refrange) match(all.haps[,refrange], paste0('R',atorrm$hapid))))
colnames(matchMatrix)=colnames(all.haps)
                                  
tefams=data.frame(RIL=rownames(all.haps))
                                               
for(tefam in colnames(atorrm)[13:ncol(atorrm)]){
 tempMat=sapply(1:ncol(all.haps), function(refrange) atorrm[matchMatrix[,refrange],tefam])
 rownames(tempMat)=rownames(all.haps)
 colnames(tempMat)=colnames(all.haps)
 tefams[,tefam]=rowSums(tempMat[,keepColnames], na.rm=T)
  ##writing output takes forever let's just get this done!!!
# write.table(tempMat, paste0('TEfamily_matrices/',tefam, '.bp.txt'), row.names=T, col.names=T, quote=F, sep='\t')
}

 
write.table(tefams, paste0('parent_TEFam_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
 
