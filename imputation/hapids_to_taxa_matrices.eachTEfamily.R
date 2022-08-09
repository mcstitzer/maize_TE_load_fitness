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


matchMatrix=data.frame(lapply(1:ncol(all.haps), function(refrange) match(all.haps[,refrange], atorrm$hapid)))
colnames(matchMatrix)=colnames(all.haps)


famMatrices=lapply(colnames(atorrm)[13:ncol(atorrm)], function(tefam){
 tempMat=sapply(1:ncol(all.haps), function(refrange) atorrm[matchMatrix[,refrange],tefam])
 rownames(tempMat)=rownames(all.haps)
 colnames(tempMat)=colnames(all.haps)
 return(tempMat)
})
names(famMatrices)=colnames(atorrm)[13:ncol(atorrm)]
#### wowowowow this might actually work!!!!


                    
saveRDS(famMatrices, 'TEfamMatrixList_bpEachRefRange.RDS')

   
                    
## get rid of b73 flavors and bad ranges
b=read.table('refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt', header=T)

rilsumKeep=cbind(id=rownames(all.haps), data.frame(lapply(famMatrices, function(x) rowSums(x[,colnames(x) %in% paste0('X', b$refrange[b$KEEPfinalFilter])], na.rm=T)))
    
write.table(rilsumKeep, paste0('ril_TEFam_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')

### okay, i know this is horrible, but renaming parents to all.haps here because I don't want to rewrite all fo this
all.haps=parents

famMatrices=lapply(colnames(atorrm)[X:ncol(atorrm)], function(x){
 tempMat=sapply(1:ncol(all.haps), function(refrange) atorrm[,x][match(all.haps[,refrange], atorrm$hapid)])
 rownames(tempMat)=rownames(all.haps)
 colnames(tempMat)=colnames(all.haps)
 return(tempMat)
})
names(famMatrices)=colnames(atorrm)[X:ncol(atorrm)]

parentsumKeep=cbind(id=rownames(all.haps), data.frame(lapply(famMatrices, function(x) rowSums(x[,colnames(x) %in% paste0('X', b$refrange[b$KEEPfinalFilter])], na.rm=T)))
                    
       
write.table(parentsumKeep, paste0('parent_TEFam_bp_repeats.', Sys.Date(), '.txt'), row.names=F, col.names=T, quote=F, sep='\t')
 
