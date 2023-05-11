library(rtracklayer)
library(dplyr)

## read in haplotypes for each ril
all.haps=read.table('../phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps[all.haps==-1]=NA

## for B73, find refranges that overlap longest sv
## then make a table of which parent they have in this region
## ideally, the majority of refranges will have one parent or the other (not expecting recombination in the sv)

lsv=read.table('longest_sv_per_nam.txt', header=T)
rr=read.csv('../maize1_reference_ranges.csv', header=T)

lsvGR=GRanges(seqnames=gsub('chr', '', lsv$V1), IRanges(start=lsv$V2, end=lsv$V3))
rrGR=GRanges(seqnames=rr$chrom, IRanges(start=rr$range_start, end=rr$range_end))

ov=findOverlaps(lsvGR, rrGR)

## for each ril
##   for its lsv
##      how many refranges have B73 vs how many have other geno as parent

## which haps go with which rr's
atorrm=read.table('../allNAM_hapids.TEbpUpdate.sup.2022-08-30.txt', header=T, comment.char='')


lsvHaps=cbind(lsv[,c(1,2,3,7,8)][queryHits(ov),], refrange=rr$ref_range_id[subjectHits(ov)])
lsvHaps$B73HapID=atorrm[atorrm$genotype=='B73',]$hapid[match(lsvHaps$refrange, atorrm[atorrm$genotype=='B73',]$refrange)]
### sort of slow but it's okay
lsvHaps$OtherParentHapID=sapply(1:nrow(lsvHaps), function(x) atorrm$hapid[atorrm$genotype==lsvHaps$genome[x] & atorrm$refrange==lsvHaps$refrange[x]])

## icky this is not working when the same region of the genome is sv across diff nams
#atorrm$SVB73genotype=lsvHaps$genome[match(atorrm$hapid, lsvHaps$B73HapID) ]
#atorrm$SVOthergenotype=lsvHaps$genome[match(atorrm$hapid, lsvHaps$OtherParentHapID) ]
                
lsvHaps$B73TEbp=atorrm$TEbp[match(lsvHaps$B73HapID, atorrm$hapid)]
lsvHaps$OtherTEbp=atorrm$TEbp[match(lsvHaps$OtherParentHapID, atorrm$hapid)]
lsvHaps$B73Size=atorrm$hapwidth[match(lsvHaps$B73HapID, atorrm$hapid)]
lsvHaps$OtherSize=atorrm$hapwidth[match(lsvHaps$OtherParentHapID, atorrm$hapid)]
     
SVsummary=lsvHaps %>% group_by(genome) %>% dplyr::summarize(B73TEbp=sum(B73TEbp), OtherTEbp=sum(OtherTEbp, na.rm=T), B73Size=sum(B73Size), OtherSize=sum(OtherSize, na.rm=T))
          
## two are bigger
SVsummary[SVsummary$B73Size>SVsummary$OtherSize,]
                            
#b73allele=atorrm %>% group_by(SVB73genotype) %>% dplyr::summarize(TEbp=sum(TEbp), size=sum(hapwidth))                                
#otherallele=atorrm %>% group_by(SVOthergenotype) %>% dplyr::summarize(TEbp=sum(TEbp), size=sum(hapwidth))                                
  ### Ohno this is not good - there are negative nonTE bp because there are more TEs than bp in the region....
  ### back to the drawing board, recheck imputation again...
                                
                                
## now, for each RIL, get proportion of SV that is each parent
namfams=read.table('../../figures/nam_fams.txt') 
svGenos=data.frame(ril=rownames(all.haps))
svGenos$zfam=namfams$V2[match(substr(svGenos$ril,1,4), namfams$V1)]
## simplify the problem!     note some of these sv were pulled out from genotyping because differences in hapsizes were too big
sv.haps=all.haps[,colnames(all.haps)[colnames(all.haps)%in%paste0('X', lsvHaps$refrange)]]
svGenos$countB73=sapply(1:nrow(svGenos), function(x) sum(lsvHaps[lsvHaps$genome==svGenos$zfam[x],]$B73HapID %in% unlist(sv.haps[x,])))
svGenos$countNonB73=sapply(1:nrow(svGenos), function(x) sum(lsvHaps[lsvHaps$genome==svGenos$zfam[x],]$OtherParentHapID %in% unlist(sv.haps[x,])))

svGenos$svGeno=ifelse(svGenos$countB73+svGenos$countNonB73>1 & (svGenos$propB73 > 0.95 | svGenos$propB73 < 0.05) & svGenos$propB73>0.95, 'B73', 'placeholder')
svGenos$svGeno[svGenos$svGeno=='placeholder']=svGenos$zfam[svGenos$svGeno=='placeholder']

                           
                           
write.table(svGenos, 'longest_sv_RIL_genotypes.txt', row.names=F, col.names=T, sep='\t', quote=F)
                           
