library(stringr)
library(data.table)
library(rtracklayer)
library(assertthat)
library(dplyr)


##################################
## Get ranges of each haplotype ##
##################################

## header lines from hvcf file, to get assembly coordinates for each haplotype id
## scp cbsufeschotte2:/workdir/mcs368/nam_te_pav/imputation/nam_unmerged_haplotypes.haplines .

a=fread('nam_unmerged_haplotypes.haplines')
colnames(a)=c('V1', 'V2')

a$hapid=str_split_fixed(a$V1, '=', 3)[,3]
a$genotype=gsub('Description="', '', str_split_fixed(a$V2, '_', 2)[,1])
a$chr=str_split_fixed( a$V2, ':', 3)[,2]
a$range=gsub('">', '', str_split_fixed( a$V2, ':', 3)[,3])
a$start=as.numeric(str_split_fixed(a$range, '-', 2)[,1])
a$end=as.numeric(str_split_fixed(a$range, '-', 2)[,2])
a$startmin=sapply(1:nrow(a), function(x) min(a$start[x], a$end[x]))
a$endmax=sapply(1:nrow(a), function(x) max(a$start[x], a$end[x]))
a$hapwidth=a$endmax-a$startmin

## switch capitalization so it works with my stuff
a$genotype[a$genotype=='IL14H']='Il14H'

## read in refrange id so that it's easy to combine with other files
b=read.table('hapids_to_refranges.txt', header=T, comment.char='')
a$refrange=b$refrange[match(a$hapid, b$hapid)]


######################################################
## Read in TE annotation to count TEs per haplotype ##
######################################################

### do this across all NAM
genomes=c('B73', 'CML103', 'CML247', 'CML52', 'CML69', 'Ky21', 'M162W', 'Mo18W', 'Oh43', 'P39', 
          'B97', 'CML228', 'CML277', 'CML322', 'CML333', 'HP301', 'Il14H', ## note here the capitalization of Il14H!!!!! phg uses IL14H
         'Ki11', 'Ki3', 'M37W', 'Ms71','NC350', 'NC358', 'Oh7B', 'Tx303', 'Tzi8')

              
## not all things in the TE file are TEs, so set up these categories
classificationTE=c('DNA/DTA', 'DNA/DTC', 'DNA/DTH', 'DNA/DTM', 'DNA/DTT', 'DNA/Helitron', 'LINE/L1', 'LINE/RTE', 'LINE/unknown', 'LTR/CRM', 'LTR/Copia', 'LTR/Gypsy', 'LTR/unknown', 'MITE/DTA', 'MITE/DTC', 'MITE/DTH', 'MITE/DTM', 'MITE/DTT')
classificationKnob=c('knob/TR-1', 'knob/knob180')
classificationCent=c('Cent/CentC')
classificationTelo=c('subtelomere/4-12-1')
classificationRibo=c('rDNA/spacer')
otherClassifications=list(classificationKnob, classificationCent, classificationTelo, classificationRibo)
names(otherClassifications)=c('knob', 'centromere', 'telomere', 'ribosomal')

tesup=c("DHH", "DTA", "DTC", "DTH", "DTM", "DTT", "RIL", "RIT", "RIX", "RLC", "RLG", "RLX")

specificRepeats=c('knob/knob180', 'knob/TR-1', 'LTR/CRM')
names(specificRepeats)=c('knob180', 'tr1', 'crm') ## keep separate knobs, centromeric retrotransposon annotation

## add numeric columns for each classification category we're going to consider
for(i in c(tesup, names(otherClassifications), names(specificRepeats))){
  a[,paste0(i, 'bp')]=0
  }

## get started!
a=data.frame(a) ## just in case it's a data table              
tes=vector(mode = "list", length = length(genomes)) ## this is a list of the te annotation of each genome
names(tes)=genomes


#######################################################3
###### OKAY, there's an updated version of EDT annotations for NAM


## loop through and populate the entry of tes
for(genome in genomes){
 ## read in EDTA TE annotation
 ### filenames have different capitalization than before :( so sad
 filegenome=genome
 if(genome == 'Oh7B'){filegenome='Oh7b'}
 edtafilename=list.files(path='../genomes_and_annotations/tes/new_edta_tes/', pattern=paste0(filegenome, '.+', 'split.gff3.gz'), full.names=T)
 te=import.gff(edtafilename)
 ## the new ones go B73_chr1
 seqlevels(te)=str_split_fixed(seqlevels(te), '_', 2)[,2]
 
 ## reformat to get superfamily for each TE annotation
 te$sup=str_split_fixed(te$Classification, '/', 2)[,2]
 te$sup[te$sup=='Gypsy']='RLG'
 te$sup[te$sup=='Copia']='RLC'
 te$sup[te$sup=='unknown' & te$Classification=='LTR/unknown']='RLX'
 te$sup[te$sup=='Helitron']='DHH'
 te$sup[te$sup=='CRM']='RLG' ## need these and following ones because not structural
 te$sup[te$sup=='L1']='RIL'
 te$sup[te$sup=='RTE']='RIT'
 te$sup[te$sup=='unknown' & te$Classification=='LINE/unknown']='RIX'
 
 te$sup[!te$Classification %in% classificationTE]=NA

## make sure we added right
  assert_that(sum(te$Classification %in% c('DNA/DTA', 'MITE/DTA'))==sum(te$sup=='DTA', na.rm=T))
  assert_that(sum(te$Classification %in% c('DNA/DTC', 'MITE/DTC'))==sum(te$sup=='DTC', na.rm=T))
  assert_that(sum(te$Classification %in% c('DNA/DTH', 'MITE/DTH'))==sum(te$sup=='DTH', na.rm=T))
  assert_that(sum(te$Classification %in% c('DNA/DTM', 'MITE/DTM'))==sum(te$sup=='DTM', na.rm=T))
  assert_that(sum(te$Classification %in% c('DNA/DTT', 'MITE/DTT'))==sum(te$sup=='DTT', na.rm=T))
  assert_that(sum(te$Classification %in% c('LTR/CRM', 'LTR/Gypsy'))==sum(te$sup=='RLG', na.rm=T))

## not an issue for these guys
## Oh7B TE annotation has "oh7b" on the chromosome names... just remove that
# if(genome=='Oh7B'){seqlevels(te)=gsub(tolower(paste0(genome, '_')), '', tolower(seqlevels(te)))}

## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
 
 ### testing! as a solution to structural variatns
          
 ## 1. flag any overlapping
 longhapov=findOverlaps(haps, drop.self=T)
 haps$overlap=1:length(haps) %in% unique(unlist(queryHits(longhapov), subjectHits(longhapov)))
 a$hapov[a$genotype==genome]=haps$overlap
          
longhaps=unique(queryHits(longhapov)[duplicated(queryHits(longhapov))])
 if(length(longhaps)>0){
    for(longhap in longhaps){ ### this should be robust to strand??? since these are all unstranded
        longhap=longhaps[1]
        insidehaps=haps[unique(subjectHits(longhapov)[queryHits(longhapov)==longhap]),]
        for(insidehap in 1:length(insidehaps)){
          ranges(haps[longhap,])=ranges(GenomicRanges::setdiff(haps[longhap,], insidehaps[insidehap,])[1,])
    }}}
          
# withinhaps=findOverlaps(haps, drop.self=T, type='within') ## these will be caught by the outer haplotype, and bp counted for them
# if(length(withinhaps)>0){
#          haps=haps[-queryHits(withinhaps),]
#           }

          

  ## two adjacent haplotype ranges that have a te spanning them don't get intersected with intersect(reduce(te), haps, 
#  ### so pintersect to collapse the sets of overlapping ranges
#  tehapintersect=GenomicRanges::intersect(reduce(te[te$Classification %in% classificationTE,]), haps, ignore.strand=T)
#  thro1=findOverlaps(tehapintersect, haps, ignore.strand=T) ## subset first and then count overlaps
#  ov=pintersect(tehapintersect[queryHits(thro1)],haps[subjectHits(thro1)], ignore.strand=T) ## pairwise, so won't collapse the 1bp apart adjacent hits
#  thro=findOverlaps(ov, haps, ignore.strand=T)
#  tebp=sapply(unique(subjectHits(thro)), function(x) sum(width(ov[queryHits(thro)[subjectHits(thro)==x],]))) ## for individual sup/fam, go through reduce with each sup/fam - this is relaly far away from reality of te inseertion
#  a$TEbp[a$genotype==genome]=0
#  a$TEbp[a$genotype==genome][unique(subjectHits(thro))]=tebp

          
### try this new findoverlappairs thing which is absolutely amazing
   tehapintersect=findOverlapPairs(  haps, reduce(te[te$Classification %in% classificationTE,]), ignore.strand=T)
   ov=pintersect(tehapintersect, ignore.strand=T)
   tebp= data.frame(ov) %>% group_by(hapid) %>% dplyr::summarize(out=sum(width))
   a$TEbp[a$genotype==genome]=0
   a$TEbp[a$genotype==genome & a$hapid %in% tebp$hapid]=tebp$out[match(a$hapid[a$genotype==genome& a$hapid %in% tebp$hapid], tebp$hapid)]

              
  ### by superfamily
  for(sup in tesup){
     # tehapintersect=GenomicRanges::intersect(reduce(te[te$sup==sup & !is.na(te$sup),]), haps, ignore.strand=T)
      tehapintersect=findOverlapPairs(haps, reduce(te[te$sup==sup & !is.na(te$sup),]), ignore.strand=T)
     # thro1=findOverlaps(tehapintersect, haps, ignore.strand=T) ## subset first and then count overlaps
      ov=pintersect(tehapintersect, ignore.strand=T)
     # ov=pintersect(tehapintersect[queryHits(thro1)],haps[subjectHits(thro1)], ignore.strand=T) ## pairwise, so won't collapse the 1bp apart adjacent hits
     # thro=findOverlaps(ov, haps, ignore.strand=T)
      tebp=data.frame(ov) %>% group_by(hapid) %>% dplyr::summarize(out=sum(width))
     # tebp=sapply(unique(subjectHits(thro)), function(x) sum(width(ov[queryHits(thro)[subjectHits(thro)==x],]))) ## for individual sup/fam, go through reduce with each sup/fam - this is relaly far away from reality of te inseertion
      a[,paste0(sup,'bp')][a$genotype==genome & a$hapid %in% tebp$hapid]=tebp$out[match(a$hapid[a$genotype==genome& a$hapid %in% tebp$hapid], tebp$hapid)]
     # a[,paste0(sup,'bp')][a$genotype==genome][unique(subjectHits(thro))]=tebp

    }
   ### for completeness, add in centromere, knob, ribosome, telomere
   for(classif in names(otherClassifications)){
      tehapintersect=findOverlapPairs(haps, reduce(te[te$Classification %in% unlist(otherClassifications[classif]),]), ignore.strand=T)
      ov=pintersect(tehapintersect, ignore.strand=T)
     tebp=data.frame(ov) %>% group_by(hapid) %>% dplyr::summarize(out=sum(width))
      colname=paste0(classif, 'bp')
      a[,colname][a$genotype==genome & a$hapid %in% tebp$hapid]=tebp$out[match(a$hapid[a$genotype==genome& a$hapid %in% tebp$hapid], tebp$hapid)]
      a[,colname]=as.numeric(a[,colname])

    }
                  
   ### for completeness, add specific counts of knob180, tr1, and crm
   for(classif in names(specificRepeats)){
       tehapintersect=findOverlapPairs(haps, reduce(te[te$Classification %in% unlist(otherClassifications[classif]),]), ignore.strand=T)
      ov=pintersect(tehapintersect, ignore.strand=T)
     tebp=data.frame(ov) %>% group_by(hapid) %>% dplyr::summarize(out=sum(width))
      colname=paste0(classif, 'bp')
      a[,colname][a$genotype==genome & a$hapid %in% tebp$hapid]=tebp$out[match(a$hapid[a$genotype==genome& a$hapid %in% tebp$hapid], tebp$hapid)]
      a[,colname]=as.numeric(a[,colname])

    }
                  
              

print(genome)
} ## running through here overnight!!!

                
### NOTE - there will be slight discrepancies when there are overlapping TE families, between rowSums of TE superfams, and TEbp!!! TEbp is always smaller
                
a$nonTEbp=a$endmax+1-a$startmin-a$TEbp
write.table(a, paste0('allNAM_hapids.TEbpUpdate.sup.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

## when nonB73 haplotypes overlap, I can end up with negative TE bp in a region
## will do more troubleshooting, but for now, set to NA and continue
#anona=a
#anona[anona$nonTEbp<0, 11:32]=NA
#write.table(anona, paste0('allNAM_hapids.TEbpUpdate.sup.overlappingRRNA.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)


## testing
#reala=a
#a=reala[reala$refrange %in% 51085:51120,]

#pdf('~/transfer/problem_rr.pdf',10,25)
## oh7b coords are way off because of translocation, skip it
#ggplot(a[a$genotype!='Oh7B',], aes(x=startmin, y=refrange, xend=endmax, yend=refrange, color=factor(refrange))) + geom_segment() + facet_wrap(~genotype, ncol=1)
#dev.off()

