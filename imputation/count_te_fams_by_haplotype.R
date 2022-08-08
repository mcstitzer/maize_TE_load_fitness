library(stringr)
library(data.table)
library(rtracklayer)
library(assertthat)


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
b=read.table('hapids_to_refranges.txt', comment.char='', header=T)
a$refrange=b$refrange[match(a$hapid, b$hapid)]


######################################################
## Read in TE annotation to count TEs per haplotype ##
######################################################

## make a different data frame for each of the NAM genomes, assign each TE in each genome into a haplotype

### do this across all NAM
genomes=c('B73', 'CML103', 'CML247', 'CML52', 'CML69', 'Ky21', 'M162W', 'Mo18W', 'Oh43', 'P39', 
          'B97', 'CML228', 'CML277', 'CML322', 'CML333', 'HP301', 'Il14H', ## note here the capitalization of Il14H!!!!! phg uses IL14H
         'Ki11', 'Ki3', 'M37W', 'Ms71','NC350', 'NC358', 'Oh7B', 'Tx303', 'Tzi8')

                
#### add TE families!

############ okay, so get the biggest families across all genotypes
## read in all gff's so that I can find biggest by bp across all genotypes

                
## read in all TE annotations, sum by bp, then ask which contribute most bp to entire NAM
allTE=lapply(genomes, function(genome){
           ## read in EDTA TE annotation
 ### filenames have different capitalization than before :( so sad
 filegenome=genome
 if(genome == 'Oh7B'){filegenome='Oh7b'}
 edtafilename=list.files(path='../genomes_and_annotations/tes/new_edta_tes/', pattern=paste0(filegenome, '.+', 'split.gff3.gz'), full.names=T)
 te=import.gff(edtafilename)
 ## the new ones go B73_chr1
 seqlevels(te)=str_split_fixed(seqlevels(te), '_', 2)[,2]
 return(te)
          })

at=do.call(c, allTE)
## need to remove _INT and _LTR from the end
at$Name=gsub('_LTR', '', gsub('_INT', '', at$Name))
mbTEs=data.frame(at) %>% group_by(Name, Classification) %>% dplyr::summarize(bp=sum(width)) %>% arrange(desc(bp)) %>% data.frame()

write.table(mbTEs, paste0('TE_content_across_NAM_genotypes_by_fam.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
             
                
                

              
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
                ## here, any te family with >10mb across all genomes
for(i in mbTEs$Name[mbTEs$bp>10000000 & !mbTEs$Classification %in% unlist(otherClassifications)]){
  a[,paste0(i)]=0
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

 ## match to what i did to collect names, removing int and ltr from families to get match
 te$Name=gsub('_LTR', '', gsub('_INT', '', te$Name))

          
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

  ## make it easier on me and only keep TEs from the families I care about!!
  te=te[te$Name %in% mbTEs$Name[mbTEs$bp>10000000],]
          
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]


  ### by "family"
  for(fam in mbTEs$Name[mbTEs$bp>10000000]){
      tehapintersect=GenomicRanges::intersect(reduce(te[te$Name==fam & !is.na(te$Name),]), haps, ignore.strand=T)
      thro1=findOverlaps(tehapintersect, haps, ignore.strand=T) ## subset first and then count overlaps
      ov=pintersect(tehapintersect[queryHits(thro1)],haps[subjectHits(thro1)], ignore.strand=T) ## pairwise, so won't collapse the 1bp apart adjacent hits
      thro=findOverlaps(ov, haps, ignore.strand=T)
      tebp=sapply(unique(subjectHits(thro)), function(x) sum(width(ov[queryHits(thro)[subjectHits(thro)==x],]))) ## for individual sup/fam, go through reduce with each sup/fam - this is relaly far away from reality of te inseertion
  a[,fam][a$genotype==genome][unique(subjectHits(thro))]=tebp

    }      

print(genome)
} ## running through here overnight!!!

write.table(a, paste0('allNAM_hapids.FamiliesUpdate.sup.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
