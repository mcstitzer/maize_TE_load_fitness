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
## Read in TE annotation to sum up deciles of feats ##
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

                

############ okay, so get the biggest families across all genotypes
## read in all gff's so that I can find biggest by bp across all genotypes
fgs=data.frame(fread('../genomes_and_annotations/genes/pan_gene_matrix_v3_cyverse.csv'))

                
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
          
## read in gene annotation and give each TE a distance to gene
          if(genome=='B73'){
    gene=import.gff('../genomes_and_annotations/genes/NAM_gene_gffs/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3')
          }else{
                    gene=import.gff(Sys.glob(paste0('../genomes_and_annotations/genes/NAM_gene_gffs/Zm-', genome, '-REFERENCE-NAM-1.0_*.gff3')))
                    }
   fgs.genome=fgs[,genome][fgs$class=='Core Gene' | fgs$class=='Near-Core Gene']
  nontandemcore=gene[gene$canonical_transcript=='1' & !is.na(gene$canonical_transcript) & gene$ID %in% fgs.genome & !is.na(gene$ID),]
  ## there are core genes without a model called, so add these
  coreregions=GRanges(seqnames=str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,1], 
                      IRanges(start=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,1]), 
                              end=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,2])))

            te$genedist=NA
  te$genedist[queryHits(distanceToNearest(te, nontandemcore))]=mcols(distanceToNearest(te, nontandemcore, ignore.strand=T))$distance
  te$coredist=NA
  te$coredist[queryHits(distanceToNearest(te, c(nontandemcore,coreregions)))]=mcols(distanceToNearest(te, c(nontandemcore, coreregions), ignore.strand=T))$distance

### read in umr counts for each TE - these should be same dimensions as te annotation (te)
          umr=read.table(paste0('../te_summaries/',genome, '-TE-UMR.bed'))
          umr$ID=gsub('ID=', '', str_split_fixed(umr$V9, ';', 2)[,1])
          te$umrCount=0
          te$umrCount[te$ID %in% umr$ID]=umr$V10[match(te$ID[te$ID %in% umr$ID], umr$ID)]
          te$umrBp=0
          te$umrBp[te$ID %in% umr$ID]=umr$V11[match(te$ID[te$ID %in% umr$ID], umr$ID)]
          te$umrCoverage=0
          te$umrCoverage[te$ID %in% umr$ID]=umr$V13[match(te$ID[te$ID %in% umr$ID], umr$ID)]
          
          
 ## fine just put in a genotype column??
 te$genotype=genome
 return(te)
          })

at=do.call(c, allTE)
## need to remove _INT and _LTR from the end
at$Name=gsub('_LTR', '', gsub('_INT', '', at$Name))
#mbTEs=data.frame(at) %>% group_by(Name, Classification) %>% dplyr::summarize(bp=sum(width)) %>% arrange(desc(bp)) %>% data.frame()

#write.table(mbTEs, paste0('TE_content_across_NAM_genotypes_by_fam.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
       
                
                
                
## sups too
 te=at
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
                
                
##define bins
telengthbins=quantile(width(te[te$sup %in% tesup]), probs = seq(0, 1, by = .1))
genedistbins=quantile(te[te$sup %in% tesup,]$genedist, probs = seq(0, 1, by = .1), na.rm=T) ## nas from contigs?
coredistbins=quantile(te[te$sup %in% tesup,]$coredist, probs = seq(0, 1, by = .1), na.rm=T)
agebins=quantile(as.numeric(te[te$sup %in% tesup,]$Identity), probs = seq(0, 1, by = .1), na.rm=T)
umrCountbins=quantile(as.numeric(te[te$sup %in% tesup,]$umrCount), probs = seq(0, 1, by = .1), na.rm=T)
umrBpbins=quantile(as.numeric(te[te$sup %in% tesup,]$umrBp), probs = seq(0, 1, by = .1), na.rm=T)
umrCoveragebins=quantile(as.numeric(te[te$sup %in% tesup,]$umrCoverage), probs = seq(0, 1, by = .1), na.rm=T)

 te$telengthdecile=cut(width(te), breaks=telengthbins, include.lowest=T)               
 te$genedistdecile=cut(te$genedist, breaks=genedistbins, include.lowest=T)               
te$coredistdecile=cut(te$coredist, breaks=coredistbins, include.lowest=T)               
te$agedecile=cut(as.numeric(te$Identity), breaks=agebins, include.lowest=T)               

  
                
te$ingene=te$coredist==0
te$onekbgene=te$coredist>0 & te$coredist<1001
te$fivekbgene=te$coredist>1000 & te$coredist<5001
te$genecat=ifelse(te$ingene, 'ingene', ifelse(te$onekbgene, 'onekb', ifelse(te$fivekbgene, 'fivekb', 'greater')))
te$umr=te$umrCount>0

## guessing on these column names, but find the youngest insertions
te$recentinsertion=F

tel=fread('../genomes_and_annotations/tes/NAM.EDTA2.0.0.MTEC02052020.TElib.contigsizes.txt')
tel$Name=str_split_fixed(tel$V1,"\\#", 2)[,1]
te$consLen=tel$V2[match(te$Name, tel$Name)]
te$consProp=width(te)/te$consLen


te$recentinsertion[te$Method=='homology' & te$Identity==1 & te$Classification!='Simple_repeat' & width(te)>500]=T
te$recentinsertion[te$Method=='structural' & te$sup %in% c('RLC', 'RLX', 'RLG') & te$Identity==1& width(te)>500]=T

                

                
#supTEs=data.frame(te[te$te,]) %>% group_by(sup) %>% dplyr::summarize(bp=sum(width)) %>% data.frame()
#write.table(supTEs, 'TE_content_across_NAM_genotypes_by_sup.txt', row.names=F, col.names=T, sep='\t', quote=F)


              
               
### NOW, check whether TE is within the quantile or decile of length, genedist, umr, age
### THIS IS WHAT YOU NEED TO DO!!!!!
               
                
## add numeric columns for each classification category we're going to consider
## each decile of gene dist, core dist, age, telength
## just nonzero umrs (last decile, because easier)
#for(i in 1:10){
#  a[,paste0('telength', 'decile',i)]=0
#  a[,paste0('coredist', 'decile',i)]=0
#  a[,paste0('genedist', 'decile',i)]=0
#  a[,paste0('age', 'decile',i)]=0
#  }
#a[,'umrCount0']=0
#a[,'umrCountnon0']=0
#a[,'umrBp0']=0
#a[,'umrBpnon0']=0
#a[,'umrCoverage0']=0
#a[,'umrCoveragenon0']=0


## get started!
a=data.frame(a) ## just in case it's a data table              


   
                
#######################################################3
######## loop through and populate the entry of tes
## omg this library is saving my life!
library(plyranges)         
                
alldeciles=lapply(genomes, function(genome){
   
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
 
          
          ### get count of age bins
          ### get count of gene dist bins
          ### get count of UMR bins
          ### get count of TE length
          

          
  ### by "family"

      tedec=te[te$genotype==genome,] %>% group_by(coredistdecile) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      coredist=data.frame(tedechaps) %>% group_by(hapid, coredistdecile) %>% summarize(coredistbp=sum(width))

      tedec=te[te$genotype==genome,] %>% group_by(genedistdecile) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      genedist=data.frame(tedechaps) %>% group_by(hapid, genedistdecile) %>% summarize(genedistbp=sum(width))
      tedec=te[te$genotype==genome,] %>% group_by(telengthdecile) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      telength=data.frame(tedechaps) %>% group_by(hapid, telengthdecile) %>% summarize(telengthbp=sum(width))
      tedec=te[te$genotype==genome & !is.na(te$Identity),] %>% group_by(agedecile) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      age=data.frame(tedechaps) %>% group_by(hapid, agedecile) %>% summarize(agedistbp=sum(width))

      tedec=te[te$genotype==genome,] %>% group_by(umrCountnonzero=umrCount>0) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      umrCountnonzero=data.frame(tedechaps) %>% group_by(hapid, umrCountnonzero) %>% summarize(umrCountnonzerobp=sum(width))
          
      tedec=te[te$genotype==genome,] %>% group_by(recentinsertion) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      recentinsertion=data.frame(tedechaps) %>% group_by(hapid, recentinsertion) %>% summarize(recentinsertionbp=sum(width))
    
      hapinfo=merge(coredist, genedist, all=T)
      hapinfo=merge(hapinfo, telength, all=T)
      hapinfo=merge(hapinfo, age, all=T)
      hapinfo=merge(hapinfo, umrCountnonzero, all=T) ## this is getting hard when there's multiple deciles per haplotype...
      hapinfo=merge(hapinfo, recentinsertion, all=T) ## this is getting hard when there's multiple deciles per haplotype...


return(hapinfo)
}) ## running through here overnight!!!

hapdec=do.call(rbind, alldeciles)
                
write.table(hapdec, paste0('allNAM_hapids.deciles.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
                
                
                
                
                
genedist=lapply(genomes, function(genome){
   
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
   ### by "family"

      tedec=te[te$genotype==genome,] %>% group_by(genecat) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      coredist=data.frame(tedechaps) %>% group_by(hapid, genecat) %>% summarize(genecatbp=sum(width))

  
return(coredist)
}) ## running through here overnight!!!

hapdec=do.call(rbind, genedist)
                
write.table(hapdec, paste0('coregene_distance_bins.', Sys.Date(), '.txt'), quote=F, row.names=F, col.names=T, sep='\t')
                
umr=lapply(genomes, function(genome){
   
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
   ### by "family"


      tedec=te[te$genotype==genome,] %>% group_by(umrCountnonzero=umrCount>0) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      umrCountnonzero=data.frame(tedechaps) %>% group_by(hapid, umrCountnonzero) %>% summarize(umrCountnonzerobp=sum(width))

return(umrCountnonzero)
}) ## running through here overnight!!!

hapumr=do.call(rbind, umr)
                
                
write.table(hapumr, paste0('umr_bins.', Sys.Date(), '.txt'), quote=F, row.names=F, col.names=T, sep='\t')
                
                
                
              
recent=lapply(genomes, function(genome){
   
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
   ### by "family"


      tedec=te[te$genotype==genome,] %>% group_by(recentinsertion) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      recentinsertion=data.frame(tedechaps) %>% group_by(hapid, recentinsertion) %>% summarize(recentinsertionbp=sum(width))

return(recentinsertion)
}) ## running through here overnight!!!

haprecent=do.call(rbind, recent)
                
                
write.table(haprecent, paste0('recentinsertion_bins.', Sys.Date(), '.txt'), quote=F, row.names=F, col.names=T, sep='\t')



#### get deleterious disasters
te$youngand1kb=te$recentinsertion & te$onekbgene & !is.na(te$onekbgene)
dds=lapply(genomes, function(genome){
   
## get the hapids that come from this genome
 haps=GRanges(seqnames=paste0('chr', a$chr[a$genotype==genome]), IRanges(start=a$startmin[a$genotype==genome], end=a$endmax[a$genotype==genome]))
 haps$hapid=a$hapid[a$genotype==genome]
   ### by "family"


      tedec=te[te$genotype==genome,] %>% group_by(youngand1kb) %>% reduce_ranges()
      tedechaps=join_overlap_intersect(tedec, haps)
      dd=data.frame(tedechaps) %>% group_by(hapid, youngand1kb) %>% dplyr::summarize(recentinsertionbp=sum(width))

return(dd)
}) ## running through here overnight!!!

youngclose=do.call(rbind, dds)
                
                
write.table(youngclose, paste0('youngclose_bins.', Sys.Date(), '.txt'), quote=F, row.names=F, col.names=T, sep='\t')

                
