library(rtracklayer)
library(stringr)
library(dplyr)

a=import.gff3('../genomes_and_annotations/tes/5_Buckler-PHZ51_mecatErrorCorrected.contigs.RepeatMasker.gff3')
a$Class=str_split_fixed(a$Target, ' ', 3)[,1]
um=read.table('../imputation/TE_content_across_NAM_genotypes_by_fam.2022-08-08.txt', header=T)
a$Classification=um$Classification[match(a$Class, um$Name)]
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


sum(width(a)[a$Classification %in% classificationTE]) ## 1952386843
sum(width(a)[a$Classification %in% classificationKnob]) ## 22369776 ## I can't find anything cytologically about phz51 knobs, but B73 this way has 36481372
sum(width(a)[a$Classification %in% classificationCent]) ## 8822077
sum(width(a)[a$Classification %in% classificationTelo]) ## 1208217
sum(width(a)[a$Classification %in% classificationRibo]) ## 2210613

## for nontenonrepeatbp
phz51=read.table('../../merritt_anchorwave/5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta.fai') ## gets the assembly size
sum(phz51$V2)-sum(width(a)) ## 215278677


## i'm gonna be lazy and just add this in the model statement


## recentinsertion
library(data.table)
b=fread('../../merritt_anchorwave/5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta.out', header=F, skip=3, fill=T)

b$recentinsertion=F
b$recentinsertion[b$V2==0 & (b$V7-b$V6)>500]=T

sum(width(a)[b$recentinsertion]) ## 988171
sum(width(a)[!b$recentinsertion]) ## 1986133372



## distance to gene

gene=import.gff3('../genomes_and_annotations/5_Buckler-PHZ51_mecatErrorCorrected.contigs.liftoff.gff3')
fgs=data.frame(fread('../genomes_and_annotations/genes/pan_gene_matrix_v3_cyverse.csv'))
fgs.genome=fgs[,'B73'][fgs$class=='Core Gene' | fgs$class=='Near-Core Gene']
nontandemcore=gene[gene$canonical_transcript=='1' & !is.na(gene$canonical_transcript) & gene$ID %in% fgs.genome & !is.na(gene$ID),]
  ## there are core genes without a model called, so add these
  coreregions=GRanges(seqnames=str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,1], 
                      IRanges(start=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,1]), 
                              end=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,2])))

  a$coredist=NA
a$coredist[queryHits(distanceToNearest(a, c(nontandemcore),ignore.strand=T))]=mcols(distanceToNearest(a, c(nontandemcore), ignore.strand=T))$distance


a$ingene=a$coredist==0
a$onekbgene=a$coredist>0 & a$coredist<1001
a$fivekbgene=a$coredist>1000 & a$coredist<5001
a$genecat=ifelse(a$ingene, 'ingene', ifelse(a$onekbgene, 'onekb', ifelse(a$fivekbgene, 'fivekb', 'greater')))

data.frame(a) %>% group_by(genecat) %>% dplyr::summarize(sum(width))
# 1 fivekb     160430082
# 2 greater   1652369149
# 3 ingene      34420637
# 4 onekb       44668864
# 5 NA          95232811

teha$fivekbPHZ51=teha$fivekb+160430082
teha$greaterPHZ51=teha$greater+1652369149
teha$onekbPHZ51=teha$onekb+44668864
teha$ingenePHZ51=teha$ingene+34420637

## okay, save bp of each TE family
a$Name=gsub('_LTR', '', gsub('_INT', '', a$Class))
## oh there are some na's - use the overall mbTEsFam from across NAM
a$Classification[is.na(a$Classification)]=mbTEsFam$Classification[match(a$collapsedFam[is.na(a$Classification)],mbTEsFam$collapsedFam)]
## okay whatever I'm not going to be the one to solve this disaster
## drop all the TEs that are in the pan-Nam library but not in any NAM individual .....
a$collapsedFam[is.na(a$Classification)]%in%mbTEsFam$collapsedFam


a$collapsedFam=a$Name
a$collapsedFam[substr(a$Classification,1,3)=='LTR' & !is.na(a$Classification) & !grepl(':|TE|chr', a$Name)]=str_split_fixed(a$collapsedFam[substr(a$Classification,1,3)=='LTR' & !is.na(a$Classification) & !grepl(':|TE|chr', a$Name)], '_', 2)[,1]

mbTEsFamPHZ51=data.frame(a) %>% group_by(collapsedFam) %>% dplyr::summarize(bp=sum(width)) %>% dplyr::arrange(desc(bp)) %>% data.frame() ## grrr htese superfamily assignements are wrong because of stupid decisions so just drop them... :(
mbTEsFamPHZ51$Classification=a$Classification[match(mbTEsFamPHZ51$collapsedFam, a$collapsedFam)] ## this will do first one it finds? so okay??
 ## reformat to get superfamily for each TE annotation
 classificationTE=c('DNA/DTA', 'DNA/DTC', 'DNA/DTH', 'DNA/DTM', 'DNA/DTT', 'DNA/Helitron', 'LINE/L1', 'LINE/RTE', 'LINE/unknown', 'LTR/CRM', 'LTR/Copia', 'LTR/Gypsy', 'LTR/unknown', 'MITE/DTA', 'MITE/DTC', 'MITE/DTH', 'MITE/DTM', 'MITE/DTT')
mbTEsFamPHZ51$sup[!is.na(mbTEsFamPHZ51$Classification)]=str_split_fixed(mbTEsFamPHZ51$Classification[!is.na(mbTEsFamPHZ51$Classification)], '/', 2)[,2] mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='Gypsy']='RLG'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='Copia']='RLC'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='unknown' & mbTEsFamPHZ51$Classification=='LTR/unknown']='RLX'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='Helitron']='DHH'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='CRM']='RLG' ## need these and following ones because not structural
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='L1']='RIL'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='RTE']='RIT'
 mbTEsFamPHZ51$sup[mbTEsFamPHZ51$sup=='unknown' & mbTEsFamPHZ51$Classification=='LINE/unknown']='RIX'
  mbTEsFamPHZ51$sup[!mbTEsFamPHZ51$Classification %in% classificationTE]=NA

write.table(mbTEsFamPHZ51, paste0('TE_content_across_PHZ51_by_collapsedfam.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
  

