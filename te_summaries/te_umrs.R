library(stringr)
library(data.table)
library(rtracklayer)
library(assertthat)
library(dplyr)



######################################################
## Read in TE annotation  ##
######################################################

### do this across all NAM
genomes=c('B73', 'CML103', 'CML247', 'CML52', 'CML69', 'Ky21', 'M162W', 'Mo18W', 'Oh43', 'P39', 
          'B97', 'CML228', 'CML277', 'CML322', 'CML333', 'HP301', 'Il14H', ## note here the capitalization of Il14H!!!!! phg uses IL14H
         'Ki11', 'Ki3', 'M37W', 'Ms71','NC350', 'NC358', 'Oh7B', 'Tx303', 'Tzi8')
              
## read in all TE annotations, sum by bp, then ask which contribute most bp to entire NAM
allTE=lapply(genomes, function(genome){
           ## read in EDTA TE annotation
 ### filenames have different capitalization than before :( so sad
 filegenome=genome
 te=read.table(paste0(genome, '-TE-UMR.bed'))
 ## the new ones go B73_chr1
 te$genome=genome
 return(te)
          })

at=do.call(rbind, allTE)
## need to remove _INT and _LTR from the end
at$Name=gsub('Name=', '', str_split_fixed(at$V9, ';', 3)[,2])
at$Name=gsub('_LTR', '', gsub('_INT', '', at$Name))
at$Classification=gsub('Classification=', '', str_split_fixed(at$V9, ';', 4)[,3])

             
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


lens=data.frame(at) %>% filter(Classification %in% unlist(classificationTE)) %>% group_by(Name, Classification) %>% 
dplyr::summarize(umrCount=sum(V10),
                 umrBp=sum(V11),
                 teBp=sum(V12),
                 umrCoverage=mean(V13))

write.table(lens, 'te_umr_NAM.txt', sep='\t', quote=F, row.names=F, col.names=T)

lensNAM=data.frame(at) %>% filter(Classification %in% unlist(classificationTE)) %>% group_by(Name, Classification, genome) %>% 
dplyr::summarize(umrCount=sum(V10),
                 umrBp=sum(V11),
                 teBp=sum(V12),
                 umrCoverage=mean(V13))



write.table(lensNAM, 'te_umr_eachNAM.txt', sep='\t', quote=F, row.names=F, col.names=T)
