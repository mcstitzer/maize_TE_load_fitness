setwd('../te_summaries/')
library(stringr)
library(data.table)
library(rtracklayer)
library(assertthat)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

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
 if(genome == 'Oh7B'){filegenome='Oh7b'}
 edtafilename=list.files(path='../genomes_and_annotations/tes/new_edta_tes/', pattern=paste0(filegenome, '.+', 'split.gff3.gz'), full.names=T)
 te=import.gff(edtafilename)
 ## the new ones go B73_chr1
 seqlevels(te)=str_split_fixed(seqlevels(te), '_', 2)[,2]
 te$genome=genome
 return(te)
          })
at=do.call(c, allTE)
## need to remove _INT and _LTR from the end
at$Name=gsub('_LTR', '', gsub('_INT', '', at$Name))
at$Identity=as.numeric(at$Identity)
             
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
head(at)
tefocusRR=read.table('../models/te_fam_ridgeregression.2022-09-13.txt', header=T, sep='\t')
head(tefocusRR)

for(chr in 1:10){
print(chr)
pdf(paste0('~/transfer/tefams_across_nam_chr', chr, '.pdf'), 20,30)
for(fam in tefocusRR$term[-1]){
print(fam)
print(ggplot(data.frame(at[seqnames(at)==paste0('chr', chr) & at$Name==tefocusRR$term[which(tefocusRR$term==fam)],]), aes(x=start, group=genome)) + geom_histogram(binwidth=1e6) + facet_wrap(~genome, ncol=1) + ggtitle(fam))
}
dev.off()
}


