library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rtracklayer)
library(plyr)
library(stringr)

source('../figures/color_palette.R')


tem=import.gff('../genomes_and_annotations/tes/new_edta_tes/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz') ## 10/25 goes up through chr4
seqlevels(tem)=gsub('B73_', '', seqlevels(tem))

gene=import.gff('../genomes_and_annotations/genes/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3')

## oh dangit I have to filter to core, nontandem
  fgs=data.frame(fread('../genomes_and_annotations/genes/pan_gene_matrix_v3_cyverse.csv'))
  fgs.genome=fgs[,'B73'][fgs$class=='Core Gene' | fgs$class=='Near-Core Gene']
  nontandemcore=gene[gene$canonical_transcript=='1' & !is.na(gene$canonical_transcript) & gene$ID %in% fgs.genome & !is.na(gene$ID),]
  ## there are core genes without a model called, so add these
  coreregions=GRanges(seqnames=str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,1], 
                      IRanges(start=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,1]), 
                              end=as.numeric(str_split_fixed(str_split_fixed(gsub('gmap_ID=', '',fgs.genome[grepl('gmap', fgs.genome)]), ':', 2)[,2], "-", 2)[,2])))


  tem$genedist=NA
  tem$genedist[queryHits(distanceToNearest(tem, nontandemcore))]=mcols(distanceToNearest(tem, nontandemcore))$distance
  tem$coredist=NA
  tem$coredist[queryHits(distanceToNearest(tem, c(nontandemcore,coreregions)))]=mcols(distanceToNearest(tem, c(nontandemcore, coreregions)))$distance

tem$Name=gsub('_LTR', '', gsub('_INT', '', tem$Name))
tem$Identity=as.numeric(tem$Identity)
             
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


lens=data.frame(tem) %>% filter(Classification %in% unlist(classificationTE)) %>% group_by(Name, Classification) %>% 
dplyr::summarize(meangenedist=mean(genedist), mingenedist=min(genedist), maxgenedist=max(genedist),
                 meancoredist=mean(coredist), mincoredist=min(coredist), maxcoredist=max(coredist))

write.table(lens, 'te_fam_gene_dist_B73.txt', sep='\t', quote=F, row.names=F, col.names=T)
