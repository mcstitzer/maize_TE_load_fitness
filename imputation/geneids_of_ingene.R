
youngin=te[te$agegenecat=='young_ingene',]
youngin$geneid=NA


############ okay, so get the biggest families across all genotypes
## read in all gff's so that I can find biggest by bp across all genotypes
fgs=data.frame(fread('../genomes_and_annotations/genes/pan_gene_matrix_v3_cyverse.csv'))



for(genome in genomes){
  younginrows=youngin$genotype==genome
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

#            te$genedist=NA
#  te$genedist[queryHits(distanceToNearest(te, nontandemcore))]=mcols(distanceToNearest(te, nontandemcore, ignore.strand=T))$distance
#  te$coredist=NA
#  te$coredist[queryHits(distanceToNearest(te, c(nontandemcore,coreregions)))]=mcols(distanceToNearest(te, c(nontandemcore, coreregions), ignore.strand=T))$distance
youngin$geneid[younginrows][queryHits(distanceToNearest(youngin[younginrows,], nontandemcore))]=nontandemcore$ID[subjectHits(distanceToNearest(youngin[younginrows,], nontandemcore))]

}
  
