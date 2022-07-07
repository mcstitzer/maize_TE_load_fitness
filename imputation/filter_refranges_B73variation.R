library(rtracklayer)
library(vcfR)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)
library(data.table)
library(dplyr)
library(stringr)

options(java.parameters = c("-Xmx32g", "-Xms32g"))

library(rPHG)

rPHG::startLogger(fullPath = NULL, fileName = NULL)

myCon <- rPHG::BrapiCon("cbsudc01.biohpc.cornell.edu")
myCon

myCon %>% rPHG::availablePHGMethods() ## need to pick a method

## use Ames genotyping to pull out B73 replicates
amesMethod <- "Ames_NonMergedReadMapping_AllLines_Haploid" ## there are several B73's - this will hopefully ask if they show differences from refgen
nam=read.table('../figures/nam_fams.txt')


## pull out all NAM parents
samples=myCon %>% 
    rPHG::PHGMethod(amesMethod) %>%
    rPHG::readSamples() %>%
    filter(str_split_fixed(toupper(sampleName), '-', 2)[,1] %in% toupper(c('B73',nam$V2[-17]))) ## don't want Mo17



## just get nam samples from Ames
ames.haps=rbind(myCon %>% 
    rPHG::PHGMethod(amesMethod) %>%
    rPHG::filterSamples(samples$sampleName[1:190]) %>% ## json error that I can't figure out, doing this way to be dumb
    rPHG::readTable(index = F) ,
    myCon %>% 
    rPHG::PHGMethod(amesMethod) %>%
    rPHG::filterSamples(samples$sampleName[191:282]) %>%
    rPHG::readTable(index = F) 
    )
    
    
## set missing haplotypes to NA
ames.haps[ames.haps==-1]=NA               

## for each refrange, ask proportion B73 allele
nsamp=nrow(ames.haps)


## we'll fill this out, adding columns for each NAM parent sample
fba=data.frame(refrange=colnames(ames.haps))

## set up which hapids belong to which reference ranges
atorr=myCon %>% rPHG::PHGMethod('mummer4_PATH') %>% rPHG::readTable()
atorrm=melt(atorr) ## Var1 is assemlby, Var2 is ref range, value is hapid
atorrm$value[atorrm$value==-1]=NA


 ## loop through each sample, check if it has been genotyped in ANY of the samples of this NAM parent
for(x in 1:nsamp){
print(x)
sample=row.names(ames.haps)[x]
inbred=str_split_fixed(sample, '-', 2)[,1]
fam.haps=ames.haps[x,] ## the haploytpes that belong to this NAM parent
fam.haps=fam.haps[fam.haps!=-1] ## missing values are messing everything up
asmhaps=atorrm$value[atorrm$Var1==paste0(inbred, '_Assembly')]

fba[,sample]=fba$refrange %in% atorrm$Var2[atorrm$Var1==paste0(inbred, '_Assembly') & atorrm$value %in% fam.haps]
## some capitalization needed, quick enough to just redo
if(inbred=='MS71'){
  fba[,sample]=fba$refrange %in% atorrm$Var2[atorrm$Var1=='Ms71_Assembly' & atorrm$value %in% fam.haps]
  }
if(inbred=='OH7B'){
  fba[,sample]=fba$refrange %in% atorrm$Var2[atorrm$Var1=='Oh7B_Assembly' & atorrm$value %in% fam.haps]
}
if(inbred=='Il14H'){
  fba[,sample]=fba$refrange %in% atorrm$Var2[atorrm$Var1=='IL14H_Assembly' & atorrm$value %in% fam.haps]
}

}


## melt to plot!
fbam=melt(fba, id.vars='refrange')

a=fread('../segregation_distortion/allNAM_hapids.TEbp.sup.2022-03-17.txt')
fbam$chrom=a$chr[match(fbam$refrange, a$refrange)]

fb=fread('fb_refranges_segregationdistortion.txt')
fbam$jointB73af=fb$jointB73af[match(fbam$refrange, fb$refrange)]

fbam$nam=str_split_fixed(fbam$variable, '-', 2)[,1]
fbam$nam[fbam$nam=='IL14H']='Il14H'

## number and proportion of support by any B73
b73correct=fbam %>% filter(nam=='B73') %>% group_by(refrange, chrom) %>% summarize(nB73correct=sum(value), pB73correct=sum(value)/n())
b73correct$jointB73af=fb$jointB73af[match(b73correct$refrange, fb$refrange)]
b73correct$start=fb$start[match(b73correct$refrange, fb$refrange)]




pdf('~/transfer/adjacentrr.ames.pdf', 28,30)
ggplot(fbam[!fbam$value,], aes(x=as.numeric(refrange), color=factor(chrom), y=variable)) + geom_point()
plot_grid(
ggplot(fbam[!fbam$value & substr(fbam$variable,1,3)=='B73',], aes(x=as.numeric(refrange), color=factor(chrom), y=variable)) + geom_point(),
    ggplot(b73correct, aes(x=as.numeric(refrange), y=nB73correct, color=factor(chrom)) )+ geom_point(),
ggplot(fbam[!fbam$value & substr(fbam$variable,1,3)=='B73',], aes(x=as.numeric(refrange), color=factor(chrom), y=jointB73af)) + geom_point(), ncol=1, align='hv')
ggplot(b73correct, aes(x=nB73correct, y=jointB73af, color=factor(chrom))) + geom_point()
dev.off()
## add centromere positions
cent=import.gff3('../annotations/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz')
cent=cent[cent$Classification%in%c('Cent/CentC', 'knob/knob180', 'knob/TR-1'),]
cent=data.frame(cent)
cent$chr=gsub('chr', '', cent$seqnames)
cent$chr[grepl('scaf', cent$chr)]=NA
cent=cent[!is.na(cent$chr),]
cent$y=0.75
cent$y[cent$Classification=='knob/knob180']=0.7
cent$y[cent$Classification=='knob/TR-1']=0.65
cent$chrom=cent$chr
                         

pdf('~/transfer/adjacentrr.b73filter.pdf', 14,10)
ggplot(b73correct,aes(x=nB73correct)) + geom_histogram(binwidth=1)
ggplot(b73correct,aes(x=nB73correct, color=factor(chrom), fill=factor(chrom))) + geom_histogram(alpha=0.5, binwidth=1)

ggplot(b73correct,aes(x=nB73correct, y=jointB73af, color=factor(chrom))) + geom_point()
ggplot(b73correct, aes(x=start, y=jointB73af, color=factor(chrom), group=chrom)) + geom_point() + ylim(0,1)+ geom_point(aes(y=y, shape=Classification, x=start), inherit.aes=F, alpha=0.1, data=cent) + facet_wrap(~factor(chrom, levels=1:10), nrow=1)
ggplot(b73correct[b73correct$nB73correct>22,], aes(x=start, y=jointB73af, color=factor(chrom), group=chrom)) + geom_point()+ ylim(0,1)+ geom_point(aes(y=y, shape=Classification, x=start), inherit.aes=F, alpha=0.1, data=cent) + facet_wrap(~factor(chrom, levels=1:10), nrow=1)
ggplot(b73correct[b73correct$nB73correct>30,], aes(x=start, y=jointB73af, color=factor(chrom), group=chrom)) + geom_point()+ ylim(0,1)+ geom_point(aes(y=y, shape=Classification, x=start), inherit.aes=F, alpha=0.1, data=cent) + facet_wrap(~factor(chrom, levels=1:10), nrow=1)

dev.off()

                  
                                   

b73correct$KEEPb73correct=b73correct$nB73correct>22

write.table(b73correct, paste0('refranges_B73correctlygenotyped.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)






