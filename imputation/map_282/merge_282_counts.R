library(data.table)
library(stringr)
library(dplyr)
a=read.csv('../phenotypes/Ames/All_Crosses_list.csv')
b=read.csv('../phenotypes/Ames/BLUEs.csv')
fq=read.table('../../282_load/282_fqnames.txt', header=F)
fq$name=str_split_fixed(fq$V1,'_',3)[,1]
a$name=str_split_fixed(a$Hybrid_pheno_name, paste0('x', a$Tester), 2)[,1]

sum(a$name %in% fq$name) ## still only 235 not 282, but a start!
sum(toupper(unique(a$name)) %in% toupper(fq$name)) ## up to 244

## 283 individuals
length(unique(fq$name))
sum(toupper(unique(fq$name)) %in% toupper(a$name)) ## 138 oh duh because there are two testers...

## what am i missing??
unique(fq$name)[!toupper(unique(fq$name)) %in% toupper(a$name)]




m=fread('multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt', header=T)
al=vector(mode='list', length=nrow(fq))
names(al)=fq$V1

for(x in 1:nrow(fq)){
  i=fread(paste0('idxstat/', fq$V1[x], '.NAMTElib.mapped.txt.idxstat.txt'))
  i$prop=i$V3/(2*m$`Unique Reads`[m$Sample==paste0(fq$V1[x], '_1')])
  i$geno=fq$name[x]
  i$sample=fq$V1[x]
  al[[x]]=i
          }

rp=do.call('rbind', al)

#a=fread('idxstat/Yu796_HLC27CCXX_L5.NAMTElib.mapped.txt.idxstat.txt')
#a$prop=a$V3/(2* m$`Unique Reads`[m$Sample=='Yu796_HLC27CCXX_L5_1'])


## milt is negatively assoociated with yield
 milt=rp[grepl('milt', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
milt$GBS=a$GBS[match(milt$geno, a$name)]
b$milt=milt$prop[match(str_split_fixed(b$genotype,'/',2)[,1], milt$GBS)]
b$tester=ifelse(grepl('B47$', b$genotype), 'B47', 'PHZ51')

## huck is negative
huck=rp[grepl('huck', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
huck$GBS=a$GBS[match(huck$geno, a$name)]
b$huck=huck$prop[match(str_split_fixed(b$genotype,'/',2)[,1], huck$GBS)]


## doke is positive
doke=rp[grepl('doke', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
doke$GBS=a$GBS[match(doke$geno, a$name)]
b$doke=doke$prop[match(str_split_fixed(b$genotype,'/',2)[,1], doke$GBS)]

## Zm_00026 is positive
zm00026=rp[grepl('DTC_ZM00026', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
zm00026$GBS=a$GBS[match(zm00026$geno, a$name)]
b$zm00026=zm00026$prop[match(str_split_fixed(b$genotype,'/',2)[,1], zm00026$GBS)]

summary(lm(b$GY~b$zm00026 + b$tester))

## mada is positive
mada=rp[grepl('mada', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
mada$GBS=a$GBS[match(mada$geno, a$name)]
b$mada=mada$prop[match(str_split_fixed(b$genotype,'/',2)[,1], mada$GBS)]

# machiavelli is positive
machiavelli=rp[grepl('machiavelli', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
machiavelli$GBS=a$GBS[match(machiavelli$geno, a$name)]
b$machiavelli=machiavelli$prop[match(str_split_fixed(b$genotype,'/',2)[,1], machiavelli$GBS)]

# sagyfy is positive
sagyfy=rp[grepl('sagyfy', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
sagyfy$GBS=a$GBS[match(sagyfy$geno, a$name)]
b$sagyfy=sagyfy$prop[match(str_split_fixed(b$genotype,'/',2)[,1], sagyfy$GBS)]


# CRM1 is negative
CRM1=rp[grepl('CRM1', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
CRM1$GBS=a$GBS[match(CRM1$geno, a$name)]
b$CRM1=CRM1$prop[match(str_split_fixed(b$genotype,'/',2)[,1], CRM1$GBS)]


# iwik is most positive
iwik=rp[grepl('iwik', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
iwik$GBS=a$GBS[match(iwik$geno, a$name)]
b$iwik=iwik$prop[match(str_split_fixed(b$genotype,'/',2)[,1], iwik$GBS)]


# eninu is negative
eninu=rp[grepl('eninu', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
eninu$GBS=a$GBS[match(eninu$geno, a$name)]
b$eninu=eninu$prop[match(str_split_fixed(b$genotype,'/',2)[,1], eninu$GBS)]

# ewib is positive
ewib=rp[grepl('ewib', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
ewib$GBS=a$GBS[match(ewib$geno, a$name)]
b$ewib=ewib$prop[match(str_split_fixed(b$genotype,'/',2)[,1], ewib$GBS)]


# DTC_ZM00081_consensus is negative
DTC_ZM00081_consensus=rp[grepl('DTC_ZM00081_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTC_ZM00081_consensus$GBS=a$GBS[match(DTC_ZM00081_consensus$geno, a$name)]
b$DTC_ZM00081_consensus=DTC_ZM00081_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTC_ZM00081_consensus$GBS)]


# DTC_ZM00004_consensus is negative
DTC_ZM00004_consensus=rp[grepl('DTC_ZM00004_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTC_ZM00004_consensus$GBS=a$GBS[match(DTC_ZM00004_consensus$geno, a$name)]
b$DTC_ZM00004_consensus=DTC_ZM00004_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTC_ZM00004_consensus$GBS)]

# DTC_ZM00012_consensus is negative
DTC_ZM00012_consensus=rp[grepl('DTC_ZM00012_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTC_ZM00012_consensus$GBS=a$GBS[match(DTC_ZM00012_consensus$geno, a$name)]
b$DTC_ZM00012_consensus=DTC_ZM00012_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTC_ZM00012_consensus$GBS)]

# DTC_ZM00102_consensus is negative
DTC_ZM00102_consensus=rp[grepl('DTC_ZM00102_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTC_ZM00102_consensus$GBS=a$GBS[match(DTC_ZM00102_consensus$geno, a$name)]
b$DTC_ZM00102_consensus=DTC_ZM00102_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTC_ZM00102_consensus$GBS)]

# DTH_ZM00024_consensus is negative
DTH_ZM00024_consensus=rp[grepl('DTH_ZM00024_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTH_ZM00024_consensus$GBS=a$GBS[match(DTH_ZM00024_consensus$geno, a$name)]
b$DTH_ZM00024_consensus=DTH_ZM00024_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTH_ZM00024_consensus$GBS)]


# DTC_ZM00063_consensus is negative
DTC_ZM00063_consensus=rp[grepl('DTC_ZM00063_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTC_ZM00063_consensus$GBS=a$GBS[match(DTC_ZM00063_consensus$geno, a$name)]
b$DTC_ZM00063_consensus=DTC_ZM00063_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTC_ZM00063_consensus$GBS)]

# DTH_ZM00032_consensus is negative
DTH_ZM00032_consensus=rp[grepl('DTH_ZM00032_consensus', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
DTH_ZM00032_consensus$GBS=a$GBS[match(DTH_ZM00032_consensus$geno, a$name)]
b$DTH_ZM00032_consensus=DTH_ZM00032_consensus$prop[match(str_split_fixed(b$genotype,'/',2)[,1], DTH_ZM00032_consensus$GBS)]





### oTHER TEs

knob=rp[grepl('knob', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
knob$GBS=a$GBS[match(knob$geno, a$name)]
b$knob=knob$prop[match(str_split_fixed(b$genotype,'/',2)[,1], knob$GBS)]

telomere=rp[grepl('telomere', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
telomere$GBS=a$GBS[match(telomere$geno, a$name)]
b$telomere=telomere$prop[match(str_split_fixed(b$genotype,'/',2)[,1], telomere$GBS)]

Cent=rp[grepl('Cent', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
Cent$GBS=a$GBS[match(Cent$geno, a$name)]
b$Cent=Cent$prop[match(str_split_fixed(b$genotype,'/',2)[,1], Cent$GBS)]

rDNA=rp[grepl('rDNA', rp$V1),] %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
rDNA$GBS=a$GBS[match(rDNA$geno, a$name)]
b$rDNA=rDNA$prop[match(str_split_fixed(b$genotype,'/',2)[,1], rDNA$GBS)]

allte=rp %>% group_by(sample, geno) %>% dplyr::summarize(prop=sum(prop)) %>% group_by(geno) %>% summarize(prop=mean(prop))
allte$GBS=a$GBS[match(allte$geno, a$name)]
b$allte=allte$prop[match(str_split_fixed(b$genotype,'/',2)[,1], allte$GBS)]


summary(lm(GY~.,data=b[,-c(1:3)]))

