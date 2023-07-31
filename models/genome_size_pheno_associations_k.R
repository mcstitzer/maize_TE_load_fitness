library(stringr)
library(dplyr)
library(lme4)
library(lme4qtl)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)

source('../figures/color_palette.R')
namfams=read.table('../figures/nam_fams.txt')

sk=read.table('../imputation/SampleToKeep.txt')


## check imputation to get tecounts object... i'm being lazy

p=read.csv('../phenotypes/all_NAM_phenos.csv')
h=read.csv('../phenotypes/BLUEs.csv')

h$ID=str_split_fixed(h$genotype, '/', 2)[,1] ## pull off the PHZ51 tester

hp=readRDS('../phenotypes/NAM_H-pheno.rds')
hp$genotype=rownames(hp)
hp$ID=str_split_fixed(rownames(hp), '/', 2)[,1] ## pull off the PHZ51 tester

## guillaume did two grain yield blues - one adjusted by flowering time (I'm calling GY), and another with the raw values (I'm calling GYraw)
h$GYraw=hp$GY[match(h$ID, hp$ID)] ## 

#gs=read.table('../imputation/ril_bp_repeats.noB73filter.2022-08-31.txt', header=T, comment.char='')

#for(b73CorrectNumber in c(0,5,10,20,22,37)){
b73CorrectNumber=22
gs=read.table(paste0('../imputation/ril_bp_repeats.greaterthan',b73CorrectNumber,'B73correct.2022-09-15.txt'), header=T, comment.char='')

gs$namRIL=substr(gs$id,1,9)
gs$namFamily=substr(gs$id,1,4)

gs$nam=namfams$V2[match(gs$namFamily, namfams$V1)]
gs$subpop=nam$subpop[match(toupper(gs$nam), toupper(nam$genome))]
gs$subpop[gs$subpop=='B73']=NA

gs$keep=gs$id %in% sk$V1


## these generate data frame with single measurment per ril, from the conservative keep list from cinta
tep=merge(gs[gs$keep,], p[,-c(2:8)], by.x='namRIL', by.y='Geno_Code')
teh=merge(gs[gs$keep,], h[,-c(1)], by.x='namRIL', by.y='ID')



kin=read.table('../phenotypes/kinship_nam_all_chroms_subset.txt', header=F, skip=3)
kinm=as.matrix(kin[,-1])
rownames(kinm)=kin$V1
colnames(kinm)=kin$V1


GYgrm=kinm[rownames(kinm)%in% teh[!is.na(teh$GY),'namRIL'], colnames(kinm) %in% teh[!is.na(teh$GY), 'namRIL']]
## weird why aren't this in kinship?
teh$namRIL[!teh$namRIL %in% colnames(GYgrm) & !is.na(teh$GY)]
#[1] "Z003E0064" "Z015E0025"


GYfit <- relmatLmer(GY ~ tebp + nontebp + b73bp + (1|namRIL) ,
                  data=teh[teh$namRIL %in% colnames(GYgrm),],
                  relmat=list(namRIL=GYgrm))
                  
VarProp(GYfit) ## prop is h2 of random effect of kinship                  
        
DTSgrm=kinm[rownames(kinm)%in% teh[!is.na(teh$DTS),'namRIL'], colnames(kinm) %in% teh[!is.na(teh$DTS), 'namRIL']]
## weird why aren't this in kinship?
teh$namRIL[!teh$namRIL %in% colnames(DTSgrm) & !is.na(teh$DTS)]
#[1] "Z003E0064" "Z015E0025"


DTSfit <- relmatLmer(DTS ~ tebp + nontebp + b73bp + (1|namRIL) ,
                  data=teh[teh$namRIL %in% colnames(DTSgrm),],
                  relmat=list(namRIL=DTSgrm))
                            
VarProp(DTSfit)



## without b73
GYfitb <- relmatLmer(GY ~ tebp + nontebp + (1|namRIL) ,
                  data=teh[teh$namRIL %in% colnames(GYgrm),],
                  relmat=list(namRIL=GYgrm))
                  
VarProp(GYfitb) ## prop is h2 of random effect of kinship  

DTSfitb <- relmatLmer(DTS ~ tebp + nontebp + (1|namRIL) ,
                  data=teh[teh$namRIL %in% colnames(DTSgrm),],
                  relmat=list(namRIL=DTSgrm))
                            
VarProp(DTSfitb)


## get variance of tebp explained by kinship!!! plus confidence interval
TEbpgrm=kinm[rownames(kinm)%in% teh[!is.na(teh$tebp),'namRIL'], colnames(kinm) %in% teh[!is.na(teh$tebp), 'namRIL']]
kinte=relmatLmer(tebp~(1|namRIL),data=teh[teh$namRIL %in% colnames(TEbpgrm),],
                 relmat=list(namRIL=TEbpgrm))
VarProp(kinte)
prof=profile(kinte)
prof_prop <- lme4qtl::varpropProf(prof) # convert to proportions
confint(prof_prop)



### first model - te vs non te
GYfit <- relmatLmer(GY ~ tebp + nontebp + (1|namRIL) ,data=teh[teh$namRIL %in% colnames(GYgrm),], relmat=list(namRIL=GYgrm))                  
VarProp(GYfit) ## prop is h2 of random effect of kinship                  

## next, repeat classes
GYfitrepeat <- relmatLmer(GY ~ tebp + allknobbp + centromerebp + telomerebp + ribosomalbp +  (1|namRIL) ,data=teh[teh$namRIL %in% colnames(GYgrm),], relmat=list(namRIL=GYgrm))                  
VarProp(GYfitrepeat) ## prop is h2 of random effect of kinship                  

## next, families 
fams=read.table('../imputation/TE_content_across_NAM_genotypes_by_collapsedfam.2023-05-16.txt', header=T)
tehfam=read.table('../imputation/ril_TEFam_bp_repeats.2022-08-12.txt', header=T)
## now have to combine columns by collapsed fam names
combofams=data.frame(lapply(fams$collapsedFam[fams$bp>10e6], function(x) rowSums(tehfam[grep(x, names(tehfam))])))
colnames(combofams)=fams$collapsedFam[fams$bp>10e6]
combofams=combofams[,colSums(combofams!=0)>0]
combofams=cbind(data.frame(namRIL=tehfam$RIL), combofams)
tehcombofams=merge(teh, combofams, by.x='id', by.y='namRIL')
GYfitfams <- relmatLmer(GY ~ . + (1|namRIL) ,data=tehcombofams[tehcombofams$namRIL %in% colnames(GYgrm),c(2,31,33:100)], relmat=list(namRIL=GYgrm))                  
GYfitfamslm=lm(GY~., data=tehcombofams[tehcombofams$namRIL %in% colnames(GYgrm),c(31,33:222)])
VarProp(GYfitfams) ## prop is h2 of random effect of kinship                  

## next, recent insertions 
tec=read.table('../imputation/ril_tecategories_bp.2023-04-10.txt', header=T)
tehcats=merge(teh, tec, by.x='id', by.y='RIL')

GYfitrecent <- relmatLmer(GY ~ recentInsertion + olderInsertion + (1|namRIL) ,data=tehcats[tehcats$namRIL %in% colnames(GYgrm),], relmat=list(namRIL=GYgrm))                  
VarProp(GYfitrecent) ## prop is h2 of random effect of kinship                  

GYfitgenic <- relmatLmer(GY ~ ingene + onekb + fivekb + greater + (1|namRIL) ,data=tehcats[tehcats$namRIL %in% colnames(GYgrm),], relmat=list(namRIL=GYgrm))                  
VarProp(GYfitgenic) ## prop is h2 of random effect of kinship                  

## this has to refit, so takes a while
anova(GYfitKin, GYfit)
