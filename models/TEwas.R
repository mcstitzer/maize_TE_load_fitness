library(stringr)
library(dplyr)
library(lme4)
library(lme4qtl)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)
library(glmnet)
library(glmnetUtils)
library(RColorBrewer)

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

fams=read.table('../imputation/ril_TEFam_bp_repeats.2022-08-12.txt', header=T)

fams$namRIL=substr(fams$RIL,1,9)
fams$namFamily=substr(fams$RIL,1,4)

fams$nam=namfams$V2[match(fams$namFamily, namfams$V1)]
fams$subpop=nam$subpop[match(toupper(fams$nam), toupper(nam$genome))]
fams$subpop[fams$subpop=='B73']=NA

fams$keep=fams$RIL %in% sk$V1


## these generate data frame with single measurment per ril, from the conservative keep list from cinta
teh=merge(fams[fams$keep,], h[,-c(1)], by.x='namRIL', by.y='ID')



## from ridge regression, each fam
fteh=teh


## regular teh
teh=read.table('geno_pheno_gphenos.greaterthan22b73correct.2022-09-15.txt', header=T)

ateh=merge(fteh, teh)

afgy=lapply(colnames(fteh)[3:308], function(x) broom::tidy(lm(ateh$GY~ateh[,x] + (ateh$tebp-ateh[,x]) + ateh$nontebp)))
afdts=lapply(colnames(fteh)[3:308], function(x) broom::tidy(lm(ateh$DTS~ateh[,x] + (ateh$tebp-ateh[,x]) + ateh$nontebp)))

r2gym=lapply(colnames(fteh)[3:308], function(x) broom::glance(lm(ateh$GY~ateh[,x] + (ateh$tebp-ateh[,x]) + ateh$nontebp)))
r2dtsm=lapply(colnames(fteh)[3:308], function(x) broom::glance(lm(ateh$DTS~ateh[,x] + (ateh$tebp-ateh[,x]) + ateh$nontebp)))

estimatesgy=unlist(lapply(afgy, function(x) x[2,2]))
pvalsgy=unlist(lapply(afgy, function(x) x[2,5]))
estimatesdts=unlist(lapply(afdts, function(x) x[2,2]))
pvalsdts=unlist(lapply(afdts, function(x) x[2,5]))
r2dts=unlist(lapply(r2dtsm, function(x) x[1,1]))
r2gy=unlist(lapply(r2gym, function(x) x[1,1]))

tewas=data.frame(fam=colnames(fteh)[3:308], estimatesgy=estimatesgy, pvalsgy=pvalsgy, r2gy=r2gy, estimatesdts=estimatesdts, pvalsdts=pvalsdts, r2dts=r2dts)

## all in one model, don't double count residual variance from allTE
atteh=merge(fteh[,3:308], teh[,c('GY','tebp', 'nontebp')])
atteh$tebp=atteh$tebp-rowSums(fteh[,3:308])
alltogether=lm(GY~., data=atteh)

                       
                       
                       
mbTEs=read.table('../imputation/TE_content_across_NAM_genotypes_by_fam.2022-08-08.txt', header=T)
mbTEs$Name=gsub('_LTR', '', gsub('_INT', '', mbTEs$Name))
tewas$Classification=mbTEs$Classification[match(tewas$fam, mbTEs$Name)]
tewas$totalbp=mbTEs$bp[match(tewas$fam, mbTEs$Name)]

tewas$superfam=str_split_fixed(tewas$Classification, '/', 2)[,2]
tewas$superfam[tewas$superfam%in%c('CRM','Gypsy')]='RLG'
tewas$sup=tewas$superfam
tewas$sup[tewas$sup=='Copia']='RLC'
tewas$sup[tewas$sup=='Helitron']='DHH'
tewas$sup[tewas$sup=='L1']='RIL'
tewas$sup[tewas$sup=='unknown']='RLX'
                                                                                            
                       
## x length, dist to gene
## y neg log10 pvalue - bonforonni at 0.000162



TESUPFACTORLEVELS=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')        


tewas$sup=factor(tewas$sup, levels=TESUPFACTORLEVELS)
                       
                       
                       

len=read.table('../te_summaries/te_lengths_NAM.txt', header=T)
umr=read.table('../te_summaries/te_umr_NAM.txt', header=T)
gene=read.table('../te_summaries/te_fam_gene_dist_B73.txt', header=T)

tewas$meanlength=len$meanlength[match(tewas$fam, len$Name)]
tewas$umrCount=umr$umrCount[match(tewas$fam, umr$Name)]
tewas$umrCount[is.na(tewas$umrCount)]=0
tewas$meangenedist=gene$meangenedist[match(tewas$fam, gene$Name)]
tewas$meancoredist=gene$meancoredist[match(tewas$fam, gene$Name)]
tewas$mincoredist=gene$mincoredist[match(tewas$fam, gene$Name)]

                       
table(tewas[tewas$pvalsgy<(0.05/307) & tewas$estimatesgy>0,]$sup)
table(tewas[tewas$pvalsgy<(0.05/307) & tewas$estimatesgy<0,]$sup)
table(tewas[tewas$pvalsdts<(0.05/307) & tewas$estimatesdts>0,]$sup)
table(tewas[tewas$pvalsdts<(0.05/307) & tewas$estimatesdts<0,]$sup)
table(tewas$estimatesgy[tewas$pvalsgy<(0.05/307)]>0 )
table(tewas$estimatesdts[tewas$pvalsdts<(0.05/307)]>0 )
                      
ab=tewas
dd.col=tecolors[TESUPFACTORLEVELS]
   dtsfam=ggplot(ab[-1,], aes(x=meancoredist, y=-log10(pvalsdts), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average Distance to Gene') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 
  gyfam=ggplot(ab[-1,], aes(x=meancoredist, y=-log10(pvalsgy), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average Distance to Gene') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 

   dtsfaml=ggplot(ab[-1,], aes(x=meanlength, y=-log10(pvalsdts), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average TE Length') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 
  gyfaml=ggplot(ab[-1,], aes(x=meanlength, y=-log10(pvalsgy), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average TE Length') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 

   dtsfamlv=ggplot(ab[-1,], aes(x=estimatesdts, y=-log10(pvalsdts), color=sup, alpha=ifelse(-log10(pvalsdts)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2))
  gyfamlv=ggplot(ab[-1,], aes(x=estimatesgy, y=-log10(pvalsgy), color=sup, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) 

   r2basedts=broom::glance(lm(ateh$DTS~ateh$tebp + ateh$nontebp))[1,1]
   r2basegy=broom::glance(lm(ateh$GY~ateh$tebp + ateh$nontebp))[1,1]
dtsfamlvr2=ggplot(ab[-1,], aes(x=estimatesdts, y=r2dts, color=sup, alpha=ifelse(-log10(pvalsdts)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=unlist(r2basedts), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('R2 of model with this TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2))
  gyfamlvr2=ggplot(ab[-1,], aes(x=estimatesgy, y=r2gy, color=sup, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=unlist(r2basegy), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('R2 of model with this TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) 

                   
   dtsfamlvs=ggplot(ab[-1,], aes(x=estimatesdts, y=-log10(pvalsdts), color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsdts)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2))
  gyfamlvs=ggplot(ab[-1,], aes(x=estimatesgy, y=-log10(pvalsgy), color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) 

   dtsfamlvsl=ggplot(ab[-1,], aes(x=estimatesdts, y=meanlength, color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsdts)>-log10(0.05/307), 1, 0.5))) +geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('Average TE length')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2))
  gyfamlvsl=ggplot(ab[-1,], aes(x=estimatesgy, y=meanlength, color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5)))+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('Average TE length')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) 
   dtsfamlvsd=ggplot(ab[-1,], aes(x=estimatesdts, y=meancoredist, color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsdts)>-log10(0.05/307), 1, 0.5))) +geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('Mean Core Gene Distance')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) + scale_y_log10()
  gyfamlvsd=ggplot(ab[-1,], aes(x=estimatesgy, y=meancoredist, color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5)))+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('Mean Core Gene Distance')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) + scale_y_log10()

                       
pdf('~/transfer/cshl_tewas.pdf', 6, 2.5)
dtsfam
gyfam
dtsfaml
gyfaml
dtsfamlv
gyfamlv
dtsfamlvr2
gyfamlvr2
dtsfamlvs
gyfamlvs
dtsfamlvsl
gyfamlvsl
dtsfamlvsd
gyfamlvsd
                       
ggplot(ab[-1,], aes(x=estimatesgy, y=estimatesdts, color=sup, size=totalbp/26, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307) & -log10(pvalsdts)>-log10(0.05/307), 1, 0.5)))+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('Effect of one bp on DTS (days)')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2))

dev.off()
