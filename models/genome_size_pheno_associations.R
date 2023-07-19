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

## get back to mean from Guillaume's paper: subtract 16.29 from GY values
h$GYgr=h$GY-16.29

#gs=read.table('../imputation/ril_bp_repeats.noB73filter.2022-08-31.txt', header=T, comment.char='')

for(b73CorrectNumber in c(0,5,10,20,22,37)){

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

### add in values from phz51 genome
teh$genomesizePHZ51=teh$genomesize+2202400220                 
teh$tebpPHZ51=teh$tebp+1952386843
teh$allknobbpPHZ51=teh$allknobbp +22369776
teh$centromerebpPHZ51=teh$centromerebp+8822077
teh$telomerebpPHZ51=teh$telomerebp+1208217
teh$ribosomalbpPHZ51=teh$ribosomalbp+2210613        
teh$nontenonrepeatbp=teh$nontebp-teh$allknobbp-teh$centromerebp-teh$telomerebp-teh$ribosomalbp
teh$nontenonrepeatbpPHZ51=teh$nontenonrepeatbp+215278677
## run a ton of linear models, start simple with pheno~geno, no covariates - split out each repeat class
## then make another one with each NAM family separately

## function for guillaume's Hybrid phenos
runAssnH=function(genocol, famSpecific=F){ ## will return df with pheno col (GY, PH, DTS), geno, intercept, effect, r2, pval
  if(famSpecific==F){
  pah=data.frame(pheno=c('DTS', 'PH', 'GY', 'GYraw'), geno=genocol)
  pah$intercept=sapply(pah$pheno, function(x) coef(lm(teh[,x]~teh[,genocol]))[1])
  pah$gsEffect=sapply(pah$pheno, function(x) coef(lm(teh[,x]~teh[,genocol]))[2])
  pah$r2=sapply(pah$pheno, function(x) summary(lm(teh[,x]~teh[,genocol]))$r.squared)
  pah$pval=sapply(pah$pheno, function(x) summary(lm(teh[,x]~teh[,genocol]))$coefficients[2,4])
  return(pah)
  }
  if(famSpecific==T){
    tmpList=lapply(unique(teh$nam), function(nam){ ## for each nam
      pah=data.frame(pheno=c('DTS', 'PH', 'GY', 'GYraw'), geno=genocol, nam=nam, intercept=NA, gsEffect=NA, r2=NA, pval=NA)
      for(ip in 1:nrow(pah)){ ## for each pheno
        if(sum(!is.na(teh[teh$nam==nam,pah$pheno[ip]]))>0){
          mod=lm(teh[teh$nam==nam,pah$pheno[ip]]~teh[,genocol][teh$nam==nam])
          pah$intercept[ip]=sapply(pah$pheno, function(x) coef(mod)[1])
          pah$gsEffect[ip]=sapply(pah$pheno, function(x) coef(mod)[2])
          pah$r2[ip]=sapply(pah$pheno, function(x) summary(mod)$r.squared)
          pah$pval[ip]=sapply(pah$pheno, function(x) summary(mod)$coefficients[2,4])
          }
       }
       return(pah)
       }
       )
       return(do.call('rbind', tmpList))
  
}}

## function for merritt's Phenos
runAssnP=function(genocol, famSpecific=F){ ## will return df with pheno col (allmerritt), geno, intercept, effect, r2, pval
  if(famSpecific==F){
   pah=data.frame(pheno=colnames(tep)[(ncol(gs)+1):ncol(tep)], geno=genocol)
   pah$intercept=sapply(pah$pheno, function(x) coef(lm(tep[,x]~tep[,genocol]))[1])
   pah$gsEffect=sapply(pah$pheno, function(x) coef(lm(tep[,x]~tep[,genocol]))[2])
   pah$r2=sapply(pah$pheno, function(x) summary(lm(tep[,x]~tep[,genocol]))$r.squared)
   pah$pval=sapply(pah$pheno, function(x) summary(lm(tep[,x]~tep[,genocol]))$coefficients[2,4])
   return(pah)
   }
  if(famSpecific==T){
      tmpList=lapply(unique(tep$nam), function(nam){ ## for each nam
      pah=data.frame(pheno=colnames(tep)[(ncol(gs)+1):ncol(tep)], geno=genocol, nam=nam, intercept=NA, gsEffect=NA, r2=NA, pval=NA)
      for(ip in 1:nrow(pah)){ ## for each pheno
        if(sum(!is.na(tep[tep$nam==nam,pah$pheno[ip]]))>0){
          mod=lm(tep[tep$nam==nam,pah$pheno[ip]]~tep[tep$nam==nam,genocol])
          pah$intercept[ip]=sapply(pah$pheno, function(x) coef(mod)[1])
          pah$gsEffect[ip]=sapply(pah$pheno, function(x) coef(mod)[2])
          pah$r2[ip]=sapply(pah$pheno, function(x) summary(mod)$r.squared)
          pah$pval[ip]=sapply(pah$pheno, function(x) summary(mod)$coefficients[2,4])
          }
       }
       return(pah)
       }
       )
       return(do.call('rbind', tmpList))

  }
}

## make a df for all geno x pheno combos!
pahList=lapply(colnames(teh)[c(3:17,20:23, 33:40)], function(x) runAssnH(x))
pahDF=do.call('rbind', pahList)

## then repeat for family-specific!
pahFList=lapply(colnames(teh)[c(3:17,20:23, 33:40)], function(x) runAssnH(x, famSpecific=T))
pahFDF=do.call('rbind', pahFList)

write.table(pahDF, paste0('lm_output_gphenos.greaterthan',b73CorrectNumber, 'b73correct.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(pahFDF, paste0('lm_output_gphenos.namFamily.greaterthan',b73CorrectNumber, 'b73correct.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

## continue for merritt's phenos
pahPList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnP(x))
pahPDF=do.call('rbind', pahPList)

## then repeat for family-specific!
pahPFList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnP(x, famSpecific=T))
pahPFDF=do.call('rbind', pahPFList)
# diff in filename is mphenos not gphenos
write.table(pahPDF, paste0('lm_output_mphenos.greaterthan',b73CorrectNumber, 'b73correct.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(pahPFDF, paste0('lm_output_mphenos.namFamily.greaterthan',b73CorrectNumber, 'b73correct.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

                 
                 
## also write out teh and tep, to use for figures - these are the merged geno and pheno
                 
write.table(teh, paste0('geno_pheno_gphenos.greaterthan',b73CorrectNumber, 'b73correct.',Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(tep, paste0('geno_pheno_mphenos.greaterthan',b73CorrectNumber, 'b73correct.',Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
}


                 
mDTS=lm(DTS~tebp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp, data=teh)
mGY=lm(GY~tebp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp, data=teh)
                 
summary(mDTS)
summary(mGY)
                 
               
mDTSPHZ51=lm(DTS~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51+b73bp, data=teh)
mGYPHZ51=lm(GYgr~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51+b73bp, data=teh)
                 
summary(mDTSPHZ51)
summary(mGYPHZ51)

                 
                 
mDTSn=lmer(DTS~tebp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp + (1|nam), data=teh)
mGYn=lmer(GY~tebp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp + (1|nam), data=teh)
                 
summary(mDTSn)
summary(mGYn)
                 
                 
mDTSnPHZ51=lmer(DTS~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51+b73bp + (1|nam), data=teh)
mGYnPHZ51=lmer(GY~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51+b73bp + (1|nam), data=teh)
                 
summary(mDTSnPHZ51)
summary(mGYnPHZ51)

### MAKE A TABLE FOR THE PAPER!!!!!
library(modelsummary)
modelsummary(list("DTS"=mDTSPHZ51, "GY"=mGYPHZ51),  output='~/transfer/te_model_table.tex', fmt=f, stars=T, statistic = NULL)


## nam family specific
library(broom)
    tmpList=lapply(unique(teh$nam), function(nam){ ## for each nam
      phenos=c('DTS', 'PH', 'GY', 'GYraw', 'GYgr')
      pahlist=lapply(phenos, function(pheno){ ## for each pheno
        if(sum(!is.na(teh[teh$nam==nam,pheno]))>0){
          tempdata=teh[teh$nam==nam,]
          tempdata$pheno=tempdata[,pheno]
          model=lm(pheno~tebpPHZ51 + allknobbpPHZ51 + centromerebpPHZ51 + telomerebpPHZ51 + ribosomalbpPHZ51 + nontenonrepeatbpPHZ51, data=tempdata)                
    
          pah=tidy(model) 
          pah$nam=nam
          pah$pheno=pheno
          pah$modelp=glance(model)$p.value 
          pah$modelR2=glance(model)$r.squared          
        }
        return(pah)
       })
       return(do.call('rbind', pahlist))
       })
       


pahFDF=do.call('rbind', tmpList)
write.table(pahFDF, paste0('lm_output_gphenos.namFamily.greaterthan',b73CorrectNumber, 'b73correct.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)


## superfamily splilt for supplement

mDTSsup=lm(DTS~dhhbp+dtabp+dtcbp+dthbp+dtmbp+dttbp+rilbp+ritbp+rlcbp+rlgbp+rlxbp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp, data=teh)
mGYsup=lm(GYgr~dhhbp+dtabp+dtcbp+dthbp+dtmbp+dttbp+rilbp+ritbp+rlcbp+rlgbp+rlxbp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp, data=teh)
modelsummary(list("DTS"=mDTSsup, "GY"=mGYsup),  output='~/transfer/supsupp_te_model_table.tex', fmt=f, stars=T, statistic = NULL)


dtslist=lapply(unique(teh$nam), function(nam){ ## for each nam
      pheno='DTS'
        if(sum(!is.na(teh[teh$nam==nam,pheno]))>0){
          tempdata=teh[teh$nam==nam,]
          tempdata$pheno=tempdata[,pheno]
          model=lm(pheno~tebpPHZ51 + allknobbpPHZ51 + centromerebpPHZ51 + telomerebpPHZ51 + ribosomalbpPHZ51 + nontenonrepeatbpPHZ51, data=tempdata)                
        }
        return(model)
       })
##dput(paste0('DTS_', dput(names(table(pahFDF$nam)))))
names(dtslist)=c("DTS_B97", "DTS_CML103", "DTS_CML228", "DTS_CML247", "DTS_CML277", 
"DTS_CML322", "DTS_CML333", "DTS_CML52", "DTS_CML69", "DTS_Il14H", 
"DTS_Ki11", "DTS_Ki3", "DTS_Ky21", "DTS_M162W", "DTS_M37W", "DTS_Mo18W", 
"DTS_Ms71", "DTS_NC350", "DTS_NC358", "DTS_Oh43", "DTS_Oh7B", 
"DTS_P39", "DTS_Tx303", "DTS_Tzi8")


gylist=lapply(unique(teh$nam), function(nam){ ## for each nam
      pheno='GYgr'
        if(sum(!is.na(teh[teh$nam==nam,pheno]))>0){
          tempdata=teh[teh$nam==nam,]
          tempdata$pheno=tempdata[,pheno]
          model=lm(pheno~tebpPHZ51 + allknobbpPHZ51 + centromerebpPHZ51 + telomerebpPHZ51 + ribosomalbpPHZ51 + nontenonrepeatbpPHZ51, data=tempdata)                
        }
        return(model)
       })
##dput(paste0('GY_', dput(names(table(pahFDF$nam)))))
names(gylist)=c("GY_B97", "GY_CML103", "GY_CML228", "GY_CML247", "GY_CML277", 
"GY_CML322", "GY_CML333", "GY_CML52", "GY_CML69", "GY_Il14H", 
"GY_Ki11", "GY_Ki3", "GY_Ky21", "GY_M162W", "GY_M37W", "GY_Mo18W", 
"GY_Ms71", "GY_NC350", "GY_NC358", "GY_Oh43", "GY_Oh7B", "GY_P39", 
"GY_Tx303", "GY_Tzi8")

                 
## nam family table for supplement
modelsummary(dtslist,  
             output='~/transfer/supsupp_namfam_dts_te_model_table.tex', fmt=f, stars=T, statistic = NULL)
modelsummary(gylist,  
             output='~/transfer/supsupp_namfam_gy_te_model_table.tex', fmt=f, stars=T, statistic = NULL)



##### gene distance bins
aa=read.table('../imputation/ril_tecategories_bp.2023-06-30.txt', header=T)
teha=merge(teh, aa, by.x='id', by.y='RIL')

## these are from imputation/get_phz51_genome_counts.R
teha$fivekbPHZ51=teha$fivekb+160430082
teha$greaterPHZ51=teha$greater+1652369149
teha$onekbPHZ51=teha$onekb+44668864
teha$ingenePHZ51=teha$ingene+34420637

gygene=lm(teha$GYgr~teha$ingenePHZ51+teha$onekbPHZ51+teha$fivekbPHZ51+teha$greaterPHZ51+teha$b73bp)
dtsgene=lm(teha$DTS~teha$ingenePHZ51+teha$onekbPHZ51+teha$fivekbPHZ51+teha$greaterPHZ51+teha$b73bp)
modelsummary(list('GY'=gygene, 'DTS'=dtsgene  ),
             output='~/transfer/te_distance_model_table.tex', fmt=f, stars=T, statistic = NULL)
           
### recent insertions
## these are from imputation/get_phz51_genome_counts.R
teha$recentInsertionPHZ51=teha$recentInsertion+988171
teha$olderInsertionPHZ51=teha$olderInsertion+1986133372

gyage=lm(teha$GYgr~teha$recentInsertionPHZ51+teha$olderInsertionPHZ51+teha$b73bp)
dtsage=lm(teha$DTS~teha$recentInsertionPHZ51+teha$olderInsertionPHZ51+teha$b73bp)
modelsummary(list('GY'=gyage, 'DTS'=dtsage  ),
             output='~/transfer/te_age_model_table.tex', fmt=f, stars=T, statistic = NULL)


## methylation - i'm being lazy and not calling phz51 
gyumr=lm(teha$GYgr~teha$noUMR+teha$hasUMR+teha$b73bp)
dtsumr=lm(teha$DTS~teha$noUMR+teha$hasUMR+teha$b73bp)
modelsummary(list('GY'=gyumr, 'DTS'=dtsumr  ),
             output='~/transfer/te_umr_model_table.tex', fmt=f, stars=T, statistic = NULL)

### model the deleterious disasters
table(te[te$recentinsertion & te$onekbgene & !is.na(te$onekbgene),]$genotype)
ddis=read.table('ril_tecategories_bp.2023-07-14.txt', header=T)
teha=merge(teha, ddis[,c('RIL', 'farorold', 'youngandclose')], by.x='id', by.y='RIL')
gydd=lm(teha$GYgr~teha$farorold+teha$youngandclose+teha$b73bp)
dtsdd=lm(teha$DTS~teha$farorold+teha$youngandclose+teha$b73bp)
modelsummary(list('GY'=gydd, 'DTS'=dtsdd  ),
             output='~/transfer/te_deldis_model_table.tex', fmt=f, stars=T, statistic = NULL)

## te fams
tefams=read.table('ril_tefamilycounts_1MbFams_bp.txt', header=T)
tefamsp=merge(tefams, teha, by.x='RIL', by.y='id')
mbTEsFamPHZ51=read.table('TE_content_across_PHZ51_by_collapsedfam.2023-07-14.txt', header=T)
tefamsp[,colnames(tefamsp)%in% mbTEsFamPHZ51$collapsedFam]=lapply(colnames(tefamsp)[colnames(tefamsp)%in% mbTEsFamPHZ51$collapsedFam], function(fam){ tefamsp[,fam]+mbTEsFamPHZ51$bp[mbTEsFamPHZ51$collapsedFam==fam]})
tefamsp$smallerFamTEbp=tefamsp$tebpPHZ51-rowSums(tefamsp[,c(mbTEsFam$collapsedFam[mbTEsFam$bp>10e6 & substr(mbTEsFam$collapsedFam,1,2)!='TE' & !mbTEsFam$collapsedFam %in% c('knob180', 'TR-1', 'CentC', 'AF013103.1', 'CL569186.1')])])
tfm=lm(GYgr~., data=tefamsp[,c('GYgr','b73bp','smallerFamTEbp', mbTEsFam$collapsedFam[mbTEsFam$bp>10e6 & substr(mbTEsFam$collapsedFam,1,2)!='TE' & !mbTEsFam$collapsedFam %in% c('knob180', 'TR-1', 'CentC', 'AF013103.1', 'CL569186.1')])])
tfmdts=lm(DTS~., data=tefamsp[,c('DTS','b73bp','smallerFamTEbp', mbTEsFam$collapsedFam[mbTEsFam$bp>10e6 & substr(mbTEsFam$collapsedFam,1,2)!='TE'& !mbTEsFam$collapsedFam %in% c('knob180', 'TR-1', 'CentC', 'AF013103.1', 'CL569186.1')])])

summary(tfm)
tfm1m=lm(GYgr~., data=tefamsp[,c('GYgr','b73bp',mbTEsFam$collapsedFam[mbTEsFam$bp>1e6 & substr(mbTEsFam$collapsedFam,1,2)!='TE'& !mbTEsFam$collapsedFam %in% c('knob180', 'TR-1', 'CentC', 'AF013103.1', 'CL569186.1')])])

 modelsummary(list('GY'=tfm, 'DTS'=tfmdts  ),
             output='~/transfer/te_fam_10m_model_table.tex', fmt=f, stars=T, statistic = NULL)


source('../figures/color_palette.R')
library(ggridges)  
library(broom)
tfmt=tidy(tfm)
tfmt$sup=mbTEsFamPHZ51$sup[match(tfmt$term, mbTEsFamPHZ51$collapsedFam)]
tfmt$sup[tfmt$sup=='Gypsy']='RLG'
tfmt$sup[tfmt$sup %in% c('RIL', 'RIT', 'RIX')]='nonLTR'
dd.col=c(dd.col, RIX='#6c4da4', RIL='#6c4da4',nonLTR='#6c4da4', nonsig='gray')



gygt=tidy(gygene)
gygt$ylab=c('intercept', 'Load within gene', 'Load 1 kb from gene', 'Load 5 kb from gene', 'Load >5 kb from gene', 'B73bp')
gyat=tidy(gyage)
gyat$ylab=c('intercept', 'Load of recent insertions', 'Load of older insertions', 'B73bp')
gyut=tidy(gyumr)
gyut$ylab=c('intercept', 'Load of no UMR', 'Load with UMR', 'B73bp')

gylim=max(abs(c(gygt$estimate[-c(1,6)], gyat$estimate[-c(1,4)], gyut$estimate[-c(1,4)])))
                 
pdf('~/transfer/fig4_tefams.pdf', 10,5)
ggplot(tfmt[-c(1:3),], aes(x=estimate, y=sup, fill=sup, color=ifelse(p.value<0.05, sup, 'nonsig'), size=ifelse(p.value<0.05,4,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_jitter(height=0.3, alpha=0.9)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(tfmt$estimate[-1])), max(abs(tfmt$estimate[-1]))) + theme(legend.position='NULL')+ ylab('TE superfamily') + xlab('Effect of one bp on GY (t/ha)')+scale_y_discrete(limits=rev)
ggplot(tfmt[-c(1:3),], aes(x=estimate, y=sup, fill=sup, color=sup, size=ifelse(p.value<0.05,4,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_jitter(height=0.3, alpha=0.9)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(tfmt$estimate[-1])), max(abs(tfmt$estimate[-1]))) + theme(legend.position='NULL')+ ylab('TE superfamily') + xlab('Effect of one bp on GY (t/ha)')+scale_y_discrete(limits=rev)
tfp=ggplot(tfmt[-c(1:3),], aes(x=estimate, y=sup, fill=sup, color=sup, size=ifelse(p.value<0.05,2,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_jitter(height=0.2, alpha=0.9)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(tfmt$estimate[-1])), max(abs(tfmt$estimate[-1]))) + theme(legend.position='NULL')+ ylab('TE superfamily') + xlab('Effect of one bp on GY (t/ha)')+scale_y_discrete(limits=rev)


gdp=ggplot(gygt[-c(1,6),], aes(x=estimate, y=ylab, color=ifelse(p.value<0.05,'black','gray'), size=ifelse(p.value<0.05,2,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_point( alpha=0.9)  + xlim(-max(gylim), max(gylim)) + theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)') + scale_color_manual(values=c(black='black', gray='gray')) + scale_y_discrete(limits = c( 'Load >5 kb from gene', 'Load 5 kb from gene',  'Load 1 kb from gene', 'Load within gene'))
gap=ggplot(gyat[-c(1,4),], aes(x=estimate, y=ylab, color=ifelse(p.value<0.05,'black','gray'), size=ifelse(p.value<0.05,2,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_point(alpha=0.9)  + xlim(-max(gylim), max(gylim)) + theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+scale_y_discrete(limits=rev) + scale_color_manual(values=c(black='black', gray='gray'))
gup=ggplot(gyut[-c(1,4),], aes(x=estimate, y=ylab, color=ifelse(p.value<0.05,'black','gray'), size=ifelse(p.value<0.05,2,1))) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_point(alpha=0.9)  + xlim(-max(gylim), max(gylim)) + theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+scale_y_discrete(limits=rev) + scale_color_manual(values=c(black='black', gray='gray'))

threem=plot_grid(gdp, gap, gup, labels='AUTO', ncol=1, align='hv')                 
plot_grid(threem, tfp, labels=c('', 'D'))
dev.off()


                 
                 
# library(rstanarm)

# model <- stan_glm(DTS~tebp+nontenonrepeatbp+allknobbp+centromerebp+telomerebp+ribosomalbp+b73bp, data=teh, refresh = 0)
# sjstats::r2(model)
# model <- stan_lmer(GY~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51 + (1|nam), data=teh, refresh = 0)

# sjstats::r2(model)

                 
# pdf('~/transfer/full_bayes.pdf', 14,8)
# result=estimate_density(model)

# plot(result)
# plot(result, stack=F)
# plot(result, stack=F, priors=T)
# result=describe_posterior(model)
# plot(result)

# result=p_direction(model)
# plot(result)                 
# dev.off()
                 
                 
# library(texreg) ## stargazer and this can't do scinot, so thin kthis through later
# mDTS=lm(DTS~tebp+allknobbp+centromerebp+telomerebp+ribosomalbp, data=teh)
# mGY=lm(GY~tebp+allknobbp+centromerebp+telomerebp+ribosomalbp, data=teh)
                 
# texreg(list(mDTS, mGY), digits=10)
                 
                 
               
                 
                 
# ## partial r2
                 
# library(rsq)
# summary(mGY)
# rsq.partial(mGY)        

                 
                 
                 
# tebp=lm(teh$GY~teh$tebp + teh$nontebp + teh$b73bp)
# genic=lm(teh$GY~teh$ingene + teh$onekb + teh$fivekb + teh$greater + teh$nontebp + teh$b73bp)
# umr=lm(teh$GY~teh$hasUMR + teh$noUMR + teh$nontebp + teh$b73bp)
# recent=lm(teh$GY~teh$recentInsertion + teh$olderInsertion + teh$nontebp + teh$b73bp)

# gs=lm(GY~genomesize+b73bp, data=teh)


# AIC(gs,tebp,genic,umr,recent) %>% arrange(AIC)
# sapply(list(gs,tebp,genic,umr,recent), function(x) summary(x)$adj.r.squared)


# models=list(gs, tebp, genic, umr, recent) 

# models=lapply(models, tidy)
# models[[1]]$name='gs'
# models[[2]]$name='tebp'
# models[[3]]$name='genic'
# models[[4]]$name='umr'
# models[[5]]$name='recent'


                 
                 
