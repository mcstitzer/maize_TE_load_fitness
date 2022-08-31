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

gs=read.table('../imputation/ril_bp_repeats.noB73filter.2022-08-31.txt', header=T, comment.char='')

gs$namRIL=substr(gs$id,1,9)
gs$namFamily=substr(gs$id,1,4)

gs$nam=namfams$V2[match(gs$namFamily, namfams$V1)]
gs$subpop=nam$subpop[match(toupper(gs$nam), toupper(nam$genome))]
gs$subpop[gs$subpop=='B73']=NA

gs$keep=gs$id %in% sk$V1


## these generate data frame with single measurment per ril, from the conservative keep list from cinta
tep=merge(gs[gs$keep,], p[,-c(2:8)], by.x='namRIL', by.y='Geno_Code')
teh=merge(gs[gs$keep,], h[,-c(1)], by.x='namRIL', by.y='ID')


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
pahList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnH(x))
pahDF=do.call('rbind', pahList)

## then repeat for family-specific!
pahFList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnH(x, famSpecific=T))
pahFDF=do.call('rbind', pahFList)

write.table(pahDF, paste0('lm_output_gphenos.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(pahFDF, paste0('lm_output_gphenos.namFamily.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

## continue for merritt's phenos
pahPList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnP(x))
pahPDF=do.call('rbind', pahPList)

## then repeat for family-specific!
pahPFList=lapply(colnames(gs)[c(2:16,19:(which(colnames(gs)=='namRIL')-1))], function(x) runAssnP(x, famSpecific=T))
pahPFDF=do.call('rbind', pahPFList)
# diff in filename is mphenos not gphenos
write.table(pahPDF, paste0('lm_output_mphenos.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(pahPFDF, paste0('lm_output_mphenos.namFamily.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

                 
                 
## also write out teh and tep, to use for figures - these are the merged geno and pheno
                 
write.table(teh, paste0('geno_pheno_gphenos.',Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
write.table(tep, paste0('geno_pheno_mphenos.',Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)
