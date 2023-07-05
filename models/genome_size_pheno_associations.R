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
mGYPHZ51=lm(GY~tebpPHZ51+nontenonrepeatbpPHZ51+allknobbpPHZ51+centromerebpPHZ51+telomerebpPHZ51+ribosomalbpPHZ51+b73bp, data=teh)
                 
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
modelsummary(list("DTS"=mDTSPHZ51, "GY"=mGYPHZ51),  output='~/transfer/te_model_table.tex', fmt=f, stars=T, statistic = 'conf.int')

                 
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


                 
                 
