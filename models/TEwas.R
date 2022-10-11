
## from ridge regression
fteh=teh
## regular teh
ateh=merge(fteh, teh)

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


afgy=lapply(colnames(fteh)[3:308], function(x) broom::tidy(lm(ateh$GY~ateh[,x] + ateh$tebp + ateh$nontebp)))
afdts=lapply(colnames(fteh)[3:308], function(x) broom::tidy(lm(ateh$DTS~ateh[,x] + ateh$tebp + ateh$nontebp)))

estimatesgy=unlist(lapply(afgy, function(x) x[2,2]))
pvalsgy=unlist(lapply(afgy, function(x) x[2,5]))
estimatesdts=unlist(lapply(afdts, function(x) x[2,2]))
pvalsdts=unlist(lapply(afdts, function(x) x[2,5]))

tewas=data.frame(fam=colnames(fteh)[3:308], estimatesgy=estimatesgy, pvalsgy=pvalsgy, estimatesdts=estimatesdts, pvalsdts=pvalsdts)


## x length, dist to gene
## y neg log10 pvalue - bonforonni at 0.000162



TESUPFACTORLEVELS=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')        

ab=merge(tefocusRR, tewas, by.x='term', by.y='fam')

ab$sup=factor(ab$sup, levels=TESUPFACTORLEVELS)
dd.col=tecolors[TESUPFACTORLEVELS]
   dtsfam=ggplot(ab[-1,], aes(x=meancoredist, y=-log10(pvalsdts), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average Distance to Gene') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 
  gyfam=ggplot(ab[-1,], aes(x=meancoredist, y=-log10(pvalsgy), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average Distance to Gene') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 

   dtsfaml=ggplot(ab[-1,], aes(x=meanlength, y=-log10(pvalsdts), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average TE Length') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 
  gyfaml=ggplot(ab[-1,], aes(x=meanlength, y=-log10(pvalsgy), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Average TE Length') + scale_x_log10() + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 

   dtsfamlv=ggplot(ab[-1,], aes(x=estimatesdts, y=-log10(pvalsdts), color=sup)) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on DTS') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') 
  gyfamlv=ggplot(ab[-1,], aes(x=estimatesgy, y=-log10(pvalsgy), color=sup, alpha=ifelse(-log10(pvalsgy)>-log10(0.05/307), 1, 0.5))) +geom_hline(yintercept=-log10(0.05/307), color='gray', lty='dashed')+geom_vline(xintercept=0, color='gray', lty='dashed') + scale_color_manual(values=dd.col) + geom_point()+ xlab('Effect of one bp on GY (t/ha)') + ylab('-log10(p) of TE family')+ labs( color='Superfamily') + guides(alpha='none', color=guide_legend(ncol=2)) 


pdf('~/transfer/cshl_tewas.pdf', 6, 2.5)
dtsfam
gyfam
dtsfaml
gyfaml
dtsfamlv
gyfamlv
dev.off()
