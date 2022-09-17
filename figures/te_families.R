library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(broom)

source('../figures/color_palette.R')
namfams=read.table('../figures/nam_fams.txt')

teh=read.table('../models/geno_pheno_gphenos.greaterthan22b73correct.2022-09-15.txt', header=T)
supTEs=read.table('../imputation/TE_content_across_NAM_genotypes_by_sup.txt', header=T)

  pah=data.frame(pheno=c('DTS', 'PH', 'GY', 'GYraw'), geno='te_sups')
models=sapply(pah$pheno, function(x) lm(x~., data=teh[,c(which(colnames(teh)==x), 6:16)]))
   dtsmod=lm(DTS~., data=teh[,c(29, 6:16)])
   gyrawmod=lm(GYraw~., data=teh[,c(32, 6:16)])
   gymod=lm(GY~., data=teh[,c(31, 6:16)])
   mc=cbind(tidy(dtsmod), cbind(tidy(gyrawmod), tidy(gymod)))  
   colnames(mc)=c('dts.term', 'dts.estimate', 'dts.se', 'dts.statistic', 'dts.p',
                  'gyraw.term', 'gyraw.estimate', 'gyraw.se', 'gyraw.statistic', 'gyraw.p',
                  'gy.term', 'gy.estimate', 'gy.se', 'gy.statistic', 'gy.p')
   mc$sup=c(NA, 'DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'RIL', 'RIT', 'RLC', 'RLG', 'RLX')
              
   mc$supbp=supTEs$bp[match(mc$sup, supTEs$sup)]
   
              
              
  


tefocusRR=read.table('../models/te_fam_ridgeregression.2022-09-13.txt', header=T, sep='\t')

len=read.table('../te_summaries/te_lengths_NAM.txt', header=T)
umr=read.table('../te_summaries/te_umr_NAM.txt', header=T)
gene=read.table('../te_summaries/te_fam_gene_dist_B73.txt', header=T)

tefocusRR$meanlength=len$meanlength[match(tefocusRR$term, len$Name)]
tefocusRR$umrCount=umr$umrCount[match(tefocusRR$term, umr$Name)]
tefocusRR$umrCount[is.na(tefocusRR$umrCount)]=0
tefocusRR$meangenedist=gene$meangenedist[match(tefocusRR$term, gene$Name)]
tefocusRR$meancoredist=gene$meancoredist[match(tefocusRR$term, gene$Name)]

classificationTE=c('DNA/DTA', 'DNA/DTC', 'DNA/DTH', 'DNA/DTM', 'DNA/DTT', 'DNA/Helitron', 'LINE/L1', 'LINE/RTE', 'LINE/unknown', 'LTR/CRM', 'LTR/Copia', 'LTR/Gypsy', 'LTR/unknown', 'MITE/DTA', 'MITE/DTC', 'MITE/DTH', 'MITE/DTM', 'MITE/DTT')

 tefocusRR$sup=str_split_fixed(tefocusRR$Classification, '/', 2)[,2]
 tefocusRR$sup[tefocusRR$sup=='Gypsy']='RLG'
 tefocusRR$sup[tefocusRR$sup=='Copia']='RLC'
 tefocusRR$sup[tefocusRR$sup=='unknown' & tefocusRR$Classification=='LTR/unknown']='RLX'
 tefocusRR$sup[tefocusRR$sup=='Helitron']='DHH'
 tefocusRR$sup[tefocusRR$sup=='CRM']='RLG' ## need these and following ones because not structural
 tefocusRR$sup[tefocusRR$sup=='L1']='RIL'
 tefocusRR$sup[tefocusRR$sup=='RTE']='RIT'
 tefocusRR$sup[tefocusRR$sup=='unknown' & tefocusRR$Classification=='LINE/unknown']='RIX'
 
 tefocusRR$sup[!tefocusRR$Classification %in% classificationTE]=NA

              
pdf(paste0('~/transfer/fig4_supAndFam.', Sys.Date(), '.pdf'),12,7)
   dtsp=ggplot(mc[-1,], aes(x=dts.estimate*supbp/26, y=sup, color=sup, size=ifelse(dts.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col) + ylab('Superfamily') + xlab('DTS estimate, scaled by average bp')
   gyrawp=ggplot(mc[-1,], aes(x=gyraw.estimate*supbp/26, y=sup, color=sup, size=ifelse(gyraw.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col) + ylab('') + xlab('Raw GY estimate, scaled by average bp')
   gyp=ggplot(mc[-1,], aes(x=gy.estimate*supbp/26, y=sup, color=sup, size=ifelse(gy.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col) + ylab('') + xlab('GY estimate, scaled by average bp')
  supscaled= plot_grid(dtsp + theme(legend.position='NULL'), gyrawp+ theme(legend.position='NULL'), gyp+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
        supscaled
              
   dtsp=ggplot(mc[-1,], aes(x=dts.estimate, y=sup, color=sup, size=ifelse(dts.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)
   gyrawp=ggplot(mc[-1,], aes(x=gyraw.estimate, y=sup, color=sup, size=ifelse(gyraw.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)
   gyp=ggplot(mc[-1,], aes(x=gy.estimate, y=sup, color=sup, size=ifelse(gy.p<0.05, 0.9, 0.3))) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)
   plot_grid(dtsp + theme(legend.position='NULL'), gyrawp+ theme(legend.position='NULL'), gyp+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
     
              
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   ## gene dist
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
   ## telen gth
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5) + ylab('Average TE length') + xlab('DTS estimate, scaled by average bp')
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5) + ylab('')+ xlab('GY estimate, scaled by average bp')
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5) + ylab('')+ xlab('Raw GY estimate, scaled by average bp')
telengthscaled=plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels=c('D', 'E', 'F'))
  telengthscaled
  
   plot_grid(supscaled, telengthscaled, ncol=1)
dev.off()
         
              

              
              
pdf('~/transfer/tefam_features.pdf')

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)


ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text()  + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text()  + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meanlength, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meangenedist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meancoredist, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount/totalbp, label=term, color=superfam)) +geom_vline(xintercept=0, color='gray', lty='dashed') + geom_text() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

dev.off()
