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
tefocusRR$mincoredist=gene$mincoredist[match(tefocusRR$term, gene$Name)]

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
              
   dtsp=ggplot(mc[-1,], aes(x=dts.estimate, y=sup, color=sup, size=supbp)) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)+ ylab('Superfamily') + xlab('Effect on DTS, per base pair')+ guides( color='none', size=guide_legend(order=2)) + labs(size='Total Superfamily\nbase pairs')
   gyrawp=ggplot(mc[-1,], aes(x=gyraw.estimate, y=sup, color=sup, size=supbp)) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)+ ylab('') + xlab('Effect on GYraw, per base pair')
   gyp=ggplot(mc[-1,], aes(x=gy.estimate, y=sup, color=sup, size=supbp)) + geom_vline(xintercept=0, lty='dashed', color='gray') + geom_point()+ scale_color_manual(values=dd.col)+ ylab('') + xlab('Effect on GY, per base pair')
   supunscaled=plot_grid(dtsp + theme(legend.position='NULL'), gyrawp+ theme(legend.position='NULL'), gyp+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
     supunscaled
   supsizelegend <- get_legend(dtsp)
              
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp, y=umrCount/totalbp, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  
   ## gene dist
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meancoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('Average Distance to Gene') + xlab('Effect on DTS, per base pair')+ guides( color='none', size=guide_legend(order=2)) + labs(size='Total Family\nbase pairs')
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meancoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GY, per base pair')
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meancoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GYraw, per base pair')
   genedistunscaled=plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  genedistunscaled
   famsizelegend <- get_legend(dtsfam)
     
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5)
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp, y=meancoredist, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5)
   plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
   ## telen gth
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meanlength, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('Average TE length') + xlab('Effect on DTS, per base pair')
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meanlength, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GY, per base pair')
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meanlength, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GYraw, per base pair')
   telengthunscaled=plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
   telengthunscaled
              
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5) + ylab('Average TE length') + xlab('DTS estimate, scaled by average bp')
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point(size=0.5) + ylab('')+ xlab('GY estimate, scaled by average bp')
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp/26, y=meanlength, color=sup)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point( size=0.5) + ylab('')+ xlab('Raw GY estimate, scaled by average bp')
telengthscaled=plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels=c('D', 'E', 'F'))
  telengthscaled
  
                 ## min gene dist
   dtsfam=ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=mincoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('Average Distance to Gene') + xlab('Effect on DTS, per base pair')
   gyfam=ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=mincoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GY, per base pair')
   gyrawfam=ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=mincoredist, color=sup, size=totalbp/26)) +geom_vline(xintercept=0, color='gray', lty='dashed') +scale_y_log10() + scale_color_manual(values=dd.col) + geom_point()+ ylab('')+ xlab('Effect on GYraw, per base pair')
   mingenedistunscaled=plot_grid(dtsfam + theme(legend.position='NULL'), gyrawfam+ theme(legend.position='NULL'), gyfam+ theme(legend.position='NULL'), ncol=3, labels='AUTO')
  mingenedistunscaled

              
   plot_grid(supscaled, telengthscaled, ncol=1)
   plot_grid(supunscaled, supsizelegend, telengthunscaled, famsizelegend, ncol=2, rel_widths=c(1,0.2,1,0.2))
   plot_grid(supunscaled, supsizelegend, genedistunscaled, famsizelegend, ncol=2, rel_widths=c(1,0.2,1,0.2))
   plot_grid(supunscaled, supsizelegend, mingenedistunscaled, famsizelegend, ncol=2, rel_widths=c(1,0.2,1,0.2))
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

              
rtf=read.table('ril_TEFam_bp_repeats.2022-08-12.txt', header=T)

rtfy=merge(teh[,c('id', 'GY', 'b73bp')], rtf, by.x='id', by.y='RIL')
tef=tidy(lm(GY~., data=rtfy[,-1]))

tef[tef$p.value<0.05,] %>% print(n=100)
              

tefocusRR=read.table('../models/te_fam_ridgeregression.2022-09-13.txt', header=T, sep='\t')

teff=merge(tef, tefocusRR)
              
teff[teff$p.value<0.05,] %>% print(n=100)
        
             
teff$sup=teff$superfam
teff$sup[teff$sup=='Copia']='RLC'
teff$sup[teff$sup=='Helitron']='DHH'
teff$sup[teff$sup=='unknown']='RLX'
teff$sup[teff$sup=='L1']='RIL'
teff$sup[teff$sup=='RTE']='RIT'
              
              
teffp=teff[teff$p.value<0.05 & substr(teff$term,1,2)!='TE' & teff$totalbp>20000000,][-1,]
 
teffd=teff[substr(teff$term,1,2)!='TE'& teff$totalbp>20000000,][-1,]
 
              
  source('../figures/color_palette.R')
library(ggridges)              
pdf('~/transfer/mm2023_tefams.pdf', 6,3)
ggplot(teff[teff$p.value<0.05,], aes(x=estimate,  y=factor(superfam),color=superfam)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+  scale_size_area(max_size=20) + xlim(-max(abs(teff[teff$p.value<0.05,]$estimate[-1]))-0.5e-8, -(min(teff[teff$p.value<0.05,]$estimate[-1])-0.5e-8))
ggplot(teff[teff$p.value<0.05,], aes(x=estimate, xend=estimate, y=superfam, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, aes(color=ifelse(pval<0.05, 'black', 'gray')), size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8)) + ylim(0,0.12)
ggplot(teff[teff$p.value<0.05,], aes(x=estimate, xend=estimate, y=superfam, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, color= 'black', size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8)) + ylim(0,0.12)

ggplot(teffp, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$estimate)), max(abs(teffp$estimate)))
ggplot(teffp, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_point(size=2)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$estimate)), max(abs(teffp$estimate)))

ggplot(teffp, aes(x=  tefocusGY, y=sup, fill=sup, color=sup)) + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$tefocusGY)), max(abs(teffp$tefocusGY)))
ggplot(teffp, aes(x=  tefocusGY, y=sup, fill=sup, color=sup)) + geom_point(size=2)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$tefocusGY)), max(abs(teffp$tefocusGY)))

              
ggplot(teffd, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$estimate)), max(abs(teffd$estimate)))
ggplot(teffd, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_point(size=2)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$estimate)), max(abs(teffd$estimate)))
  
ggplot(teffd, aes(x=  tefocusGY, y=sup, fill=sup, color=sup)) + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$tefocusGY)), max(abs(teffd$tefocusGY)))
ggplot(teffd, aes(x=  tefocusGY, y=sup, fill=sup, color=sup)) + geom_point(size=2)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$tefocusGY)), max(abs(teffd$tefocusGY)))
   
ggplot(teffd, aes(x=  tefocusGYraw, y=sup, fill=sup, color=sup)) + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$tefocusGYraw)), max(abs(teffd$tefocusGYraw)))
ggplot(teffd, aes(x=  tefocusGYraw, y=sup, fill=sup, color=sup)) + geom_point(size=2)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$tefocusGYraw)), max(abs(teffd$tefocusGYraw)))
  
              
ggplot(teffp, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=4, jittered_points=T, scale=0.7)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$estimate)), max(abs(teffp$estimate)))
ggplot(teffp, aes(x=  tefocusGY, y=sup, fill=sup, color=sup))+ geom_vline(xintercept=0, color='gray', lty='dashed') + geom_density_ridges(aes(pointcolor=sup, point_fill=sup), alpha=0.2, point_alpha=1, point_size=2, jittered_points=T, scale=0.8)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$tefocusGY)), max(abs(teffp$tefocusGY)))

ggplot(teffp, aes(x=  estimate, y=sup, fill=sup, color=sup)) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_jitter(size=4, height=0.3, alpha=0.9)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$estimate)), max(abs(teffp$estimate)))
ggplot(teffp, aes(x=  tefocusGY, y=sup, fill=sup, color=sup)) + geom_vline(xintercept=0, color='gray', lty='dashed')+ geom_jitter(size=4, height=0.3, alpha=0.9)       + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffp$tefocusGY)), max(abs(teffp$tefocusGY)))

 ggplot(teffd, aes(x=  tefocusGY))+ geom_vline(xintercept=0, color='gray', lty='dashed') + geom_density(alpha=0.2) + geom_jitter(height=1, aes(y=1, color=sup))      + scale_color_manual(values=dd.col)     + scale_fill_manual(values=dd.col)   +scale_point_color_hue(l = 40)+ xlim(-max(abs(teffd$tefocusGY)), max(abs(teffd$tefocusGY)))
             
dev.off()
              
                  
                  
ggplot(iris, aes(x = Sepal.Length, y = Species, fill = Species)) +
  geom_density_ridges(
    aes(point_color = Species, point_fill = Species, point_shape = Species),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23))
              
        
