library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)


source('../figures/color_palette.R')


teh=read.table('../models/geno_pheno_gphenos.2022-08-31.txt', header=T)

pahDF=read.table('../models/lm_output_gphenos.2022-08-31.txt', header=T)
pahFam=read.table('../models/lm_output_gphenos.namFamily.2022-08-31.txt', header=T)
pahFam$effect=pahFam$gsEffect ## go back and change this in generating file in models dir
pahFam$subpop=nam$subpop[match(toupper(pahFam$nam), toupper(nam$genome))] ## also go back and do this there


legend=readRDS('subpopulation_legend.RDS')


pdf(paste0('~/transfer/fig2_phenoassns.', Sys.Date(), '.pdf'), 15,8)

dts=ggplot(teh, aes(x=genomesize/1e6, y=DTS, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Days to Silking (BLUE, days)') + xlab('Imputed genome size (Mbp)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='genomesize' & pahDF$pheno=='DTS']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='genomesize' & pahDF$pheno=='DTS'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='genomesize' & pahDF$pheno=='DTS'], digits=2)), vjust=1, hjust=1, color='#99195E')

gy=ggplot(teh, aes(x=genomesize/1e6, y=GY, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Grain Yield Corrected for\nFlowering Time (BLUE, t/ha)') + xlab('Imputed genome size (Mbp)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='genomesize' & pahDF$pheno=='GY']*1e6, digits=2), '\m',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='genomesize' & pahDF$pheno=='GY'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='genomesize' & pahDF$pheno=='GY'], digits=2)), vjust=1, hjust=1, color='#99195E')

gyraw=ggplot(teh, aes(x=genomesize/1e6, y=GYraw, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Grain Yield (BLUE, t/ha)') + xlab('Imputed genome size (Mbp)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#            '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='genomesize' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='genomesize' & pahDF$pheno=='GYraw'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='genomesize' & pahDF$pheno=='GYraw'], digits=2)), vjust=1, hjust=1, color='#99195E')


ph=ggplot(teh, aes(x=genomesize/1e6, y=PH, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Plant Height (BLUE, cm)') + xlab('Imputed genome size (Mbp)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='genomesize' & pahDF$pheno=='PH']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='genomesize' & pahDF$pheno=='PH'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='genomesize' & pahDF$pheno=='PH'], digits=2)), vjust=1, hjust=1, color='#99195E')

genomesize=plot_grid(dts, gyraw, gy, ncol=3, labels='AUTO')
                             
dtst=ggplot(teh, aes(x=tebp/1e6, y=DTS, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + ylab('Days to Silking (BLUE, days)')+xlab('Imputed TE content (Mbp)')+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='DTS']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='tebp' & pahDF$pheno=='DTS'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='tebp' & pahDF$pheno=='DTS'], digits=2)), vjust=1, hjust=1, color='#99195E')

gyt=ggplot(teh, aes(x=tebp/1e6, y=GY, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + ylab('Grain Yield Corrected for\nFlowering Time (BLUE, t/ha)')+xlab('Imputed TE content (Mbp)')+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GY']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='tebp' & pahDF$pheno=='GY'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='tebp' & pahDF$pheno=='GY'], digits=2)), vjust=1, hjust=1, color='#99195E')


gyrawt=ggplot(teh, aes(x=tebp/1e6, y=GY, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + ylab('Grain Yield (BLUE, t/ha)')+xlab('Imputed TE content (Mbp)')+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='tebp' & pahDF$pheno=='GYraw'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='tebp' & pahDF$pheno=='GYraw'], digits=2)), vjust=1, hjust=1, color='#99195E')

pht=ggplot(teh, aes(x=tebp/1e6, y=PH, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='gray', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')+ ylab('Plant Height (BLUE, cm)')+xlab('Imputed TE content (Mbp)')+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='PH']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahDF$r2[pahDF$geno=='tebp' & pahDF$pheno=='PH'], digits=2), 
             '\np = ', signif(pahDF$pval[pahDF$geno=='tebp' & pahDF$pheno=='PH'], digits=2)), vjust=1, hjust=1, color='#99195E')

tecontent=plot_grid(dtst, gyrawt, gyt, ncol=3, labels=c('D', 'E', 'F'))

plot_grid(genomesize, tecontent, ncol=1)


plot_grid(plot_grid(genomesize, tecontent, ncol=1), legend, rel_widths=c(1,0.2))



dtsE=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one bp on DTS')
gyE=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY corrected for DTS')
gyrawE=ggplot(pahFam[pahFam$pheno=='GYraw',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY')

effects=plot_grid(dtsE, gyrawE, gyE, ncol=3, labels=c('G', 'H', 'I'))

plot_grid(genomesize, tecontent, effects,ncol=1, rel_heights=c(1,1,0.5))

dev.off()



pdf(paste0('~/transfer/supp2_phenovpheno.', Sys.Date(), '.pdf'), 12,5)

dtsgy=ggplot(teh, aes(x=DTS, y=GY, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + xlab('Days to Silking (BLUE, days)') + ylab('Grain Yield Corrected for\nFlowering Time (BLUE, t/ha)')+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(teh$GY~teh$DTS))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(teh$GY~teh$DTS))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')


dtsgyraw=ggplot(teh, aes(x=DTS, y=GYraw, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Grain Yield (BLUE, t/ha)') + xlab('Days to Silking (BLUE, days)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(teh$GYraw~teh$DTS))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(teh$GYraw~teh$DTS))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')


gygyraw=ggplot(teh, aes(x=GY, y=GYraw, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Grain Yield (BLUE, t/ha)') + xlab('Grain Yield Corrected for\nFlowering Time (BLUE, t/ha)') +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(teh$GYraw~teh$GY))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(teh$GYraw~teh$GY))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')


phenocorr=plot_grid(dtsgyraw, dtsgy, gygyraw, ncol=3, labels='AUTO')
plot_grid(phenocorr, legend, rel_widths=c(1,0.2))

dev.off()


#pdf(paste0('~/transfer/fig3_phenoeffects.', Sys.Date(), '.pdf'), 14,3)
#
### or effects as separate figure, to emphasize a totally different thing b/c by family?
#plot_grid(dtsE, gyE, phE, ncol=3, labels='AUTO')
#
### we probably want families labeled, at least the sig ones!
#
#pahFam$sigLabels=pahFam$nam
#pahFam$sigLabels[pahFam$pval>=0.05]=''
#
#dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
#gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=effect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
#phEs=ggplot(pahFam[pahFam$pheno=='PH',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
#plot_grid(dtsEs, gyEs, phEs, rel_widths=c(1.25,1,1), ncol=3, labels='AUTO')
#
#legend <- get_legend(
#  # create some space to the left of the legend
#  ggplot(pahFam[pahFam$pheno=='PH',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
# + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
#)
#
#plot_grid(plot_grid(dtsEs, gyEs, phEs, rel_widths=c(1.25,1,1), ncol=3, labels='AUTO'), legend, rel_widths=c(1,0.1))
#
#dev.off()      




### add in ril inbred phenos from merritt

pahPDF=read.table('../models/lm_output_mphenos.2022-08-31.txt', header=T)
tep=read.table('../models/geno_pheno_mphenos.2022-08-31.txt', header=T)


pdf('~/transfer/flowering_inbreds.pdf')

ggplot(pahPDF[grepl('GDD_DTS', pahPDF$pheno) & pahPDF$geno=='genomesize',], aes(x=intercept, y=gsEffect)) + ylim(700,1400) + xlim(min(tep$genomesize), max(tep$genomesize)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + geom_point() + xlab('Genome size (bp)') + ylab('DTS Diff GDD phenos')
#  geom_abline(data=df,aes(slope=slope,intercept=intercept,color=factor(wf)))
ggplot(pahPDF[grepl('GDD_DTA', pahPDF$pheno) & pahPDF$geno=='genomesize',], aes(x=intercept, y=gsEffect, slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + ylim(700,1400) + xlim(min(tep$genomesize), max(tep$genomesize)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05))+ geom_point() + xlab('Genome size (bp)') + ylab('DTA Diff GDD phenos')

ggplot(pahPDF[grepl('GDD_DTS', pahPDF$pheno) & pahPDF$geno=='tebp',], aes(x=intercept, y=gsEffect)) + ylim(700,1400) + xlim(min(tep$tebp), max(tep$tebp)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + geom_point() + xlab('te bp (bp)') + ylab('DTS Diff GDD phenos')
#  geom_abline(data=df,aes(slope=slope,intercept=intercept,color=factor(wf)))
ggplot(pahPDF[grepl('GDD_DTA', pahPDF$pheno) & pahPDF$geno=='tebp',], aes(x=intercept, y=gsEffect, slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + ylim(700,1400) + xlim(min(tep$tebp), max(tep$tebp)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05))+ geom_point() + xlab('te bp (bp)') + ylab('DTA Diff GDD phenos')

### nongdd
ggplot(pahPDF[grepl('^DTS', pahPDF$pheno) & pahPDF$geno=='genomesize',], aes(x=intercept, y=gsEffect)) + ylim(60,92) + xlim(min(tep$genomesize), max(tep$genomesize)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + geom_point() + xlab('Genome size (bp)') + ylab('DTS Diff phenos') + geom_abline(data=pahDF[pahDF$geno=='genomesize' & pahDF$pheno=='DTS',], aes(intercept=intercept, slope=gsEffect), color='red')
#  geom_abline(data=df,aes(slope=slope,intercept=intercept,color=factor(wf)))
ggplot(pahPDF[grepl('^DTA', pahPDF$pheno) & pahPDF$geno=='genomesize',], aes(x=intercept, y=gsEffect, slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + ylim(60,92) + xlim(min(tep$genomesize), max(tep$genomesize)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05))+ geom_point() + xlab('Genome size (bp)') + ylab('DTA Diff phenos')

ggplot(pahPDF[grepl('^DTS', pahPDF$pheno) & pahPDF$geno=='tebp',], aes(x=intercept, y=gsEffect)) + ylim(60,92) + xlim(min(tep$tebp), max(tep$tebp)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + geom_point() + xlab('te bp (bp)') + ylab('DTS Diff phenos')+ geom_abline(data=pahDF[pahDF$geno=='tebp' & pahDF$pheno=='DTS',], aes(intercept=intercept, slope=gsEffect), color='red')
#  geom_abline(data=df,aes(slope=slope,intercept=intercept,color=factor(wf)))
ggplot(pahPDF[grepl('^DTA', pahPDF$pheno) & pahPDF$geno=='tebp',], aes(x=intercept, y=gsEffect, slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05)) + ylim(60,92) + xlim(min(tep$tebp), max(tep$tebp)) + geom_abline(inherit.aes=T, aes(slope=gsEffect, intercept=intercept, color=r2, lty=pval<0.05))+ geom_point() + xlab('te bp (bp)') + ylab('DTA Diff phenos')



dev.off()
  


## do all the genotypes and all the phenotypes!!!
pdf(paste0('~/transfer/ALL_phenoassns.', Sys.Date(), '.pdf'), 12,8)

for(i in colnames(teh)[3:20]){

dts=ggplot(teh, aes_string(x=i, y='DTS', color='subpop'))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Days to Silking (BLUE, days)') + xlab(i) 

gy=ggplot(teh, aes_string(x=i, y='GY', color='subpop'))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + ylab('Grain Yield (BLUE, t/ha)') + xlab(i) 

ph=ggplot(teh, aes_string(x=i, y='PH', color='subpop'))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none')  + ylab('Plant Height (BLUE, cm)') + xlab(i) 

genomesize=plot_grid(dts, gy, ph, ncol=3, labels='AUTO')
print(genomesize)
}
dev.off()


