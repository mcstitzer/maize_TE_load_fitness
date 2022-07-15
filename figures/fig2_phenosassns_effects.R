library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)


source('../figures/color_palette.R')


teh=read.table('../models/geno_pheno_gphenos.2022-07-14.txt', header=T)

pahDF=read.table('../models/lm_output_gphenos.2022-07-14.txt ', header=T)
pahFam=read.table('../models/lm_output_gphenos.namFamily.2022-07-14.txt', header=T)
pahFam$effect=pahFam$gsEffect ## go back and change this in generating file in models dir
pahFam$subpop=nam$subpop[match(toupper(pahFam$nam), toupper(nam$genome))] ## also go back and do this there

pdf(paste0('~/transfer/fig2_phenoassns.', Sys.Date(), '.pdf'), 12,8)

dts=ggplot(teh, aes(x=genomesize/1e9, y=DTS, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.825,1.9) + ylab('Days to Silking (BLUE, days)') + xlab('Imputed genome size (Gbp)') 

gy=ggplot(teh, aes(x=genomesize/1e9, y=GY, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.825,1.9) + ylab('Grain Yield (BLUE, t/ha)') + xlab('Imputed genome size (Gbp)') 

ph=ggplot(teh, aes(x=genomesize/1e9, y=PH, color=subpop))+ geom_point(alpha=0.7)  + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.825,1.9) + ylab('Plant Height (BLUE, cm)') + xlab('Imputed genome size (Gbp)') 

genomesize=plot_grid(dts, gy, ph, ncol=3, labels='AUTO')
                             
dtst=ggplot(teh, aes(x=tebp/1e9, y=DTS, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.55,1.65)+ ylab('Days to Silking (BLUE, days)')+xlab('Imputed TE content (Gbp)')

gyt=ggplot(teh, aes(x=tebp/1e9, y=GY, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.55,1.65)+ ylab('Grain Yield (BLUE, t/ha)')+xlab('Imputed TE content (Gbp)')

pht=ggplot(teh, aes(x=tebp/1e9, y=PH, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='gray', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + xlim(1.55,1.65)+ ylab('Plant Height (BLUE, cm)')+xlab('Imputed TE content (Gbp)')

tecontent=plot_grid(dtst, gyt, pht, ncol=3, labels=c('D', 'E', 'F'))

plot_grid(genomesize, tecontent, ncol=1)

dtsE=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one bp on DTS')
gyE=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY')
phE=ggplot(pahFam[pahFam$pheno=='PH',], aes(x=effect, y=geno, col=subpop, shape=pval<0.05, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on PH')

effects=plot_grid(dtsE, gyE, phE, ncol=3, labels=c('G', 'H', 'I'))

plot_grid(genomesize, tecontent, effects,ncol=1, rel_heights=c(1,1,0.5))

dev.off()

pdf(paste0('~/transfer/fig3_phenoeffects.', Sys.Date(), '.pdf'), 14,3)

## or effects as separate figure, to emphasize a totally different thing b/c by family?
plot_grid(dtsE, gyE, phE, ncol=3, labels='AUTO')

## we probably want families labeled, at least the sig ones!

pahFam$sigLabels=pahFam$nam
pahFam$sigLabels[pahFam$pval>=0.05]=''

dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=effect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
phEs=ggplot(pahFam[pahFam$pheno=='PH',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
plot_grid(dtsEs, gyEs, phEs, rel_widths=c(1.25,1,1), ncol=3, labels='AUTO')

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(pahFam[pahFam$pheno=='PH',], aes(x=effect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
)

plot_grid(plot_grid(dtsEs, gyEs, phEs, rel_widths=c(1.25,1,1), ncol=3, labels='AUTO'), legend, rel_widths=c(1,0.1))

dev.off()      




### add in ril inbred phenos from merritt

pahPDF=read.table('../models/lm_output_pphenos.2022-07-14.txt ', header=T)
tep=read.table('../models/geno_pheno_pphenos.2022-07-14.txt', header=T)


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


