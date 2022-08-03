library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)

## fig3
# a) 

source('../figures/color_palette.R')


teh=read.table('../models/geno_pheno_gphenos.2022-08-02.txt', header=T)

pahFam=read.table('../models/lm_output_gphenos.namFamily.2022-08-02.txt', header=T)

## we probably want families labeled, at least the sig ones!

pahFam$sigLabels=pahFam$nam
pahFam$sigLabels[pahFam$pval>=0.05]=''

pahFam$subpop=nam$subpop[match(pahFam$nam, nam$genome)]


pahFamOrig=pahFam

pdf(paste0('~/transfer/fig3_phenoeffects.', Sys.Date(), '.pdf'), 14,8)
pahFam=pahFamOrig
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
)

plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=1, labels='AUTO'), legend, rel_widths=c(1,0.1))

## just knob and repeats
pahFam=pahFamOrig[pahFamOrig$geno %in% c('tr1bp', 'knob180bp', 'allknobbp', 'telomerebp'),]
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
)

plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=1, labels='AUTO'), legend, rel_widths=c(1,0.1))

## just bigish TEs
pahFam=pahFamOrig[pahFamOrig$geno %in% c('rlxbp', 'rlgbp', 'rlcbp', 'dtmbp', 'dthbp', 'dtcbp', 'dtabp', 'dhhbp'),]
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
)

plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=1, labels='AUTO'), legend, rel_widths=c(1,0.1))

## just genome size te and nonte
pahFam=pahFamOrig[pahFamOrig$geno %in% c('tebp', 'genomesize', 'nontebp'),]
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one Mb on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one Mb on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect*1e6, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 12)) + guides(color='none', alpha='none') + labs(size='R2')
)

plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=1, labels='AUTO'), legend, rel_widths=c(1,0.1))



dev.off()      


