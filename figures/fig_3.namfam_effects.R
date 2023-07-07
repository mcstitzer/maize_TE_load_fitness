library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)

## fig3
# a) 

source('../figures/color_palette.R')
legend=readRDS('subpopulation_legend.RDS')


teh=read.table('../models/geno_pheno_gphenos.greaterthan22b73correct.2022-09-15.txt', header=T)

pahFam=read.table('../models/lm_output_gphenos.namFamily.greaterthan22b73correct.2023-05-19.txt', header=T)

## we probably want families labeled, at least the sig ones!

pahFam$sigLabels=pahFam$nam
pahFam$sigLabels[pahFam$p.value>=0.05]=''

pahFam$subpop=nam$subpop[match(pahFam$nam, nam$genome)]


pahFamOrig=pahFam



pdf(paste0('~/transfer/fig3_phenoeffects.', Sys.Date(), '.pdf'), 13,3)
pahFam=pahFamOrig
## just genome size te and nonte
pahFam=pahFamOrig[pahFamOrig$term %in% c('tebpPHZ51'),]
pahFam$label='TE Content' ## haah, first one goes at bottom because of 1,2 :)
dtsRange=c(max(abs(pahFam$estimate[pahFam$pheno=='DTS'])))
gygrRange=c(max(abs(pahFam$estimate[pahFam$pheno=='GYgr'])))

dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=estimate, y=label, col=subpop,  alpha=ifelse(p.value<0.05, 0.9, 0.5), size=modelR2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one bp on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) + xlim(-dtsRange,dtsRange)
gyEs=ggplot(pahFam[pahFam$pheno=='GYgr',], aes(x=estimate, y=label, col=subpop, alpha=ifelse(p.value<0.05, 0.9, 0.5), size=modelR2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) + xlim(-gygrRange,gygrRange) 
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')
legend <- get_legend(
  # create some space to the left of the legend
  ## also the nampal is subset 2:5 because we don't want NA, and we also don't have a popcorn in this data!!!!
  ggplot(pahFam[pahFam$pheno %in% c('DTS', 'GYgr'),], aes(x=estimate, y=term, col=subpop,  alpha=ifelse(p.value<0.05, 0.9, 0.3), size=modelR2)) + geom_point() + scale_color_manual(values=nampal[2:5], na.translate=F)+ theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) +scale_size_area(max_size=20)
 + theme(legend.box.margin = margin(0, 0, 0, 30)) + guides( alpha='none',  color= guide_legend(order = 1), size=guide_legend(order=2)) + labs(color='Subpopulation', size='R2')
)
plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO'), legend, ncol=2, rel_widths=c(1,0.15))



pahFam=pahFamOrig[pahFamOrig$term %in% c('tebpPHZ51', 'allknobbpPHZ51', 'centromerebpPHZ51', 'telomerebpPHZ51', 'ribosomalbpPHZ51', 'nontenonrepeatPHZ51'),]
pahFam$label=factor(pahFam$term) ## haah, first one goes at bottom because of 1,2 :)
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=estimate, y=label, col=subpop,  alpha=ifelse(p.value<0.05, 0.9, 0.5), size=modelR2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one bp on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) + xlim(-dtsRange,dtsRange)
gyEs=ggplot(pahFam[pahFam$pheno=='GYgr',], aes(x=estimate, y=label, col=subpop, alpha=ifelse(p.value<0.05, 0.9, 0.5), size=modelR2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) + xlim(-gygrRange,gygrRange) 
dtsRange=c(max(abs(pahFam$estimate[pahFam$pheno=='DTS'])))
gygrRange=c(max(abs(pahFam$estimate[pahFam$pheno=='GYgr'])))

plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')
legend <- get_legend(
  # create some space to the left of the legend
  ## also the nampal is subset 2:5 because we don't want NA, and we also don't have a popcorn in this data!!!!
  ggplot(pahFam[pahFam$pheno %in% c('DTS', 'GYgr'),], aes(x=estimate, y=term, col=subpop,  alpha=ifelse(p.value<0.05, 0.9, 0.3), size=modelR2)) + geom_point() + scale_color_manual(values=nampal[2:5], na.translate=F)+ theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$modelR2), max(pahFam$modelR2))) +scale_size_area(max_size=20)
 + theme(legend.box.margin = margin(0, 0, 0, 30)) + guides( alpha='none',  color= guide_legend(order = 1), size=guide_legend(order=2)) + labs(color='Subpopulation', size='R2')
)
plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO'), legend, ncol=2, rel_widths=c(1,0.15))

dev.off()



#### old
pdf(paste0('~/transfer/fig3_phenoeffects.', Sys.Date(), '.pdf'), 13,3)
pahFam=pahFamOrig
## just genome size te and nonte
pahFam=pahFamOrig[pahFamOrig$geno %in% c('tebp', 'genomesize'),]
pahFam$label=factor(ifelse(pahFam$geno=='genomesize', 'Genome Size', 'TE Content'), levels=c('TE Content','Genome Size')) ## haah, first one goes at bottom because of 1,2 :)
dtsEs=ggplot(pahFam[pahFam$pheno=='DTS',], aes(x=gsEffect, y=label, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal) + theme(legend.position='NULL') + ylab('') + xlab('Effect of one bp on DTS (days)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8, box.padding = 0.5, max.overlaps = Inf, show.legend=F, direction='y') + scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) 
gyEs=ggplot(pahFam[pahFam$pheno=='GY',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2))) + scale_alpha(range=c(0.5,0.9)) +scale_size_area(max_size=20) 
gyrawEs=ggplot(pahFam[pahFam$pheno=='GYraw',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY corrected for DTS (t/ha)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2))) + scale_alpha(range=c(0.5,0.9))

#plot_grid(dtsEs, gyEs, gyrawEs, rel_widths=c(1.25,1, 1), ncol=3, labels='AUTO')
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')
legend <- get_legend(
  # create some space to the left of the legend
  ## also the nampal is subset 2:5 because we don't want NA, and we also don't have a popcorn in this data!!!!
  ggplot(pahFam[pahFam$pheno %in% c('DTS', 'GY', 'GYraw'),], aes(x=gsEffect, y=geno, col=subpop,  alpha=ifelse(pval<0.05, 0.9, 0.3), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal[2:5], na.translate=F)+ ylab('') + xlab('Effect of one Mb on PH (cm)')+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y') + theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2)))
 + theme(legend.box.margin = margin(0, 0, 0, 30)) + guides( alpha='none',  color= guide_legend(order = 1), size=guide_legend(order=2)) + labs(color='Subpopulation', size='R2')
)
#plot_grid(plot_grid(dtsEs, gyEs, gyrawEs, rel_widths=c(1.25,1, 1), ncol=3, labels='AUTO'), legend, ncol=2, rel_widths=c(1,0.15))
plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO'), legend, ncol=2, rel_widths=c(1,0.15))


## lollipop
pahFam=pahFamOrig
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop, size=r2)) + geom_point() +geom_segment() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+  scale_size_area(max_size=20) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8))
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, aes(color=ifelse(pval<0.05, 'black', 'gray')), size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8)) + ylim(0,0.12)
gylollipop=ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, color= 'black', size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8)) + ylim(0,0.12)

gylollipop

ggplot(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop, size=r2)) + geom_point() +geom_segment() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on DTS (days)')+  scale_size_area(max_size=20) + xlim(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8))
ggplot(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, aes(color=ifelse(pval<0.05, 'black', 'gray')), size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on DTS (days)')+ xlim(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8)) 
dtslollipop=ggplot(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, color= 'black', size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on DTS (days)')+ xlim(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='DTS' & pahFam$geno=='tebpPHZ51',]$gsEffect)-0.5e-8)) 
dtslollipop

plot_grid(plot_grid(dtslollipop, gylollipop, rel_widths=c(1.25,1), ncol=2, labels='AUTO'), legend, ncol=2, rel_widths=c(1,0.15))


gylollipop+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y')
dtslollipop+geom_label_repel(aes(label=sigLabels), size=3, alpha=0.8,box.padding = 0.5, max.overlaps = Inf,  show.legend=F, direction='y')

dev.off()

pdf(paste0('~/transfer/cshl_phenoeffects.', Sys.Date(), '.pdf'), 18,2.5)
plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2, labels='AUTO')
plot_grid(plot_grid(dtsEs, gyEs, rel_widths=c(1.25,1), ncol=2), legend, ncol=2, rel_widths=c(1,0.2))

dev.off()


pdf(paste0('~/transfer/fig3_phenoeffects.explore.', Sys.Date(), '.pdf'), 14,8)
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
pahFam=pahFamOrig[pahFamOrig$geno %in% c('tr1bp', 'knob180bp', 'allknobbp', 'telomerebp', 'centromerebp', 'ribosomalbp'),]
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


pdf(paste0('~/transfer/mm2023_yieldeffects.', Sys.Date(), '.pdf'), 8,2)

gyEs=ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+  scale_size_area(max_size=20) + scale_alpha(range=c(0.5,0.9)) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))
gyEs

## mo18w first
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp' & pahFam$nam=='Mo18W',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+  scale_size_area(max_size=20) + scale_alpha(range=c(0.5,0.9)) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))

## then nc350
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp' & pahFam$nam%in%c('Mo18W', 'NC350'),], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+  scale_size_area(max_size=20) + scale_alpha(range=c(0.5,0.9)) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))

## then b97
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp' & pahFam$nam%in%c('Mo18W', 'NC350', 'B97'),], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+  scale_size_area(max_size=20) + scale_alpha(range=c(0.5,0.9)) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))


## get a legend?
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_point() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal, guide = "none")+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+  scale_size_area(max_size=20) + scale_alpha(range=c(0.5,0.9)) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))


gyEs=ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, y=label, col=subpop, alpha=ifelse(pval<0.05, 0.9, 0.5), size=r2)) + geom_text(aes(label=nam)) + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ theme(axis.text.y=element_blank())+ scale_size_continuous(limits=c(min(pahFam$r2), max(pahFam$r2))) + scale_alpha(range=c(0.5,0.9))
gyEs

## lollipop
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop, size=r2)) + geom_point() +geom_segment() + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+  scale_size_area(max_size=20) + xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8))
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, aes(color=ifelse(pval<0.05, 'black', 'gray')), size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8)) + ylim(0,0.12)
ggplot(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',], aes(x=gsEffect, xend=gsEffect, y=r2, yend=0, col=subpop,fill=subpop))+geom_segment(size=2) + geom_point(pch=21, color= 'black', size=6)   + geom_vline(xintercept=0, color='gray', lty='dashed')+ scale_color_manual(values=nampal)+ scale_fill_manual(values=nampal)+ theme(legend.position='NULL')+ ylab('') + xlab('Effect of one bp on GY (t/ha)')+ xlim(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8, -(min(pahFam[pahFam$pheno=='GY' & pahFam$geno=='tebp',]$gsEffect)-0.5e-8)) + ylim(0,0.12)


dev.off()


pahDF=read.table('../models/lm_output_gphenos.greaterthan22b73correct.2022-09-15.txt', header=T)

teh=read.table('../models/geno_pheno_gphenos.greaterthan22b73correct.2022-09-15.txt', header=T)

legend=readRDS('subpopulation_legend.RDS')

pdf(paste0('~/transfer/mm2023_yieldbyfam.', Sys.Date(), '.pdf'), 3,3)

for(nam in unique(teh$nam)){
  
  
gyrawt=ggplot(teh[teh$nam==nam,], aes(x=tebp/1e6, y=GYraw, color=subpop)) + geom_point(alpha=0.7) + stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8) + scale_color_manual(values=nampal) + theme(legend.position='none') + ylab('Grain Yield (BLUE, t/ha)')+xlab('Imputed TE content (Mbp)')+ xlim(min(teh$tebp/1e6), max(teh$tebp/1e6)) + ylim(min(teh$GYraw, na.rm=T), max(teh$GYraw, na.rm=T)) +
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(pahFam$r2[pahFam$geno=='tebp' & pahFam$pheno=='GYraw' & pahFam$nam==nam], digits=2), 
             '\np = ', signif(pahFam$pval[pahFam$geno=='tebp' & pahFam$pheno=='GYraw' & pahFam$nam==nam], digits=2)), vjust=1, hjust=1, color='#99195E') + ggtitle(nam)
print(gyrawt)


  }

  dev.off()
