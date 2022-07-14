library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)


source('../figures/color_palette.R')

source('../figures/color_palette.R')
namfams=read.table('../figures/nam_fams.txt')

rilsumKeep=read.table('ril_bp_repeats.2022-07-14.txt', header=T)
old=read.table('teh_2022-04-14.txt', header=T)

ro=merge(rilsumKeep, old, by.x='id', by.y='sample')

cor.test(ro$tebp, ro$TEbpFilter)

sk=read.table('SampleToKeep.txt')

pdf('~/transfer/oldnew.pdf', 8,4)
ggplot(ro, aes(x=TEbpFilter/1e9, y=tebp/1e9, color=subpop)) + geom_point(alpha=0.7) +scale_color_manual(values=nampal) + theme(legend.position='none') + xlab('Old TE Mbp') + ylab('New TE Mbp')
ggplot(ro[ro$id %in% sk$V1,], aes(x=TEbpFilter/1e9, y=tebp/1e9, color=subpop)) + geom_point(alpha=0.7) +scale_color_manual(values=nampal) + theme(legend.position='none') + xlab('Old TE Mbp') + ylab('New TE Mbp')

ggplot(ro[ro$id %in% sk$V1,], aes(x=tebp/1e9, y=GY, color=subpop)) + geom_point(alpha=0.7) +scale_color_manual(values=nampal) + theme(legend.position='none') + xlab('New TE Mbp') + ylab('GY') +stat_smooth(method='lm', se=F)

dev.off()


