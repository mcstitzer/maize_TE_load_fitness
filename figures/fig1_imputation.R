library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

## figure 1 is :
# a) genome size vs NAM parent, ranked by parent
# b) TE content vs NAM parent
# c) TE content vs Genome size
## message: we impute genome size and te content, and these are correlated but not identical

source('../figures/color_palette.R')
namfams=read.table('../figures/nam_fams.txt')
sk=read.table('../imputation/SampleToKeep.txt')

gs=read.table('../imputation/ril_bp_repeats.2022-07-15.txt', header=T, comment.char='')
gs$namRIL=substr(gs$id,1,9)
gs$namFamily=substr(gs$id,1,4)

gs$nam=namfams$V2[match(gs$namFamily, namfams$V1)]
gs$subpop=nam$subpop[match(toupper(gs$nam), toupper(nam$genome))]
gs$subpop[gs$subpop=='B73']=NA

gs=gs[gs$id %in% sk$V1,]


## need to redo this part and it's a pain...
parents=read.table('../imputation/parent_bp_repeats.2022-08-04.txt', header=T)
parents$nam=
## add parents subpop

## unify parental naming -  il14h and hp301

## size of point 2, make a and b have points instead of just name

## minmax range of GSminmax, TEminmax
GSminmax=c(min(c(gs$genomesize, parents$genomesize)), max(c(gs$genomesize, parents$genomesize)))/1e9
TEminmax=c(min(c(gs$tebp, parents$tebp)), max(c(gs$tebp, parents$tebp)))/1e9

pdf(paste0('~/transfer/fig1_imputation.', Sys.Date(), '.pdf'), 10,9)
famTE=ggplot(gs, aes(x=genomesize/1e9, y=tebp/1e9, color=subpop))+ scale_fill_manual(values=nampal, name='Subpopulation', na.translate=F) + scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ facet_wrap(~nam, strip.position='top') + 
#                             geom_label(data=parents[1,-which(colnames(parents)=='nam')], aes(label='B73'), alpha=0.8)+ 
#                             geom_hline(data=parents[1,-which(colnames(parents)=='nam')], aes(yintercept=TEbpFilter/1e9), alpha=0.8, col='gray')+ 
#                             geom_vline(data=parents[1,-which(colnames(parents)=='nam')], aes(xintercept=bpFilter/1e9), alpha=0.8, col='gray')+ 
                             geom_point() +
                             geom_point(data=parents[1,-which(colnames(parents)=='nam')], aes(x=genomesize/1e9, y=tebp/1e9), col='gray') +
                             geom_point(data=parents[-1,], aes(x=genomesize/1e9, y=tebp/1e9, fill=subpop), shape=21, col='black') +
                             geom_label_repel(data=parents[-1,], aes(label=sample), size=3, alpha=0.8, nudge_y = 0.3, nudge_x=-0.5, show.legend=F) + 
                             ylab('Imputed TE content (Gbp)') + xlab('Imputed genome size (Gbp)') +
                             scale_x_continuous(n.breaks=4, limits=GSminmax) +
                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='left', legend.justification = "center", strip.background = element_blank(),strip.text.x = element_blank())+
                             guides(colour = guide_legend(ncol = 1)) +scale_y_continuous(n.breaks = 3, limits=TEminmax)

                             
hbp=ggplot(gs, aes(y=factor(nam, levels=parents$nam[-1][order(parents$genomesize[-1])]), x=genomesize, color=subpop)) + geom_vline(xintercept=as.numeric(parents$bpFilter[1]), lty='dashed') + geom_jitter(height=0.2) + scale_color_manual(values=nampal) + xlab('All base pairs') + ylab('NAM Family') + geom_text(data=parents[-1,], inherit.aes=F, aes(y=sample, x=bpFilter, label=sample)) + scale_y_discrete(limits=parents$sample[-1][order(parents$bpFilter[-1])])
tbp=ggplot(gs, aes(y=factor(nam, levels=parents$nam[-1][order(parents$tebp[-1])]), x=tebp, color=subpop)) + geom_vline(xintercept=as.numeric(parents$TEbpFilter[1]), lty='dashed') + geom_jitter(height=0.2) + scale_color_manual(values=nampal) + xlab('TE base pairs') + ylab('NAM Family') + geom_text(data=parents[-1,], inherit.aes=F, aes(y=sample, x=TEbpFilter, label=sample)) + scale_y_discrete(limits=parents$sample[-1][order(parents$TEbpFilter[-1])])

plot_grid(plot_grid(hbp, tbp, nrow=1, labels=c('A', 'B')), famTE, ncol=1,labels=c('','C'))



dev.off()
