

                       
source('../figures/color_palette.R')

gs=read.table('../imputation/ril_bp_repeats.2022-07-14.txt', header=T, comment.char='')

## need to redo this part and it's a pain...
parents=read.table('../imputation/nam_parents_totalbp.2022-03-28.txt', header=T)

pdf('~/transfer/fig1_imputation.maizemeeting.pdf', 10/1.3,5/1.3)
famTE=ggplot(psg, aes(x=bpFilter/1e9, y=TEbpFilter/1e9, color=subpop))+ scale_fill_manual(values=nampal, name='Subpopulation', na.translate=F) + scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ facet_wrap(~nam, strip.position='top') + 
#                             geom_label(data=parents[1,-which(colnames(parents)=='nam')], aes(label='B73'), alpha=0.8)+ 
#                             geom_hline(data=parents[1,-which(colnames(parents)=='nam')], aes(yintercept=TEbpFilter/1e9), alpha=0.8, col='gray')+ 
#                             geom_vline(data=parents[1,-which(colnames(parents)=='nam')], aes(xintercept=bpFilter/1e9), alpha=0.8, col='gray')+ 
                             geom_point() +
                             geom_point(data=parents[1,-which(colnames(parents)=='nam')], aes(x=bpFilter/1e9, y=TEbpFilter/1e9), col='gray') +
                             geom_point(data=parents[-1,], aes(x=bpFilter/1e9, y=TEbpFilter/1e9, fill=subpop), shape=21, col='black') +
                             geom_label_repel(data=parents[-1,], aes(label=sample), size=3, alpha=0.8, nudge_y = 0.3, nudge_x=-0.5, show.legend=F) + 
                             ylab('Imputed TE content (Gbp)') + xlab('Imputed genome size (Gbp)') +
                             scale_x_continuous(n.breaks=4, limits=c(1.825,1.9)) +
                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='left', legend.justification = "center", strip.background = element_blank(),strip.text.x = element_blank())+
                             guides(colour = guide_legend(ncol = 1)) +scale_y_continuous(n.breaks = 3, limits=c(1.55,1.64))
famTE
                             
ggplot(psg, aes(x=bpFilterB73, y=TEbpFilterB73, color=subpop)) + geom_point() + scale_color_manual(values=nampal)+ facet_wrap(~nam) + 
#                             geom_label(data=parents[-1,], aes(label=sample), alpha=0.8) + 
#                             geom_label(data=parents[1,-which(colnames(parents)=='nam')], aes(label='B73'), alpha=0.8)+ 
                             ylab('Imputed TE content (bp)') + xlab('Imputed genome size (bp)')

dev.off()
