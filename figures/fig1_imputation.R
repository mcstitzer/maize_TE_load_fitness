library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(stringr)
library(ggrepel)
library(dplyr)

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

## remove RILs that are >3 standard deviations away from mean
famsd=gs %>% group_by(nam) %>% dplyr::summarize(sd=sd(genomesize), mean=mean(genomesize))
gs$sd=famsd$sd[match(gs$nam, famsd$nam)]
gs$mean=famsd$mean[match(gs$nam, famsd$nam)]
gs$threesd=abs(gs$mean-gs$genomesize)/gs$sd >3
gs$foursd=abs(gs$mean-gs$genomesize)/gs$sd >4

gs=gs[!gs$threesd,]


## need to redo this part and it's a pain...
parents=read.table('../imputation/parent_bp_repeats.2022-08-04.txt', header=T)
parents$nam=namfams$V2[match(toupper(str_split_fixed(parents$id, '_', 2)[,1]), toupper(namfams$V2))]
parents$subpop=nam$subpop[match(toupper(parents$nam), toupper(nam$genome))]
parents$subpop[parents$subpop=='B73']=NA


## size of point 2, make a and b have points instead of just name

## minmax range of GSminmax, TEminmax
GSminmax=c(min(c(gs$genomesize, parents$genomesize)), max(c(gs$genomesize, parents$genomesize)))/1e6
TEminmax=c(min(c(gs$tebp, parents$tebp)), max(c(gs$tebp, parents$tebp))+25*1e6)/1e6 ## bump up te max so that labels don't overlap??

pdf(paste0('~/transfer/fig1_imputation.', Sys.Date(), '.pdf'), 10,9)
famTE=ggplot(gs, aes(x=genomesize/1e6, y=tebp/1e6, color=subpop))+ scale_fill_manual(values=nampal, name='Subpopulation', na.translate=F) + scale_color_manual(values=nampal[1:5], name='Subpopulation', na.translate=F)+ facet_wrap(~nam, strip.position='top') + 
#                             geom_label(data=parents[1,-which(colnames(parents)=='nam')], aes(label='B73'), alpha=0.8)+ 
#                             geom_hline(data=parents[1,-which(colnames(parents)=='nam')], aes(yintercept=TEbpFilter/1e9), alpha=0.8, col='black')+ 
#                             geom_vline(data=parents[1,-which(colnames(parents)=='nam')], aes(xintercept=bpFilter/1e9), alpha=0.8, col='black')+ 
                             geom_point() +
                             geom_point(data=parents[1,-which(colnames(parents)=='nam')], aes(x=genomesize/1e6, y=tebp/1e6), col='black', size=2) +
                             geom_point(data=parents[-1,], aes(x=genomesize/1e6, y=tebp/1e6, fill=subpop), shape=21, col='black', size=2, show.legend=F) +
#                             geom_label_repel(data=parents[-1,], aes(label=nam), size=3, alpha=0.8, nudge_y = 0.4, nudge_x=-0.6, show.legend=F) + 
                             geom_label_repel(data=parents[-1,], aes(label=nam), size=3, alpha=0.8, nudge_y = 1500, nudge_x=-800, show.legend=F) + 
                             ylab('Imputed TE Content (Mbp)') + xlab('Imputed Genome Size (Mbp)') +
                             scale_x_continuous(n.breaks=4, limits=GSminmax) +
                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='left', legend.justification = "center", strip.background = element_blank(),strip.text.x = element_blank())+
                             guides(colour = guide_legend(ncol = 1)) +scale_y_continuous(n.breaks = 4, limits=TEminmax)

                             
genomesize=ggplot(gs, aes(y=factor(nam, levels=parents$nam[-1][order(parents$genomesize[-1])]), x=genomesize/1e6, color=subpop))  + geom_jitter(height=0.2)+ geom_vline(xintercept=as.numeric(parents$genomesize[1])/1e6, lty='dashed', color='black') + scale_color_manual(values=nampal) + scale_fill_manual(values=nampal, name='Subpopulation', na.translate=F)+ xlab('Imputed Genome Size (Mbp)') + ylab('NAM Family') + 
                geom_point(data=parents[-1,], aes(x=genomesize/1e6, y=nam, fill=subpop), shape=21, col='black', size=2) + scale_y_discrete(limits=parents$nam[-1][order(parents$genomesize[-1])]) + theme(legend.position='NULL')
tebp=ggplot(gs, aes(y=factor(nam, levels=parents$nam[-1][order(parents$tebp[-1])]), x=tebp/1e6, color=subpop))  + geom_jitter(height=0.2) + geom_vline(xintercept=as.numeric(parents$tebp[1])/1e6, lty='dashed', color='black')+ scale_color_manual(values=nampal) + scale_fill_manual(values=nampal, name='Subpopulation', na.translate=F)+ xlab('Imputed TE Content (Mbp)') + ylab('NAM Family') + 
                geom_point(data=parents[-1,], aes(x=tebp/1e6, y=nam, fill=subpop), shape=21, col='black', size=2)+ scale_y_discrete(limits=parents$nam[-1][order(parents$tebp[-1])])+ theme(legend.position='NULL')

legend <- get_legend(famTE)

plot_grid(plot_grid(plot_grid(genomesize, tebp, nrow=1, labels=c('A', 'B')), famTE + theme(legend.position='NULL'), ncol=1,labels=c('','C')), legend, rel_widths=c(1,0.2))



dev.off()
saveRDS(legend, 'subpopulation_legend.RDS')


## supplement 1 is validation of genome size vs flow cyt
## compare parent imputed length to flow cytometry to pick which is best
fc=read.table('../imputation/maize_flow_cytometry_chia.txt', header=T, sep='\t')
fc$genomesize=parents$genomesize[match(toupper(fc$line), toupper(c('B73', parents$nam[-1])))]
fc$tebp=parents$tebp[match(toupper(fc$line), toupper(c('B73', parents$nam[-1])))]

fc$subpop=nam$subpop[match(toupper(fc$line), toupper(nam$genome))]
fc$subpop[fc$subpop=='B73']=NA
## remove the tils and non-nam things
fc=fc[!is.na(fc$bp_in_assembly),]

cor.test(fc$assembly10chr, fc$standardized.genome.size, use='complete', method='spearman')
cor.test(fc$bp_in_assembly, fc$standardized.genome.size, use='complete', method='spearman')
cor.test(fc$genomesize, fc$standardized.genome.size, use='complete', method='spearman')
cor.test(fc$tebp, fc$standardized.genome.size, use='complete', method='spearman')

pdf(paste0('~/transfer/supp1_flow_cytometry_comparison.', Sys.Date(), '.pdf'),10,10)
assembly=ggplot(fc, aes(x=bp_in_assembly/1e6, y=standardized.genome.size, label=line, color=subpop)) + geom_text() + ylab('Flow cytometry (standardized to B73)')+xlab('Genome Assembly (Mb)')+ scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$standardized.genome.size~fc$bp_in_assembly))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$standardized.genome.size~fc$bp_in_assembly))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')

imputedgs=ggplot(fc, aes(x=genomesize/1e6, y=standardized.genome.size, label=line, color=subpop)) + geom_text()+ ylab('Flow cytometry (standardized to B73)')+ xlab('Imputed Genome Size (Mb)')+scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$standardized.genome.size~fc$genomesize))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$standardized.genome.size~fc$genomesize))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')

imputedte=ggplot(fc, aes(x=tebp/1e6, y=standardized.genome.size, label=line, color=subpop)) + geom_text()+ ylab('Flow cytometry (standardized to B73)')+ xlab('Imputed TE Content (Mb)')+scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$standardized.genome.size~fc$tebp))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$standardized.genome.size~fc$tebp))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')


imputedgsnice=ggplot(fc, aes(x=genomesize/1e6, y=standardized.genome.size, label=line, color=subpop))+ geom_errorbar(aes(ymin=standardized.genome.size-SE, ymax=standardized.genome.size+SE), width=.1, color='gray') + geom_text()+ ylab('Flow cytometry (standardized to B73)') +xlab('Imputed Genome Size (Mb)')+ scale_color_manual(values=c(nampal[1:5], 'gray'), name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$standardized.genome.size~fc$genomesize))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$standardized.genome.size~fc$genomesize))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')
assemblynice=ggplot(fc, aes(x=bp_in_assembly/1e6, y=standardized.genome.size, label=line, color=subpop))+ geom_errorbar(aes(ymin=standardized.genome.size-SE, ymax=standardized.genome.size+SE), width=.1, color='gray') + geom_text()+ ylab('Flow cytometry (standardized to B73)') +xlab('Genome Assembly (Mb)')+ scale_color_manual(values=c(nampal[1:5], 'gray'), name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$standardized.genome.size~fc$bp_in_assembly))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$standardized.genome.size~fc$bp_in_assembly))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')
assemblyvgs=ggplot(fc, aes(x=bp_in_assembly/1e6, y=genomesize/1e6, label=line, color=subpop))+ geom_text()+ ylab('Imputed Genome Size (Mb)') +xlab('Genome Assembly (Mb)')+ scale_color_manual(values=nampal, name='Subpopulation', na.translate=F)+ stat_smooth(geom='line', lwd=1.5, method='lm', se=F, color='#99195E', alpha=0.8)+
           annotate("text",  x=Inf, y = Inf, label = paste0(
#             '\u03B2 = ', signif(pahDF$gsEffect[pahDF$geno=='tebp' & pahDF$pheno=='GYraw']*1e6, digits=2), '\n',
             'R\U00B2 = ', signif(summary(lm(fc$genomesize~fc$bp_in_assembly))$r.squared, digits=2), 
             '\np = ', signif(summary(lm(fc$genomesize~fc$bp_in_assembly))$coefficients[2,4], digits=2)), vjust=1, hjust=1, color='#99195E')
#plot_grid(plot_grid(assemblynice + theme(legend.position='NULL'), imputedgsnice + theme(legend.position='NULL'), assemblyvgs + theme(legend.position='NULL'), labels='AUTO', nrow=3), legend, rel_widths=c(1,0.2))
 plot_grid(plot_grid(assemblynice + theme(legend.position='NULL'), imputedgsnice + theme(legend.position='NULL'), labels='AUTO', nrow=2), legend, rel_widths=c(1,0.2))
       
dev.off()




