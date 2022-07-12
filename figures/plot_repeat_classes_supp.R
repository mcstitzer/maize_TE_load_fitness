library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)


anona=read.table('../imputation/allNAM_hapids.TEbpUpdate.sup.overlappingRRNA.2022-07-12.txt', header=T, comment.char='')

source('../figures/color_palette.R') ## this loads object nam, with cols genome and subpop
namfams=read.table('../figures/nam_fams.txt')

anona$subpop=nam$subpop[match(toupper(anona$genotype), toupper(nam$genome))]

par=anona %>% group_by(genotype, subpop) %>% summarize_at(vars(hapwidth:nonTEbp), sum, na.rm=T)

parm=melt(par, id.vars=c('genotype', 'subpop'))

pdf('~/transfer/repeat_classes.pdf',14,8)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs')
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,100e6)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,20e6)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,1e6)
dev.off()

