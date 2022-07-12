library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)


anona=read.table('../imputation/allNAM_hapids.TEbpUpdate.sup.overlappingRRNA.2022-07-12.txt', header=T, comment.char='')
keep=read.table('refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt', header=T)

source('../figures/color_palette.R') ## this loads object nam, with cols genome and subpop
namfams=read.table('../figures/nam_fams.txt')

anona$keep=anona$refrange %in% keep$refrange[keep$KEEPfinalFilter]

anona$subpop=nam$subpop[match(toupper(anona$genotype), toupper(nam$genome))]

par=anona %>% group_by(genotype, subpop) %>% summarize_at(vars(hapwidth:nonTEbp), sum, na.rm=T)
parm=melt(par, id.vars=c('genotype', 'subpop'))

park=anona[anona$keep,] %>% group_by(genotype, subpop) %>% summarize_at(vars(hapwidth:nonTEbp), sum, na.rm=T)
parkm=melt(park, id.vars=c('genotype', 'subpop'))


pdf('~/transfer/repeat_classes.pdf',14,8)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs')
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,100e6)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,20e6)
ggplot(parm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - \nTHESE DO NOT HAVE LONG RR REMOVED)\n stop worrying') + ylab('base pairs') + ylim(0,1e6)

ggplot(parkm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - reconsider long rr removal after imputing!!!') + ylab('base pairs')
ggplot(parkm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - reconsider long rr removal after imputing!!!') + ylab('base pairs') + ylim(0,100e6)
ggplot(parkm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - reconsider long rr removal after imputing!!!') + ylab('base pairs') + ylim(0,20e6)
ggplot(parkm, aes(x=factor(genotype, levels=par$genotype[order(par$hapwidth)]), y=value, color=variable)) + geom_point()  + xlab('Parent (sorted by imputed genome - reconsider long rr removal after imputing!!!') + ylab('base pairs') + ylim(0,1e6)



dev.off()

