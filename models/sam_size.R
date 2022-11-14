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
gs=read.table('../imputation/ril_bp_repeats.greaterthan22B73correct.2022-09-15.txt', header=T, comment.char='')
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
parents=read.table('../imputation/parent_bp_repeatsgreaterthan22B73correct.2022-09-16.txt', he
ader=T)
parents$nam=namfams$V2[match(toupper(str_split_fixed(parents$id, '_', 2)[,1]), toupper(namfams$V2))]
parents$subpop=nam$subpop[match(toupper(parents$nam), toupper(nam$genome))]
parents$subpop[parents$subpop=='B73']=NA
## size of point 2, make a and b have points instead of just name
## minmax range of GSminmax, TEminmax
GSminmax=c(min(c(gs$genomesize, parents$genomesize)), max(c(gs$genomesize, parents$genomesize)))/1e6
TEminmax=c(min(c(gs$tebp, parents$tebp)), max(c(gs$tebp, parents$tebp))+25*1e6)/1e6 ## bump up te max so that labels don't overlap??
head(gs)
cs=read.table('thompson2015_sammeans_rils_tables8.txt', header=T)


cs$nam=paste0(ifelse(substr(cs$line,1,1)=='Z', ifelse(substr(cs$line,1,2)=='Z2', 'Z024', 'Z005'), NA), 'E', str_pad(gsub('E', '', str_split_fixed(cs$line, '-', 2)[,2]), 4, side='left', pad=0))

gs$SAM_height=cs$SAM_height[match(gs$namRIL, cs$nam)]
gs$SAM_width=cs$SAM_width[match(gs$namRIL, cs$nam)]
gs$SAM_height.width.ratio=cs$SAM_height.width.ratio[match(gs$namRIL, cs$nam)]
gs$SAM_volume=cs$SAM_volume[match(gs$namRIL, cs$nam)]
gs$SAM_arc.length=cs$SAM_arc.length[match(gs$namRIL, cs$nam)]
gs$SAM_midpoint.width=cs$SAM_midpoint.width[match(gs$namRIL, cs$nam)]


cor.test(gs$SAM_width, gs$genomesize, use='complete.obs')
cor.test(gs$SAM_height, gs$genomesize, use='complete.obs')
cor.test(gs$SAM_height.width.ratio, gs$genomesize, use='complete.obs')
cor.test(gs$SAM_arc.length, gs$genomesize, use='complete.obs')
cor.test(gs$SAM_midpoint.width, gs$genomesize, use='complete.obs')
cor.test(gs$SAM_volume, gs$genomesize, use='complete.obs')


cor.test(gs$SAM_width, gs$tebp, use='complete.obs')
cor.test(gs$SAM_height, gs$tebp, use='complete.obs')
cor.test(gs$SAM_height.width.ratio, gs$tebp, use='complete.obs')
cor.test(gs$SAM_arc.length, gs$tebp, use='complete.obs')
cor.test(gs$SAM_midpoint.width, gs$tebp, use='complete.obs')
cor.test(gs$SAM_volume, gs$tebp, use='complete.obs')




