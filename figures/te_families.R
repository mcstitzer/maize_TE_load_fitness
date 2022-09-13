
tefocusRR=read.table('te_fam_ridgeregression.2022-09-13.txt', header=T)

len=read.table('../te_summaries/te_lengths_NAM.txt', header=T)
umr=read.table('../te_summaries/te_umr_NAM.txt', header=T)
gene=read.table('../te_summaries/te_fam_gene_dist_B73.txt', header=T)

tefocusRR$meanlength=len$meanlength[match(tefocusRR$term, len$Name)]
tefocusRR$umrCount=umr$umrCount[match(tefocusRR$term, umr$Name)]
tefocusRR$umrCount[is.na(tefocusRR$umrCount)]=0
tefocusRR$meangenedist=gene$meangenedist[match(tefocusRR$term, gene$Name)]

pdf('~/transfer/tefam_features.pdf')

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=totalbp, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=totalbp, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=totalbp, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meanlength, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meanlength, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meanlength, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=umrCount, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=umrCount, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=umrCount, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=meangenedist, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=meangenedist, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=meangenedist, label=term, color=superfam)) + geom_text() +scale_y_log10() + scale_color_brewer(palette='Set1') + geom_point(color='black', size=0.5)

dev.off()
