library(corrplot)


teh=read.table('../imputation/ril_tecategories_bp.2023-04-10.txt', header=T)

teh2=read.table('../imputation/ril_bp_repeats.greaterthan22B73correct.2022-09-15.txt', header=T)



te=merge(teh, teh2, by.x='RIL', by.y='id', all=T)
te$knob180bp=NULL
te$tr1bp=NULL

pdf('~/transfer/correlations_tecats.pdf',  8,8)
corrplot(cor(te[,-1], use='pairwise.complete.obs'))
corrplot(cor(te[,-1], method='spearman', use='pairwise.complete.obs'))
dev.off()
