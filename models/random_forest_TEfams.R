library(stringr)
library(dplyr)
library(lme4)
library(lme4qtl)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)
library(glmnet)
library(glmnetUtils)
library(RColorBrewer)
library(xgboost)

source('../figures/color_palette.R')
namfams=read.table('../figures/nam_fams.txt')

sk=read.table('../imputation/SampleToKeep.txt')


## check imputation to get tecounts object... i'm being lazy

p=read.csv('../phenotypes/all_NAM_phenos.csv')
h=read.csv('../phenotypes/BLUEs.csv')

h$ID=str_split_fixed(h$genotype, '/', 2)[,1] ## pull off the PHZ51 tester

hp=readRDS('../phenotypes/NAM_H-pheno.rds')
hp$genotype=rownames(hp)
hp$ID=str_split_fixed(rownames(hp), '/', 2)[,1] ## pull off the PHZ51 tester

## guillaume did two grain yield blues - one adjusted by flowering time (I'm calling GY), and another with the raw values (I'm calling GYraw)
h$GYraw=hp$GY[match(h$ID, hp$ID)] ## 

fams=read.table('../imputation/ril_TEFam_bp_repeats.2022-08-12.txt', header=T)

fams$namRIL=substr(fams$RIL,1,9)
fams$namFamily=substr(fams$RIL,1,4)

fams$nam=namfams$V2[match(fams$namFamily, namfams$V1)]
fams$subpop=nam$subpop[match(toupper(fams$nam), toupper(nam$genome))]
fams$subpop[fams$subpop=='B73']=NA

fams$keep=fams$RIL %in% sk$V1


## these generate data frame with single measurment per ril, from the conservative keep list from cinta
teh=merge(fams[fams$keep,], h[,-c(1)], by.x='namRIL', by.y='ID')

## # columns got to be: 1=sample name, 2=pheno value, 3:end=genotypes
tefocus=teh[,c(1,313:316,3:308)]



r=sample(1:nrow(tefocus), nrow(tefocus)*0.8)
trainDTS=tefocus[r,-c(3:5)] %>% na.omit()
testDTS=tefocus[-r,-c(3:5)] %>% na.omit()
trainGY=tefocus[r,-c(2:3,5)] %>% na.omit()
testGY=tefocus[-r,-c(2:3,5)] %>% na.omit()
trainGYraw=tefocus[r,-c(2:4)] %>% na.omit()
testGYraw=tefocus[-r,-c(2:4)] %>% na.omit()


###### run ml


DTSxgb <-
  xgboost(
data = data.matrix(trainDTS[,-c(1:2)]), 
 label = trainDTS$DTS, 
    nrounds = 1000,
    objective = "reg:squarederror",
    early_stopping_rounds = 3,
    max_depth = 6,
    eta = .25
  )   
# predict values in test set
y_pred <- predict(DTSxgb, data.matrix(testDTS[,-c(1:2)]))
cor(y_pred, testDTS$DTS)^2


DTSmodel <- xgb.dump(DTSxgb, with.stats = T)
names <- colnames(trainDTS)[-c(1:2)]
# Compute feature importance matrix
DTSimportance_matrix <- xgb.importance(names, model = DTSxgb)




GYxgb <-
  xgboost(
data = data.matrix(trainGY[,-c(1:2)]), 
 label = trainGY$GY, 
    nrounds = 1000,
    objective = "reg:squarederror",
    early_stopping_rounds = 3,
    max_depth = 6,
    eta = .25
  )   
# predict values in test set
y_pred <- predict(GYxgb, data.matrix(testGY[,-c(1:2)]))
cor(y_pred, testGY$GY)^2


GYmodel <- xgb.dump(GYxgb, with.stats = T)
names <- colnames(trainGY)[-c(1:2)]
# Compute feature importance matrix
GYimportance_matrix <- xgb.importance(names, model = GYxgb)


cor(tefocus$GY, tefocus$opie_AC198924_6206, use='pairwise.complete.obs')

cor(tefocus$GY, tefocus$DTC_ZM00063_consensus, use='pairwise.complete.obs')


GYrawxgb <-
  xgboost(
data = data.matrix(trainGYraw[,-c(1:2)]), 
 label = trainGYraw$GYraw, 
    nrounds = 1000,
    objective = "reg:squarederror",
    early_stopping_rounds = 3,
    max_depth = 6,
    eta = .25
  )   
# predict values in test set
y_pred <- predict(GYrawxgb, data.matrix(testGYraw[,-c(1:2)]))
cor(y_pred, testGYraw$GYraw)^2


GYrawmodel <- xgb.dump(GYrawxgb, with.stats = T)
names <- colnames(trainGYraw)[-c(1:2)]
# Compute feature importance matrix
GYrawimportance_matrix <- xgb.importance(names, model = GYrawxgb)


#####
mbTEs=read.table('../imputation/TE_content_across_NAM_genotypes_by_fam.2022-08-08.txt', header=T)
mbTEs$Name=gsub('_LTR', '', gsub('_INT', '', mbTEs$Name))


allImportant=merge(merge(DTSimportance_matrix, GYimportance_matrix, by='Feature'), GYrawimportance_matrix, by='Feature')
colnames(allImportant)=c('Feature', 'Gain.DTS', 'Cover.DTS', 'Frequency.DTS', 'Gain.GY', 'Cover.GY', 'Frequency.GY', 
                          'Gain.GYraw', 'Cover.GYraw', 'Frequency.GYraw')


allImportant$Classification=mbTEs$Classification[match(allImportant$Feature, mbTEs$Name)]
allImportant$totalbp=mbTEs$bp[match(allImportant$Feature, mbTEs$Name)]

allImportant$superfam=str_split_fixed(allImportant$Classification, '/', 2)[,2]
allImportant$superfam[allImportant$superfam%in%c('CRM','Gypsy')]='RLG'
allImportant$superfam[allImportant$superfam%in%c('Copia')]='RLC'
allImportant$superfam[allImportant$superfam%in%c('Helitron')]='DHH'
allImportant$superfam[allImportant$superfam%in%c('L1')]='RIT'
allImportant$superfam[allImportant$superfam%in%c('unknown')]='RLX'


tecolors=c("#c04673","#826061","#c6773c","#9783b5","#6c4da4","#a34bd6","#a9ce40","#8fa657","#657563","#68c29b","#6ea5c0","#4771be","#ed3725")
names(tecolors)=c('RLC', 'RIL', 'RLX', 'RST', 'RIT', 'RLG', 					'DTA', 'DTC', 'DTH', 'DTM', 'DTX', 'DTT', 'DHH')
dd.col=tecolors  

pdf('~/transfer/TEfam_rf.pdf', 8,10)
allImportant %>% top_n(10, Gain.DTS) %>% ggplot( aes(x=Gain.DTS, y=reorder(Feature, Gain.DTS), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)
allImportant %>% top_n(10, Gain.GY) %>% ggplot(aes(x=Gain.GY, y=reorder(Feature, Gain.GY), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)
allImportant %>% top_n(10, Gain.GYraw) %>% ggplot(aes(x=Gain.GYraw, y=reorder(Feature, Gain.GYraw), fill=superfam)) + geom_bar(stat='identity')+ scale_fill_manual(values=dd.col)

allImportant %>% top_n(100, Gain.DTS) %>% ggplot( aes(x=Gain.DTS, y=reorder(Feature, Gain.DTS), fill=superfam)) + geom_bar(stat='identity')+ scale_fill_manual(values=dd.col)
allImportant %>% top_n(100, Gain.GY) %>% ggplot(aes(x=Gain.GY, y=reorder(Feature, Gain.GY), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)
allImportant %>% top_n(100, Gain.GYraw) %>% ggplot(aes(x=Gain.GYraw, y=reorder(Feature, Gain.GYraw), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)


allImportant %>% top_n(50, Gain.DTS) %>% ggplot( aes(x=Gain.DTS, y=reorder(Feature, Gain.DTS), fill=superfam)) + geom_bar(stat='identity')+ scale_fill_manual(values=dd.col)
allImportant %>% top_n(50, Gain.GY) %>% ggplot(aes(x=Gain.GY, y=reorder(Feature, Gain.GY), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)
allImportant %>% top_n(50, Gain.GYraw) %>% ggplot(aes(x=Gain.GYraw, y=reorder(Feature, Gain.GYraw), fill=superfam)) + geom_bar(stat='identity') + scale_fill_manual(values=dd.col)


dev.off()
