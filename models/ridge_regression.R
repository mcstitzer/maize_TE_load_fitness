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


## set up 
## te focus
## all repeat focus
## # columns got to be: 1=sample name, 2=pheno value, 3:end=genotypes
tefocus=teh[,c(1,313:316,3:308)]



### nam sampling approach from anju
## https://bitbucket.org/bucklerlab/haplotype_prediction/src/master/HARE_redone/src/4c_phenotype_prediction_outofpopulation_NAMsample_sampletraits.R

#function for phenotype prediction

predictPhenotype <- function(train_set, ## columns got to be: 1=sample name, 2=pheno value, 3:end=genotypes
                                             test_set,
                                             method = c("ridge", "LASSO"),
                                             nFolds = 10,
                                             returnModel=F,
                                             ...) {
  method <- method[1]
  
  if (method == "ridge") {
    message("Running ridge-regression...")
    
#   #select common x in train and test set
#   common_x <- intersect(colnames(train_set), colnames(test_set))
#   train <- train_set[, ..common_x] %>% na.omit()
#   test <- test_set[, ..common_x] %>% na.omit()
    
    # transform both sets into matrices using `data.matrix`
    x_train_mat <- data.matrix(train_set[, 3:ncol(train_set)])
    x_test_mat  <- data.matrix(test_set[, 3:ncol(test_set)])
    
    # Model the matrices using glmnet
    cv <- glmnet::cv.glmnet(
      x      = x_train_mat,
      y      = unlist(train[, 2]),
      alpha  = 0,
      nfolds = nFolds
    )
    #fit model with optimum lambda
    fit <- glmnet::glmnet(
      x      = x_train_mat,
      y      = unlist(train[, 2]),
      alpha  = 0,
      lambda = cv$lambda.min, #10^seq(4, -3, length = 100), # grid
      thresh = 1e-12
    )
    # Predict
    pred <- predict(
      object = fit,
      newx   = x_test_mat,
      type   = "response",
      s      = cv$lambda.min
    )
  } else {
    stop("Method not implemented.")
  }
  
  # Phenotype class / data frame example
  if(!returnModel){
   return(
     data.frame(
       name = test[,1],
       predicted = pred,
       observed = test[,2]
     )
     )
     }else{
       return(fit)
       }
}


error <- function(obs, pred){
  temp_df2 <- data.frame(obs, pred) %>% na.omit()
  
  rmse <- sqrt(sum((temp_df2$pred-temp_df2$obs)^2)/nrow(temp_df2))
  mae <- sum(abs(temp_df2$pred-temp_df2$obs))/nrow(temp_df2)
  return(c(rmse, mae))
}

#------------------------------------------------------- 
#cross population prediction accuracy across NAM and 282 random values
#-------------------------------------------------------

#loop starts here


r=sample(1:nrow(tefocus), nrow(tefocus)*0.8)
train=tefocus[r,-c(3:5)] %>% na.omit()
test=tefocus[-r,-c(3:5)] %>% na.omit()
tefocusDTScor=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=F)
tefocusDTS=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)
#cor.test(temp[,2], temp[,3])
train=tefocus[r,-c(2:3,5)] %>% na.omit()
test=tefocus[-r,-c(2:3,5)] %>% na.omit()
tefocusGYcor=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=F)
tefocusGY=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)
train=tefocus[r,-c(2:4)] %>% na.omit()
test=tefocus[-r,-c(2:4)] %>% na.omit()
tefocusGYrawcor=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=F)
tefocusGYraw=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)

cor.test(tefocusDTScor[,2], tefocusDTScor[,3])
cor.test(tefocusGYcor[,2], tefocusGYcor[,3])
cor.test(tefocusGYrawcor[,2], tefocusGYrawcor[,3])

b=broom::tidy(tefocusDTS)
b=broom::tidy(tefocusGY)
b=broom::tidy(tefocusGYraw)

tefocusRR=data.frame(term=broom::tidy(tefocusDTS)$term, tefocusDTS=broom::tidy(tefocusDTS)$estimate, tefocusGY=broom::tidy(tefocusGY)$estimate, tefocusGYraw=broom::tidy(tefocusGYraw)$estimate)



mbTEs=read.table('../imputation/TE_content_across_NAM_genotypes_by_fam.2022-08-08.txt', header=T)
mbTEs$Name=gsub('_LTR', '', gsub('_INT', '', mbTEs$Name))
tefocusRR$Classification=mbTEs$Classification[match(tefocusRR$term, mbTEs$Name)]
tefocusRR$totalbp=mbTEs$bp[match(tefocusRR$term, mbTEs$Name)]

tefocusRR %>%group_by(Classification) %>% dplyr::summarize(posDTS=sum(tefocusDTS>0), negDTS=sum(tefocusDTS<0), posGY=sum(tefocusGY>0), negGY=sum(tefocusGY<0))




pdf('~/transfer/ridgeregression_tefams.pdf', 14,6)
ggplot(tefocusRR[-1,], aes(x=tefocusDTS, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()
ggplot(tefocusRR[-1,], aes(x=tefocusGY, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()

ggplot(tefocusRR[-1,], aes(x=tefocusDTS*totalbp, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()
ggplot(tefocusRR[-1,], aes(x=tefocusGY*totalbp, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()
ggplot(tefocusRR[-1,], aes(x=tefocusGYraw*totalbp, y=totalbp, label=term, color=Classification)) + geom_text() +scale_y_log10()


dev.off()


pred_accuracy=c()
#run for each different sets of 250 RILs (10 RILs from each family)
#foreach(sample = 1:20) %dopar% {
#  sample_NAM_RILs <- fread(paste0(pheno_dir, '/sample_NAM_RILs/sample_NAM_RILs', sample, '.csv'))
#  
#  geno_NAM_all <- geno_NAM_all[geno_NAM_all$taxa %in% sample_NAM_RILs$V2, ]
#  geno_NAM_wona <- geno_NAM_all[,which(unlist(lapply(geno_NAM_all, function(x)!all(is.na(x))))),with=F]
#  
#  common_colnames_NAM <- intersect(colnames(geno_NAMfounders), colnames(geno_NAM_wona))
#  geno_NAM <- rbind(geno_NAMfounders[,..common_colnames_NAM], geno_NAM_wona[, ..common_colnames_NAM])
  
#  accuracy_trait <- c()
#  for (trait in c('DTS', 'GY', 'GYraw')){
    
#    pheno_geno_282 <- merge(pheno_282[,c('taxa',trait),with = FALSE], geno_282, by = 'taxa') 
#    pheno_geno_NAM <- merge(pheno_NAM[,c('taxa',trait),with = FALSE], geno_NAM, by = 'taxa')
    
#    temp_df <- predictPhenotype_crosspopulation(train_set = pheno_geno_NAM, 
 #                                               test_set = pheno_geno_282, 
 #                                               method = 'ridge')
  
  
  r=sample(1:nrow(tefocus), nrow(tefocus)*0.8)

    train=tefocus[r,-c(2:3,5)] %>% na.omit()
    test=tefocus[-r,-c(2:3,5)] %>% na.omit()
    tefocusGYcor=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=F)

#    write.csv(temp_df, paste0(result_dir, '/predictedvalues_exp_haplotypes_sampleNAMtaxa/', tissue, '_', trait, '_train_NAMsample',sample, '_test_282.csv'))
    
    ra <- cor(tefocusGYcor[,2], tefocusGYcor[,3])
    pred_accuracy <- rbind(pred_accuracy, ra)
#    message('done with trait ', trait)
    
#  }
#  rownames(accuracy_trait) <- traits
  
#  write.csv(accuracy_trait, paste0(result_dir, '/', tissue, '_prediction_accuracy_train_NAMsample', sample, '_test_282.csv'))
  
#  message('done with sample ', sample)
#}




