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
train=tefocus[r,] %>% na.omit()
dim(train)
test=tefocus[-r,] %>% na.omit()
tefocusDTS=predictPhenotype(train_set=train[,-c(3:5)], test_set=test[,-c(3:5)], method='ridge', returnModel=T)
#cor.test(temp[,2], temp[,3])
tefocusGY=predictPhenotype(train_set=train[,-c(2:3,5)], test_set=test[,-c(2:3,5)], method='ridge', returnModel=T)
tefocusGYraw=predictPhenotype(train_set=train[,-c(2:4)], test_set=test[,-c(2:4)], method='ridge', returnModel=T)


b=broom::tidy(tefocusDTS)
b=broom::tidy(tefocusGY)
b=broom::tidy(tefocusGYraw)

tefocusRR=data.frame(term=broom::tidy(tefocusDTS)$term, tefocusDTS=broom::tidy(tefocusDTS)$estimate, tefocusGY=broom::tidy(tefocusGY)$estimate, tefocusGYraw=broom::tidy(tefocusGYraw)$estimate)

## repeat focus

r=sample(1:nrow(repeatfocus), nrow(repeatfocus)*0.8)
train=repeatfocus[r,] %>% na.omit()
dim(train)
test=repeatfocus[-r,] %>% na.omit()
trainDTS=train[,-c(3:5)]
testDTS=test[,-c(3:5)]
repeatfocusDTS=predictPhenotype(train_set=trainDTS, test_set=testDTS, method='ridge', returnModel=T)
#cor.test(temp[,2], temp[,3])
repeatfocusGY=predictPhenotype(train_set=train[,-c(2:3,5)], test_set=test[,-c(2:3,5)], method='ridge', returnModel=T)
repeatfocusGYraw=predictPhenotype(train_set=train[,-c(2:4)], test_set=test[,-c(2:4)], method='ridge', returnModel=T)


b=broom::tidy(repeatfocusDTS)
b=broom::tidy(repeatfocusGY)
b=broom::tidy(repeatfocusGYraw)

repeatfocusRR=data.frame(term=broom::tidy(repeatfocusDTS)$term, repeatfocusDTS=broom::tidy(repeatfocusDTS)$estimate, repeatfocusGY=broom::tidy(repeatfocusGY)$estimate, repeatfocusGYraw=broom::tidy(repeatfocusGYraw)$estimate)




## try with full matrix of each repeat class

sk=read.table('../imputation/SampleToKeep.txt')
h=read.csv('../phenotypes/BLUEs.csv')
h$ID=str_split_fixed(h$genotype, '/', 2)[,1] ## pull off the PHZ51 tester
hp=readRDS('../phenotypes/NAM_H-pheno.rds')
hp$genotype=rownames(hp)
hp$ID=str_split_fixed(rownames(hp), '/', 2)[,1] ## pull off the PHZ51 tester
## guillaume did two grain yield blues - one adjusted by flowering time (I'm calling GY), and another with the raw values (I'm calling GYraw)
h$GYraw=hp$GY[match(h$ID, hp$ID)] ## 

a=readRDS('../imputation/matrixList_bpEachRefRange.RDS')
tempMat=a[[1]]
tempMat=data.frame(tempMat[rownames(tempMat) %in% sk$V1,])
tempMat$namRIL=substr(rownames(tempMat),1,9)

newMat=merge(h[,-c(1)], tempMat, by.x='ID', by.y='namRIL')
             
## oh dang i've got a lot of missing data... mean impute (should fix to mean of two parents!)
for(i in 6:ncol(newMat)){
  newMat[is.na(newMat[,i]), i] <- mean(newMat[,i], na.rm = TRUE)
}

r=sample(1:nrow(newMat), nrow(newMat)*0.8)

train=newMat[r,-c(3:5)] %>% na.omit()
dim(train)
test=newMat[-r,-c(3:5)] %>% na.omit()
dts=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)
#cor.test(dts[,2], dts[,3])

## second 
train=newMat[r,-c(2:3,5)] %>% na.omit()
dim(train)
test=newMat[-r,-c(2:3,5)] %>% na.omit()
gy=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)
#cor.test(dts[,2], dts[,3])

## third
train=newMat[r,-c(2:4)] %>% na.omit()
dim(train)
test=newMat[-r,-c(2:4)] %>% na.omit()
gyraw=predictPhenotype(train_set=train, test_set=test, method='ridge', returnModel=T)
#cor.test(dts[,2], dts[,3])


## just run fit function from above because i want the model to get coefficients!!!

library(broom)

rr=read.table('../imputation/maize1_reference_ranges.csv', header=T, sep=',')
head(rr)
rr$genic=ifelse(rr$name=='FocusRegion', T, F)


b=broom::tidy(dts)

rr$dts=b$estimate[match(paste0('X', rr$ref_range_id), b$term)]
b=broom::tidy(gy)

rr$gy=b$estimate[match(paste0('X', rr$ref_range_id), b$term)]
b=broom::tidy(gyraw)

rr$gyraw=b$estimate[match(paste0('X', rr$ref_range_id), b$term)]


pdf('~/transfer/ridgeregression_gs.pdf', 14,6)
ggplot(rr, aes(x=range_start, y=estimate)) + geom_point() + facet_wrap(~chrom, ncol=10)
ggplot(rr[rr$ref_range_id!=40979,], aes(x=range_start, y=estimate)) + geom_point() + facet_wrap(~chrom, ncol=10)
dev.off()
sum(rr$estimate>0)
sum(rr$estimate>0, na.rm=T)
sum(rr$estimate<0, na.rm=T)
sum(rr$estimate==0, na.rm=T)
sum(rr$estimate<0, na.rm=T)
sum(rr$estimate[rr$genic]<0, na.rm=T)
sum(rr$estimate[!rr$genic]<0, na.rm=T)
sum(rr$estimate[!rr$genic]>0, na.rm=T)
sum(rr$estimate[rr$genic]>0, na.rm=T)



#run for each different sets of 250 RILs (10 RILs from each family)
foreach(sample = 1:20) %dopar% {
#  sample_NAM_RILs <- fread(paste0(pheno_dir, '/sample_NAM_RILs/sample_NAM_RILs', sample, '.csv'))
#  
#  geno_NAM_all <- geno_NAM_all[geno_NAM_all$taxa %in% sample_NAM_RILs$V2, ]
#  geno_NAM_wona <- geno_NAM_all[,which(unlist(lapply(geno_NAM_all, function(x)!all(is.na(x))))),with=F]
#  
#  common_colnames_NAM <- intersect(colnames(geno_NAMfounders), colnames(geno_NAM_wona))
#  geno_NAM <- rbind(geno_NAMfounders[,..common_colnames_NAM], geno_NAM_wona[, ..common_colnames_NAM])
  
  accuracy_trait <- c()
  for (trait in traits){
    
    pheno_geno_282 <- merge(pheno_282[,c('taxa',trait),with = FALSE], geno_282, by = 'taxa') 
    pheno_geno_NAM <- merge(pheno_NAM[,c('taxa',trait),with = FALSE], geno_NAM, by = 'taxa')
    
    temp_df <- predictPhenotype_crosspopulation(train_set = pheno_geno_NAM, 
                                                test_set = pheno_geno_282, 
                                                method = 'ridge')
    
    write.csv(temp_df, paste0(result_dir, '/predictedvalues_exp_haplotypes_sampleNAMtaxa/', tissue, '_', trait, '_train_NAMsample',sample, '_test_282.csv'))
    
    r <- cor(temp_df[,2], temp_df[,3])
    accuracy_trait <- rbind(accuracy_trait, r)
    message('done with trait ', trait)
    
  }
  rownames(accuracy_trait) <- traits
  
  write.csv(accuracy_trait, paste0(result_dir, '/', tissue, '_prediction_accuracy_train_NAMsample', sample, '_test_282.csv'))
  
  message('done with sample ', sample)
}




pdf('~/transfer/rr_test.pdf')
plot(cv$glmnet.fit, 'lambda')
betas=as.matrix(cv$glmnet.fit$beta)
lambdas=cv$lambda
names(lambdas)=colnames(betas)
library(tidyr)
library(dplyr)
library(ggrepel)
as.data.frame(betas) %>% 
tibble::rownames_to_column("variable") %>% 
pivot_longer(-variable) %>% 
mutate(lambda=lambdas[name]) %>% 
ggplot(aes(x=lambda,y=value,col=variable)) + 
geom_line() + 
geom_label_repel(data=~subset(.x,lambda==min(lambda)),
aes(label=variable),nudge_x=-0.5) +
scale_x_log10()
as.data.frame(betas) %>% 
tibble::rownames_to_column("variable") %>% 
pivot_longer(-variable) %>% 
mutate(lambda=lambdas[name]) %>% 
ggplot(aes(x=lambda,y=value,col=variable)) + 
geom_line() + 
geom_label_repel(data=~subset(.x,lambda==min(lambda)),
aes(label=variable),nudge_x=-0.5)
dev.off()



