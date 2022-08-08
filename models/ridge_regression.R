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


teh=read.table('geno_pheno_gphenos.2022-08-05.txt', header=T)


## set up 
## te focus
## all repeat focus
## # columns got to be: 1=sample name, 2=pheno value, 3:end=genotypes
tefocus=teh[,c(1,27:30,6:16)]
repeatfocus=teh[,c(1,27:30,6:16, 18:22)]



### nam sampling approach from anju
## https://bitbucket.org/bucklerlab/haplotype_prediction/src/master/HARE_redone/src/4c_phenotype_prediction_outofpopulation_NAMsample_sampletraits.R

#function for phenotype prediction

predictPhenotype <- function(train_set, ## columns got to be: 1=sample name, 2=pheno value, 3:end=genotypes
                                             test_set,
                                             method = c("ridge", "LASSO"),
                                             nFolds = 10,
                                             ...) {
  method <- method[1]
  
  if (method == "ridge") {
    message("Running ridge-regression...")
    
#   #select common x in train and test set
#   common_x <- intersect(colnames(train_set), colnames(test_set))
#   train <- train_set[, ..common_x] %>% na.omit()
#   test <- test_set[, ..common_x] %>% na.omit()
    
    # transform both sets into matrices using `data.matrix`
    x_train_mat <- data.matrix(train[, 3:ncol(train)])
    x_test_mat  <- data.matrix(test[, 3:ncol(test)])
    
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
   return(
     data.frame(
       name = test[,1],
       predicted = pred,
       observed = test[,2]
     )
   )
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




train=tefocus[r,-c(2:3,5)] %>% na.omit()
dim(train)
test=tefocus[-r,-c(2:3,5)] %>% na.omit()
temp=predictPhenotype(train_set=train, test_set=test, method='ridge')
cor.test(temp[,2], temp[,3])



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

   # transform both sets into matrices using `data.matrix`
    x_train_mat <- data.matrix(train[, 3:ncol(train)])
    x_test_mat  <- data.matrix(test[, 3:ncol(test)])
    
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



