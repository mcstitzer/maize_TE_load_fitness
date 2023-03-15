

## from Guillaume
## https://bitbucket.org/bucklerlab/ames_nam_hybrid/src/master/src/8d-R1-GY_adjustments.R


# Traits whose BLUEs to estimate
selected_traits <- c("GY")
covariates <- list(c("DTS"))

# Output
trait.summary_file <- "Ames_NAM_hybrid/data/GY_adjustments.csv"



##############################################################
# Functions
##############################################################
library(asreml)
library(ggplot2)

## had to change asreml.control() to asreml.options() for shift from v3 to v4
asreml.options(maxiter = 100, workspace = 2048e+06)


getBLUPs <- function (pheno.data, traits, covariates = NULL, pdf.out = NULL,
                      s1.model = "ar1v(row):ar1(range)", s2.model = NULL) {
  
  require(asreml)
  
  # Updating the levels at factors
  pheno.data$env <- factor(pheno.data$env)
  pheno.data$genotype <- factor(pheno.data$genotype)
  
  if ("rep" %in% colnames(pheno.data)) pheno.data$rep <- factor(pheno.data$rep)
  if ("field" %in% colnames(pheno.data)) pheno.data$field <- factor(pheno.data$field)
  if ("set" %in% colnames(pheno.data)) pheno.data$set <- factor(pheno.data$set)
  
  # Considering only the traits with some non-missing values
  pheno.data <- pheno.data[, apply( pheno.data, 2, function(x) any(!is.na(x)) )]
  
  traits <- intersect( traits, names(pheno.data) )
  
  # Formatting the list of covariates by trait
  if (is.null(covariates)) covariates <- lapply(traits, function(trait) NULL)
  
  if (length(traits) != length(covariates)) stop("'traits' and 'covariates' must be of same length")
  
  names(covariates) <- traits
  
  # Output graphics?
  if (!is.null(pdf.out)) pdf(pdf.out)
  
  # BLUP file
  trait_BLUPs <- data.frame(genotype=as.character(unique(pheno.data$genotype)),
                            stringsAsFactors=FALSE)
  
  # r2 file
  trait_r2 <- data.frame(genotype=as.character(unique(pheno.data$genotype)),
                         stringsAsFactors=FALSE)
  
  # Trait summary
  trait_summary <- data.frame()
  
  # Model diagnostics
  model_summary <- data.frame(trait = traits, sigma2 = NA, loglik = NA, BIC = NA, model_type = NA, converge = NA)
  
  # Merging genotype BLUPs to the BLUP file
  for ( trait in traits ) {
    
    cat(paste("\n\t\t\t", trait,":\n",sep=""))
    
    if (is.null(covariates[[trait]])) {
      to_keep <- is.finite(pheno.data[, trait])
    } else {
      to_keep <- is.finite(pheno.data[, trait]) & is.finite(pheno.data[, covariates[[trait]]])
    }
    
    # Table of panels by env
    tab <- table(pheno.data[to_keep,"env"], pheno.data[to_keep,"field"])
    print(tab)
    cat("\n\n")
    
    # Number of envs with data for the trait
    env_levels <- unique(pheno.data[to_keep, "env"])
    
    # Fixed-effect formula
    fixed <- paste(trait, "~ 1")
    fixed <- ifelse(is.null(covariates[[trait]]), fixed, paste(fixed, paste(covariates[[trait]], collapse="+"), sep=" + "))
    
    # Random-effect formula
    if (length(env_levels) > 1) {
      random_0 <- "~ genotype + env + env:field"
      random_s1 <- ifelse( !is.null(s1.model), paste( random_0, s1.model, sep = " + "), random_0 )
      random_s2 <- ifelse( !is.null(s2.model), paste( random_0, s2.model, sep = " + "), random_0 )
    } else {
      random_0 <- "~ genotype + field"
      random_s1 <- ifelse( !is.null(s1.model), paste( random_0, s1.model, sep = " + "), random_0 )
      random_s2 <- ifelse( !is.null(s2.model), paste( random_0, s2.model, sep = " + "), random_0 )
    }
    
    # Model fitting
    cat(paste("Model:\n\tFixed:", fixed, "\n\tRandom:", random_s1,"\n"))
    
    fit <- asreml(fixed = as.formula(fixed),
                  random = as.formula(random_s1),
                  na.action= na.method(x = "omit", y="omit"),
  #                na.method.X = "omit",
  #                na.method.Y = "omit",
                  data = pheno.data#,
          #        control = asreml.options(maxiter = 100, workspace = 2048e+06)
    )
    
    model_type <- "S1"
    
    if (!fit$converge) {
      
      cat("\n\t----> Fitting alternate spatial model:\n")
      cat(paste("Model:\n\tFixed:", fixed, "\n\tRandom:", random_s2,"\n"))
      
      fit <- asreml(fixed = as.formula(fixed),
                    random = as.formula(random_s2),
                    na.action=na.method(x = "omit", y="omit"),
    #                na.method.X = "omit",
    #                na.method.Y = "omit",
                    data = pheno.data)#,
              #      control = asreml.options(maxiter = 100, workspace = 2048e+06))
      
      model_type <- "S2"
      
    }
    
    if (!fit$converge) {
      
      cat("\n\tFitting simple model:\n")
      cat(paste("Model:\n\tFixed:", fixed, "\n\tRandom:", random_0,"\n"))
      
      fit <- asreml(fixed = as.formula(fixed),
                    random = as.formula(random_0),
                    na.action=na.method(x = "omit", y="omit"),
     #               na.method.X = "omit",
     #               na.method.Y = "omit",
                    data = pheno.data)#,
                  #  control = asreml.options(maxiter = 100, workspace = 2048e+06))
      
      model_type <- "simple"
      
    }
    
    # Summary of fitting process
    N <- sum( !is.na(fit$residuals) )
    BIC <- -2*fit$loglik + fit$nwv*log(N)
    
    model_summary[model_summary$trait == trait, -1] <- c(fit$sigma2, fit$loglik, BIC, model_type, fit$converge)
    
    # Extracting BLUPs
    trait_BLUPs[, trait] <- sapply(as.character(trait_BLUPs$genotype), function(x) {
      
      effect_id <- paste("genotype", x, sep="_")
      
      if ( (effect_id %in% rownames(coef(fit)$random)) & (x %in% unique(pheno.data[to_keep, "genotype"])) ) {
        out <- coef(fit)$random[effect_id, 1]
      } else { 
        out <- NA
      }
      
      return(out)
      
    })
    
    # Variance and heritability
    V_G <- fit$gammas["genotype!genotype.var"] * fit$sigma2
    V_R <- fit$gammas["R!variance"] * fit$sigma2
    
    V_P <- V_G + V_R
    
    H2 <- V_G / V_P
    
    # Reliabilities
    pev <- fit$vcoeff$random * fit$sigma2
    names(pev) <- rownames(coef(fit)$random)
    
    r2 <- sapply(as.character(trait_r2$genotype), function(x) {
      
      effect_id <- paste("genotype", x, sep="_")
      
      if ( (effect_id %in% rownames(coef(fit)$random)) & (x %in% unique(pheno.data[!is.na(pheno.data[, trait]), "genotype"])) ) {
        pevG <- pev[effect_id]
        out <- 1 - pevG/V_G
      } else {
        out <- NA
      }
      
      return(out)
      
    })
    
    trait_r2[, trait] <- r2
    
    # Fixed effects
    fixed_effects <- coef(fit)$fixed[, "effect", drop=FALSE]
    
    # Summary
    trait_summary <- rbind(trait_summary,
                           cbind(
                             data.frame(trait = trait,
                                        envs = paste(unique(pheno[to_keep, "env_name"]), collapse = ", "),
                                        V_P = V_P,
                                        V_G = V_G,
                                        H2 = H2,
                                        r2.mean = mean(r2, na.rm=TRUE),
                                        r2.sd = sd(r2, na.rm=TRUE),
                                        r2.min = min(r2, na.rm=TRUE),
                                        r2.max = max(r2, na.rm=TRUE)
                             ),
                             data.frame(t(fixed_effects), check.names=FALSE)
                           ))
    
    # Diagnostic plots
    fitted <- fitted(fit)[to_keep]
    resid <- pheno.data[to_keep, trait] - fitted
    
    if (!is.null(pdf.out)) {
      
      data_e <- data.frame(y_hat = fitted,
                           e = scale(resid),
                           env = pheno.data[to_keep, "env"]
      )
      
      data_u <- data.frame(u = scale(trait_BLUPs[,trait]))
      
      plot.new()
      text(x=0.5, y=0.5, trait)
      
      # Residual plot
      plot(
        ggplot(data = data_e, mapping = aes(x = y_hat, y = e)) + 
          geom_point() + stat_smooth(method = loess) + geom_abline(intercept = 0, slope = 0, linetype = "longdash") +
          labs(x = "Predicted values", y = "Standardized residuals", title = paste(trait,"- Residual plot")) + theme_bw(base_size = 16)
      )
      
      # Normal Q-Q plot
      plot(
        ggplot(data = data_e, mapping = aes(sample = e)) + 
          stat_qq() + 
          geom_abline(intercept = 0, slope = 1, linetype = "longdash") +
          labs(x = "Expected quantiles", y = "Observed quantiles (standardized)", title = paste(trait,"- Normal Q-Q plot of residuals")) + theme_bw(base_size = 16)
      )
      
      # Histogram of BLUPs
      range_u <- max(data_u$u, na.rm=TRUE) - min(data_u$u, na.rm=TRUE)
      
      plot(
        ggplot(data = data_u, mapping = aes(x = u)) + 
          geom_histogram(aes(y = ..density..), fill="lightblue", colour = "blue", binwidth = range_u/25) + geom_density(color = "darkblue") +
          labs(x = "Standardized genotype BLUPs", y = "Density", title = paste(trait,"- Histogram of BLUPs")) + theme_bw(base_size = 16)
      )
      
    }
    
  }
  
  if (!is.null(pdf.out)) graphics.off()
  
  # Returning the output files
  output <- list(BLUP = trait_BLUPs,
                 r2 = trait_r2,
                 trait.summary = trait_summary,
                 model.summary = model_summary
  )
  
  return(output)
  
}


##############################################################
# NAM
##############################################################
#-------------------------------------------------------------
# Raw data
#-------------------------------------------------------------
# Phenotypic file
pheno_file <- "../phenotypes/NamTC_asreml160422.csv"

# Input
pheno <- read.csv(pheno_file)
colnames(pheno)[colnames(pheno) == "Yield"] <- "GY"

# Formatting
for (variable in c("Order", "field", "pblock", "pop", "entrynum", "env")) {
  pheno[, variable] <- as.factor(pheno[, variable])
}
pheno$row_factor <- as.factor(pheno$row)
pheno$range_factor <- as.factor(pheno$range)

pheno$env_field <- as.factor(pheno$env):as.factor(pheno$field)

# Genotype names
pheno$genotype <- as.character(pheno$female_pedigree)
pheno$genotype <- paste(pheno$genotype, "PHZ51", sep="/")

# Location variable
pheno$loc <- as.factor(substr(as.character(pheno$env_name), 3, 4))

#-------------------------------------------------------------
# Description
#-------------------------------------------------------------
n_traits <- length(selected_traits)

# Structure by env
n_env <- numeric(length(selected_traits))
names(n_env) <- selected_traits

for (trait in selected_traits) {
  
  cat(paste(sep="", "--------------\n", trait, "\n--------------\n"))
  
  tmp <- pheno[!is.na(pheno[,trait]), ]
  
  print( tab <- table(tmp$field, tmp$env_name) )
  
  n_env[trait] <- ncol(tab)
  
}

#-------------------------------------------------------------
# Estimating BLUPs
#-------------------------------------------------------------
MMA.BLUP <- getBLUPs(pheno.data=pheno,
                     traits=selected_traits,
                     covariates=covariates,
                     s1.model="idv(env_field):ar1(row_factor):ar1(range_factor)")

write.table(MMA.BLUP$trait.summary, trait.summary_file, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
