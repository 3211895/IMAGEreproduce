library('Rcpp')
library('Matrix')
sourceCpp('AI.cpp')

source('simulation_function.R')

load(file)  #######load the simulation data

n=nrow(geno$hap1)  ###sample size
m=ncol(geno$hap1)  #####number of SNP-CpG pairs

summary<-data.frame()

for(iVar in 1:m)
{
  cat(iVar, '\n')
  
  genotypes = c(geno$hap1[, iVar] + geno$hap2[, iVar])
  genotypes=genotypes/2
  
  RelatednessMatrix <- K
  CountData <- c(data$y[, iVar])
  LibSize <- c(data$r[, iVar])
  ratio <- CountData/LibSize
  ratio[is.na(ratio)] <- 1.1
  idx <- which(ratio > 1.0)
  if (length(idx) > 0){
    CountData <- CountData[-idx]
    LibSize   <- LibSize[-idx]
    Phenotypes <- Phenotypes[-idx]
    ratio <- ratio[-idx]
    RelatednessMatrix <- RelatednessMatrix[-idx, -idx]
  }
  numIDV <- length(CountData)
  t1 <- system.time( model0 <- glm(formula = ratio~genotypes,
                                   family = binomial(link = "logit"),
                                   weights = LibSize) )
  coef <- coef(summary(model0))
  
  eig               <- eigen(RelatednessMatrix)
  eigval            <- eig$value
  eigvector         <- eig$vectors
  if(any(eigval<1e-10)){ 
    #   warning("PQLseq::the relatedness matrix is singular, it has been modified!")
    RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T)$mat)	
  }
  rm(eig)
  rm(eigval)
  rm(eigvector)
  numIDV <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
  numIDV=length(CountData)
  model0$numTotal <- LibSize
  model0$numSucc  <- CountData
  tmpRelatednessMatrix <- list(RelatednessMatrix, diag(nrow=numIDV))
  names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
  t1 <- system.time(model1 <- try( PQLseq.fit(model0, tmpRelatednessMatrix) ))
  if(class(model1) != "try-error"){
    # numIDV <- length(idx)
    beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
    se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
    pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
    sigma2      <- model1$theta[2]+model1$theta[3]
    h2          <- model1$theta[2]/(sigma2)
    tau1        <- model1$theta[2]
    tau2        <- model1$theta[3]
    converged   <- model1$converged
  }else{converged <- FALSE}
  summary <- rbind(summary,data.frame(Loc=iVar,numIDV = numIDV, beta = beta, se_beta = se_beta, 
                              pvalue = pvalue, h2 = h2, sigma2 = sigma2, 
                              converged = converged,time=t1[3],
                              Method='Add_pql_g'))

}

saveRDS(summary, paste0('IMAGE-I.rds'))