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
  
  geno1 = geno$hap1[, iVar]
  geno2 = geno$hap2[, iVar]
  
  heter  = which(geno1 + geno2 == 1)
  
  if(length(heter)>5) ######only applied to analyze sites for which contain at least 5 heterozygotes 
  {
  
    ratio1 <- c(data$y1[heter, iVar])/c(data$r1[heter, iVar])
    ratio1[is.na(ratio1)] <- 1.1
    ratio2 <- c(data$y2[heter, iVar])/c(data$r2[heter, iVar])
    ratio2[is.na(ratio2)] <- 1.1
    idx <- union( which(ratio1 > 1.0), which(ratio2 > 1.0) )
    if (length(idx) > 0) {
      heter = heter[-idx]
    }
    n2 <- length(heter)
    genotypes = c(rbind(geno1[heter], geno2[heter]))
    CountData <-  c(rbind(data$y1[heter, iVar], data$y2[heter, iVar]))
    LibSize <- c(rbind(data$r1[heter, iVar], data$r2[heter, iVar]))
    
    ratio <- CountData/LibSize
    
    numIDV <- length(CountData)
    
    t1 <- system.time( model0 <- glm(formula = ratio~genotypes,
                                     family = binomial(link = "logit"),
                                     weights = LibSize) )
    coef <- coef(summary(model0))
    
    K <- K[heter, heter]
    
    RelatednessMatrix <-kronecker(K, matrix(1, nrow=2, ncol=2))
    
    eig               <- eigen(RelatednessMatrix)
    eigval            <- eig$value
    eigvector         <- eig$vectors
    if(any(eigval<1e-10)){ 
      RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T)$mat)	
    }
    rm(eig)
    rm(eigval)
    rm(eigvector)
    RelatednessMatrix <- list( RelatednessMatrix, diag(rep(1, 2*n2)))
    
    numIDV <- beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA
    numIDV=length(LibSize)
    model0$numTotal <- LibSize
    model0$numSucc  <- CountData
    tmpRelatednessMatrix=RelatednessMatrix
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
                                converged = converged,time=t1[3]))
  }
  
  
  
}

saveRDS(summary, paste0('IMAGE-A.rds'))