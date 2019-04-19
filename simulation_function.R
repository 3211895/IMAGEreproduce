
#' @export
PQLseq.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {
  
  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  # if((method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
  if(method.optim == "AI") {
    fixtau.old <- rep(0, length(RelatednessMatrix)+1)
    # to use average information method to fit alternative model
    model1 <- PQLseq.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new <- 1*(model1$theta < 1.01 * tol)
    
    while(any(fixtau.new != fixtau.old)) {
      fixtau.old <- fixtau.new
      # to use average information method to fit alternative model
      model1 <- PQLseq.AI(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(model1$theta < 1.01 * tol)
    }
    return(model1)
  }
}



##########################################################
#       PQLseq FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################
#' @export
PQLseq.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {
  
  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
  numIDV <- length(y)
  offset <- model0$offset
  if(is.null(offset)) {offset <- rep(0, numIDV)}
  
  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  D <- mu.eta/sqrt(model0$family$variance(mu))
  
  if(family$family %in% c("binomial")){
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
    mu <- model0$numTotal*mu
  }
  
  Y <- eta - offset + (y - mu)/mu.eta	
  X <- model.matrix(model0)
  alpha <- model0$coef
  
  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  numK <- length(RelatednessMatrix)
  idxtau <- which(fixtau == 0)
  numK2 <- sum(fixtau == 0)
  if(numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)  
    
    H <- tau[1]*diag(1/D^2)
    for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}
    
    Hinv <- chol2inv(chol(H))
    HinvX <- crossprod(Hinv, X)
    XHinvX <- crossprod(X, HinvX)
    
    P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
    
    if(class(P) == "try-error"){
      stop("Error in P matrix calculation!")
    }
    
    PY <- crossprod(P, Y)
    tau0 <- tau
    for(ik in 1:numK2) {
      if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
      else {
        PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
        tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
      }
    }
  } 
  
  for (iter in seq_len(maxiter)) {	
    alpha0 <- alpha
    tau0 <- tau
    model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)
    
    tau <- as.numeric(model1$tau)
    cov <- as.matrix(model1$cov)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta) + offset
    
    
    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))
    
    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
      mu <- model0$numTotal*mu
    }
    
    Y <- eta - offset + (y - mu)/mu.eta
    
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
    if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {
      
      iter <- maxiter
      break
    }
  }
  
  converged <- ifelse(iter < maxiter, TRUE, FALSE) 

  
  res <- y - mu
  P <- model1$P
  return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, P = P, residuals = res, cov = cov, converged = converged))
}# end function

# fit beta-binomial models
# alternative
calc_logl_alt<-function(param) {
  phi=exp(param[1])/(exp(param[1])+1)
  mu=param[2]
  gamma=param[3]
  theta=exp(c*gamma+mu)/(exp(c*gamma+mu)+1)
  
  alpha=theta*(1-phi)/phi
  beta=(1-theta)*(1-phi)/phi
  logl=sum(lbeta(alpha+x, beta+y-x))-sum(lbeta(alpha, beta))
  return(logl)
}

# null
calc_logl_null<-function(param) {
  phi=exp(param[1])/(exp(param[1])+1)
  mu=param[2]
  gamma=0
  theta=exp(mu)/(exp(mu)+1)
  
  alpha=theta*(1-phi)/phi
  beta=(1-theta)*(1-phi)/phi
  
  logl=sum(lbeta(alpha+x, beta+y-x))-length(x)*lbeta(alpha, beta)
  return(logl)
}

# fit the beta-binomial for each site, where mcounts is a vector of methylated counts, tcounts is a vector of total counts and cvt is a vector of covariates
bbfit<-function (mcounts, tcounts, cvt) {
  x=mcounts; y=tcounts; c=cvt
  param=c(-1,0.5,0)
  phi=exp(param[1])/(exp(param[1])+1)
  mu=param[2]
  gamma=param[3]
  theta=exp(c*gamma+mu)/(exp(c*gamma+mu)+1)
  
  alpha=theta*(1-phi)/phi
  beta=(1-theta)*(1-phi)/phi
  
  phi=1/(alpha+beta+1); mu=log(alpha/beta); gamma=0;
  par=c(-2, mu)
  
  for (i_iter in 1:100) {
    if (i_iter==1) {
      null=optim(par,calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
    } else {
      null=optim(alt$par[1:2],calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
    }
    alt=optim(c(null$par,0),calc_logl_alt, control=list(fnscale=-1),  method = c("CG"))
    
    while (null$value==0 | alt$value==0 | alt$value<(null$value-0.01)) {
      alpha=runif(1); beta=runif(1)
      phi=1/(alpha+beta+1); mu=log(alpha/beta); gamma=0;
      par=c(log(phi), mu)
      
      null=optim(par,calc_logl_null, control=list(fnscale=-1),  method = c("CG"))
      alt=optim(c(null$par,0),calc_logl_alt, control=list(fnscale=-1),  method = c("CG"))
    }
    
    lr=2*(alt$value-null$value)
    if (lr<0) {lr=0}
    
    if (i_iter==1) {
      lr_old=lr
      lr_new=lr
    } else {
      lr_new=lr
      if (abs(lr_new-lr_old)<0.01 & lr!=0) {break}
      lr_old=lr
    }
  }
  
  return(1-pchisq(lr,df=1))
}
