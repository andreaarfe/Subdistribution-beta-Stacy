
# setup R session
library(mvtnorm)   # for multivariate normal
library(compiler)  # to compile the functions
library(MCMCpack)  # for the inverse-gamma distribution

#================================================================================#
# Goal: Sample small-shape gamma random variables via accept-reject              #
# Paper: arXiv:1302.1884                                                         #
# Function: rgamss                                                               #
# Authors: C. Liu, R. Martin (www.math.uic.edu/~rgmartin), and N. Syring         #
# Version: 2nd (04/07/2015)                                                      #
#================================================================================#

# INPUT
# n = sample size
# shape = shape parameter (same for all samples)
# scale = scale parameter (same for all samples; default = 1)
# do.log = logical -- return values on log scale

# OUTPUT
# vector of length n containing gamma samples;
# on log-scale depending on do.log or magnitude of shape

rgamss <- function(n, shape, scale=1, do.log=TRUE) {
  a <- shape
  if(a > 0.2){
    oo <- rgamma(n, a, 1 / scale)
    oo <- ifelse(do.log,log(oo),oo)
  } else {
    e1 <- 2.71828182845905
    L <- 1 / a - 1
    w <- a / e1 / (1 - a)
    ww <- 1 / (1 + w)
    eta <- function(z) if(z >= 0) exp(-z) else w * L * exp(L * z)
    h <- function(z) exp(-z - exp(-z / a))
    rh <- function(a) {
      repeat {
        U <- runif(1)
        z <- if(U <= ww) -log(U / ww) else log(runif(1)) / L
        if(h(z) / eta(z) > runif(1)) return(z)
      }
    }
    Z <- numeric(n)
    for(i in 1:n) Z[i] <- rh(a)
    o <- log(scale) - Z / a
    if(!do.log) {
      oo <- exp(o)
      if(any(oo == 0)) {
        oo <- o
        warning("Output given on log-scale since shape is small")
      }
    } else oo <- o
  }
  return(oo)
}


####################################################################################
# Function to simulate from the dirichlet distribution using the rgamss function
# to allow for small parameters values.
####################################################################################

myrdirichlet <- function(n,alpha){
  x <- matrix(0,ncol=length(alpha),nrow=n)
  for(i in 1:n){
    x[i,] <- sapply(1:length(alpha),function(i) rgamss(1,alpha[i]))
    x[i,] <- x[i,] - max(x[i,]) # to avoid Inf values in exp(x[i,])
    x[i,] <- exp(x[i,])/sum(exp(x[i,]))
  }
  return(x)
}

#####################################################################################
# Definition of the functions required for the analyses
#####################################################################################

# Function computing the logarithm of the centering subdistribution function F0 and the precision
# parameters omega for given values of the hyperparameters. For F0, it returns a nxk matrix
# with all possible values of F0. For omega, it return a n-vector with the values of the
# precision weight corresponding to each specified time point.
# INPUTS: theta   = p(2k-1)+k vec of parameters (b1,b2,...,bk-1,v1,...,vk,u1,...,uk)
#         times   = n vec of event times at which to compute F0
#         W       = nxp matrix of regressors
#         k       = number of event types
#         log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
#         dist    = centering distribution of event times; either "weibull" or "lognormal"
#         omega_m = scalar used to control the prior precision; divides the values of omega
getF0omega <- function(theta=NULL,
                       times=NULL,
                       delta=1,
                       W=NULL,
                       k=NULL,
                       omega_m=1,
                       log_u=FALSE,
                       dist="weibull"){
  p <- ncol(W)   # number of regressors
  # scomposition of theta 
  if(dist=="weibull"){
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p) # multinomial logistic regression parameters
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p) # weibull regression parameters
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)] # weibull shape
    if(log_u) u <- exp(u)
  }
  if(dist=="lognormal"){
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p) # multinomial logistic regression parameters
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p) # lognormal regression parameters
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)] # lognormal scale parameter
    if(log_u) u <- exp(u)
  }
  
  # multinomial logistic model for the type of event, on the log-scale
  # F01 is a nxk matrix
  F01      <- cbind(W%*%b,0)
  F01.norm <- rowSums(exp(F01))
  F01      <- F01-log(F01.norm)
  
  # Evaluates the cumulative probabilities by suppressing R warnings.
  # These may be generated if a location or scale parameter <=0 is provided.
  # Be careful with the inputs!!!
  # Works on the log-scale. F02 and F02delta are nxk matrices.
  if(dist=="weibull"){
    F02      <- suppressWarnings(pweibull(times,matrix(u,ncol=length(u),nrow=nrow(W),byrow=TRUE),exp(t(t(-W%*%v)/u)),log=TRUE)) 
    F02lag   <- suppressWarnings(pweibull(times-delta,matrix(u,ncol=length(u),nrow=nrow(W),byrow=TRUE),exp(t(t(-W%*%v)/u)),log=TRUE))
  }
  if(dist=="lognormal"){
    F02      <- suppressWarnings(plnorm(times,W%*%v,matrix(u,ncol=length(u),nrow=nrow(W),byrow=TRUE),log=TRUE)) 
    F02lag   <- suppressWarnings(plnorm(times-delta,W%*%v,matrix(u,ncol=length(u),nrow=nrow(W),byrow=TRUE),log=TRUE))
  }
  F02delta <- log(exp(F02) - exp(F02lag))
  
  # Final mean subdistribution function and increments, expressed on the log-scale
  F0      <- F01 + F02 
  F0delta <- F01 + F02delta
  
  # Precision parameters
  omega <- delta/rowSums(exp(F0delta))
  omega[is.infinite(omega)] <- 1e6 # threshold to avoid overflows
  omega <- omega/omega_m # applies the multiplier omega_m
  
  return(list(
    omega     = omega,    # precision parameters
    F0delta   = F0delta,  # log-increments of the mean subdistribution function
    F0        = F0        # log-mean subdistribution function
  ))
}

# function computing the log-posterior distribution for the parametric models
# INPUTS: theta   = p(2k-1)+k vec of parameters (b1,b2,...,bk-1,v1,...,vk,u1,...,uk)
#         times   = n vec of event times
#         types   = n vec of event types (0 for censored)
#         k       = n of event types
#         W       = nxp matrix of regressors
#         b0      = px(k-1) matrix of prior means of b
#         Sb      = pxpx(k-1) array of prior cov of b
#         v0      = pxk matrix of prior means of v
#         Sv      = pxpxk array of prior cov of v
#         pu      = pxk vec of hyper alphas for u
#         qu      = pxk vec of hyper betas for u
#         log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
#         verbose = if FALSE (default) only return the value of the log-posterior distribution (the kernel of),
#                   if TRUE also returns the log-likelihood and log-hyperpriors in a list
#         dist    = centering distribution of event times; either "weibull" or "lognormal"
logposterior_param <- function(theta=NULL,times=NULL,delta=1,types=NULL,k=NULL,W=NULL,b0=NULL,Sb=NULL,
                        v0=NULL,Sv=NULL,pu=NULL,qu=NULL,verbose=FALSE,log_u=FALSE,
                        dist="weibull"){
  p <- ncol(W)   # number of regressors
  
  # logprior
  if(dist=="weibull"){
    # scomposition of theta in (b,v,u)
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p)
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p)
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)]
    if(log_u){ u <- exp(u) }
    # log prior for b
    logpriorb <- sum(sapply(1:(k-1),function(i)
      mvtnorm::dmvnorm(b[,i],b0[,i],as.matrix(Sb[,,i]),log=TRUE)))
    # log prior for v
    logpriorv <- sum(sapply(1:k,function(i)
      mvtnorm::dmvnorm(v[,i],v0[,i],as.matrix(Sv[,,i]),log=TRUE)))
    # log prior for u
    logprioru <- sum(dgamma(u,pu,qu,log=TRUE))
    # total weibull log prior
    logprior <- logpriorb + logpriorv + logprioru
  }
  if(dist=="lognormal"){
    # scomposition of theta in (b,v,u)
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p)
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p)
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)]
    if(log_u){ u <- exp(u) }
    # log prior for b
    logpriorb <- sum(sapply(1:(k-1),function(i)
      mvtnorm::dmvnorm(b[,i],b0[,i],as.matrix(Sb[,,i]),log=TRUE)))
    # log prior for v
    logpriorv <- sum(sapply(1:k,function(i)
      mvtnorm::dmvnorm(v[,i],v0[,i],as.matrix(Sv[,,i]),log=TRUE)))
    # log prior for u
    dinvg <- sapply(1:length(u), function(i) dinvgamma(u[i]^2,pu[i],qu[i]))
    logprioru <- sum(log(dinvg)+log(u))
    # total lognormal log prior
    logprior <- logpriorb + logpriorv + logprioru
  }
  
  # loglikelihood
  z <- types!=0  # no censorship indicator
  n <- nrow(W)   # sample size
  F0omega <- getF0omega(theta=theta,W=W,k=k,times=times,
                        delta=delta,log_u=log_u,dist=dist)
  F0      <- F0omega$F0 
  F0delta <- F0omega$F0delta
  F0delta.full <- F0delta
  for(i in 1:n){
    delete <- ((1:k)!=types[i])
    F0delta[i,delete] = 0
  }
  F0delta <- rowSums(F0delta)
  cumF0  <- rowSums(exp(F0))
  # computes the loglikelihood; warnings about NaN being 
  # generated by the log are suppressed; a NaN may be 
  # generate by the log if it receives a negative argument;
  # this may happen for numerical approximation when a 
  # cumF0 has a value very close to 1.
  loglik <- suppressWarnings(sum(z*F0delta + (1-z)*log(1-cumF0))) 
  if(is.na(loglik)){ loglik <- -Inf } # to ensure that illicit values are rejected
  
  # log posterior
  logpost <- loglik + logprior
  if(log_u){ logpost <- logpost + sum(log(u)) } #adds the Jacobian when considering with log-transformed u's
  
  # output
  if(verbose==FALSE) return(logpost = logpost)  
  if(verbose==TRUE) return(list(
    logpost   = logpost,
    loglik    = loglik,
    logpriorb = logpriorb,
    logpriorv = logpriorv,
    logprioru = logprioru
  ))
}

# Function to compute the conditional posterior expected cause-specific cumulative hazard
# functions (A), overall survival function (S), precision parameters (omega.post)
# and subdistribution function (F.post) over the time interval from 0 to tau.
# Estimates are conditional on the observed times and event types and the given value 
# of the covariates and hyperparameters
# INPUTS
# tau     = (scalar) upper bound for the event times
# times   = vector of observed event times
# k       = n of event types
# Wpred   = 1xp matrix of regressors of the target subdistribution function
# Wobs    = nxp matrix of observed regressors
# log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
# dist    = centering distribution of event times; either "weibull" or "lognormal"
# omega_m = scalar used to control the prior precision; divides the values of omega
postest <- function(tau=NULL,
                    times=NULL,
                    delta=1,
                    types=NULL,
                    k=NULL,
                    Wpred=NULL,
                    Wobs=NULL,
                    theta=NULL,
                    omega_m=1,
                    log_u=FALSE,
                    dist="weibull"){
  p <- ncol(Wpred) #number of predictors
  times_est <- seq(0,tau,by=delta)
  L <- length(times_est)
  # for safety, checks that Wobs and Wpred have the same number of colums
  if(ncol(Wobs)!=p) stop("Wpred and Wobs must have the same number of columns.") 
  # expands Wpred to be L x p
  Wpred <- matrix(Wpred,ncol=p,nrow=L,byrow = TRUE)
  F0omega    <- getF0omega(theta=theta,times=times_est,W=Wpred,k=k,
                           delta=delta,log_u=log_u,dist=dist,omega_m=omega_m)
  F0         <- exp(F0omega$F0)
  F0delta    <- exp(F0omega$F0delta)
  omega      <- F0omega$omega
  cumF0      <- rowSums(F0)
  cumF0delta <- rowSums(F0delta)
  cumF0lag   <- c(0,cumF0[1:(length(cumF0)-1)])
  # selects only the observations with regressors equal to Wpred
  I <- which(apply(Wobs,1,function(r){ all(r==Wpred[1,]) }))
  times <- times[I]
  types <- types[I]
  # computes the posterior expectation
  A      <- matrix(0,nrow=L,ncol=k)
  S      <- matrix(0,nrow=L,ncol=1)
  F.post <- matrix(0,nrow=L,ncol=k)
  l      <- rep(0,times=L)
  m      <- matrix(0,ncol=(k+1),nrow=L)
  den    <- matrix(0,nrow=L,ncol=1)
  S_ <- 1
  for(tt in 1:L){
    m[tt,] <- sapply(0:k,function(d) sum((times==times_est[tt])&(types==d)) )
    l[tt] <- sum(times>times_est[tt])
    # accounts for the fact that 1-cumF0lag may be negative 
    # due to numerical errors if cumF0lag is very close to 1
    # Also accounts that omega may be Inf for when one tries to evaluate the 
    # posterior estimate of F at times<=0
    numS <- ifelse(is.infinite(omega[tt]) & cumF0delta[tt]==0,0,omega[tt]*cumF0delta[tt]) + sum(m[tt,-1])
    temp <- 1-cumF0lag[tt]; temp[temp<0] <- 0
    den[tt]  <- ifelse(is.infinite(omega[tt]) & temp==0,0,omega[tt]*temp) + l[tt] + sum(m[tt,]) 
    # Accounts for the fact that den[tt] may be 0 and that 
    # 1-numS/den[tt] may be negative for large tt due to numerical approximations
    S_ <- S_ * ifelse(den[tt]>0,max(0,(1 - numS/den[tt])),0)
    S[tt] <- S_
  }
  Slag   <- c(1,S[1:length(S)-1])
  for(C in 1:k){
    A_ <- 0
    for(tt in 1:L){
      # Accounts that omega may be Inf for when one tries to evaluate the 
      # posterior estimate of F at times<=0
      numA <- ifelse(is.infinite(omega[tt]) & F0delta[tt,C]==0,0,omega[tt]*F0delta[tt,C]) + m[tt,C+1]
      # Accounts for the fact that den[tt] may be 0 for large tt
      A_ <- A_ + ifelse(den[tt]>0,numA/den[tt],0)
      A[tt,C] <- A_
    }
    Adelta <- c(A[1,C],diff(A[,C]))
    F.post[,C] <- cumsum(Slag*Adelta)
  }
  # Computes the posterior precision parameters.
  # accounts for the fact that 1-rowSums(F.post) or 1-rowSums(F.post) may be negative 
  # due to numerical errors if rowSums(F.post) or rowSums(F.post) are very close to 1
  # Additionally, performs a "safe division" to avoid Inf values when 1-rowSums(F.post)
  # is close to 0.
  f1 <- 1-rowSums(F0);     f1[f1<0] <- 0
  f2 <- 1-rowSums(F.post); f2[f2<0] <- 0
  omega.post <- (omega*f1+l+m[,1])/f2
  omega.post[is.infinite(omega.post)] <- 1e6
  
  # output
  return(list(A=A,S=S,F.post=F.post,omega.post=omega.post))
}

# function computing the marginal likelihood for the subdistribution-beta-Stacy model
# INPUTS: 
#         times   = n vec of event times
#         types   = n vec of event types (0 for censored)
#         k       = n of event types
#         W       = nxp matrix of regressors
#         tau     = 
#         delta   =
#         theta   = p(2k-1)+k vec of parameters (b1,b2,...,bk-1,v1,...,vk,u1,...,uk)
#         log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
#         dist    = centering distribution of event times; either "weibull" or "lognormal"
#         omega_m = scalar used to control the prior precision; divides the values of omega
loglik <- function(times=NULL,
                   types=NULL,
                   k=NULL,
                   W=NULL,
                   tau=NULL,
                   delta=1,
                   theta=NULL,
                   omega_m=1,
                   log_u=TRUE,
                   dist="weibull"){
  
  p        = ncol(W)              # number of predictors
  n        = nrow(W)              # num of observations
  W.unique = as.matrix(unique(W)) # unique values of predictors
  m        = nrow(W.unique)       # num of unique predictors
  tt       = seq(delta,tau,delta) # discrete time points
  
  # find how many reps of each config of W
  ni <- sapply(1:m, function(i) sum(apply(W,1,function(x) all(x == W.unique[i,]))) )
  
  # Computes the marginal log-likelihood; uses the formulas in Corollary 4.1 to 
  # compute the predictive distributions.
  loglik <- 0
  for(i in 1:m){
    # select obs with equal regressors
    I       = apply(W,1,function(x)all(x == W.unique[i,]))
    times.i = times[I]
    types.i = types[I]
    ni.i    = ni[i]   # num of obs with equal regressors
    
    # compute the initial values of the parameters (initial urns' composition)
    F0omega = getF0omega(times=tt,
                         W=matrix(W.unique[i,],nr=length(tt),ncol=ncol(W),byrow = TRUE),
                         k=k,
                         theta=theta,
                         delta=delta,
                         log_u=log_u,
                         dist=dist,
                         omega_m=omega_m)
    F0delta    <- exp(F0omega$F0delta)
    F0         <- exp(F0omega$F0)
    cumF0      <- rowSums(F0)
    omega      <- F0omega$omega; 
    A          <- t(cbind(1-cumF0,F0delta)*omega)
    
    # Start computing the predictive probabilities and reinforcing the urns,
    # each time adding the relevant contribution to the log-likelihood
    for(j in 1:ni.i){
      # Computes the contribution to the marginal log-likelihood
      sumA <- colSums(A)
      divA <- t(A)/sumA
      deltaFstar <- t(divA[,2:(k+1)]*c(1,cumprod(divA[1:(nrow(divA)-1),1])))
      t.index <- match(times.i[j],tt)
      if(types.i[j]!=0){ # Case: uncensored observation
        loglik <- loglik + log(deltaFstar[types.i[j],t.index])
      } else { # Case: censored observation
        Fstar <- cumsum(colSums(deltaFstar))[t.index]
        loglik <- loglik + log(1-Fstar)
      }
      
      # Updates the urn compositions according to the update rules of Theorem 4.1
      A[types.i[j]+1,t.index] <- A[types.i[j]+1,t.index]+1
      A[1,] <- A[1,] + as.numeric(times.i[j]>tt)
    }
  }
  if(is.na(loglik) | is.infinite(loglik)) loglik = -Inf
  return(loglik)
}


# function computing the log-posterior distribution for the subdistribution beta-stacy model
# INPUTS: theta   = p(2k-1)+k vec of parameters (b1,b2,...,bk-1,v1,...,vk,u1,...,uk)
#         times   = n vec of event times
#         types   = n vec of event types (0 for censored)
#         k       = n of event types
#         W       = nxp matrix of regressors
#         b0      = px(k-1) matrix of prior means of b
#         Sb      = pxpx(k-1) array of prior cov of b
#         v0      = pxk matrix of prior means of v
#         Sv      = pxpxk array of prior cov of v
#         pu      = pxk vec of hyper alphas for u
#         qu      = pxk vec of hyper betas for u
#         log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
#         verbose = if FALSE (default) only return the value of the log-posterior distribution (the kernel of),
#                   if TRUE also returns the log-likelihood and log-hyperpriors in a list
#         dist    = centering distribution of event times; either "weibull" or "lognormal"
#         omega_m = scalar used to control the prior precision; divides the values of omega
logposterior <-      function(theta=NULL,
                              times=NULL,
                              delta=1,
                              tau=NULL,
                              types=NULL,
                              k=NULL,
                              W=NULL,
                              b0=NULL,
                              Sb=NULL,
                              v0=NULL,
                              Sv=NULL,
                              pu=NULL,
                              qu=NULL,
                              verbose=FALSE,
                              log_u=FALSE,
                              omega_m=1,
                              dist="weibull"){
  p <- ncol(W)   # number of regressors
  #logprior
  if(dist=="weibull"){
    # scomposition of theta in (b,v,u)
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p)
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p)
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)]
    if(log_u){ u <- exp(u) }
    # log prior for b
    logpriorb <- sum(sapply(1:(k-1),function(i)
      mvtnorm::dmvnorm(b[,i],b0[,i],as.matrix(Sb[,,i]),log=TRUE)))
    # log prior for v
    logpriorv <- sum(sapply(1:k,function(i)
      mvtnorm::dmvnorm(v[,i],v0[,i],as.matrix(Sv[,,i]),log=TRUE)))
    # log prior for u
    logprioru <- sum(dgamma(u,pu,qu,log=TRUE))
    # total weibull logprior
    logprior <- logpriorb + logpriorv + logprioru
  }
  if(dist=="lognormal"){
    # scomposition of theta in (b,v,u)
    b <- theta[1:(p*(k-1))]; b <- matrix(b,nr=p)
    v <- theta[(p*(k-1)+1):(p*(2*k-1))]; v <- matrix(v,nr=p)
    u <- theta[(p*(2*k-1)+1):(p*(2*k-1)+k)]
    if(log_u){ u <- exp(u) }
    # log prior for b
    logpriorb <- sum(sapply(1:(k-1),function(i)
      mvtnorm::dmvnorm(b[,i],b0[,i],as.matrix(Sb[,,i]),log=TRUE)))
    # log prior for v
    logpriorv <- sum(sapply(1:k,function(i)
      mvtnorm::dmvnorm(v[,i],v0[,i],as.matrix(Sv[,,i]),log=TRUE)))
    # log prior for u
    dinvg <- sapply(1:length(u), function(i) dinvgamma(u[i]^2,pu[i],qu[i]))
    logprioru <- sum(log(dinvg)+log(u))
    # total lognormal log prior
    logprior <- logpriorb + logpriorv + logprioru
  }
    
  # loglikelihood
  loglik <- loglik(times,types,k,W,tau,delta,theta,log_u,dist=dist,omega_m=omega_m)

  # log posterior
  logpost <- loglik + logprior
  if(log_u){ logpost <- logpost + sum(log(u)) } #adds the Jacobian when considering with log-transformed u's
  
  # output
  if(verbose==FALSE) return(logpost = logpost)  
  if(verbose==TRUE) return(list(
    logpost   = logpost,
    loglik    = loglik,
    logpriorb = logpriorb,
    logpriorv = logpriorv,
    logprioru = logprioru
  ))
}


# Function to compute the Discrete-time Nelson-Aalen (A) and Kaplan-Meier (K) estimators
# of the cause-specific cumulative hazard and overall survival functions, as well
# as the corresponding (frequentist) non-parametric estimator (F.hat) of the subdistribution
# function over the time interval from 0 to tau.
# INPUTS
# tau    = (scalar) upper bound for the event times
# times  = vector of observed event times
# k      = n of event types
# Wpred  = 1xp matrix of regressors of the target subdistribution function
# Wobs   = nxp matrix of observed regressors
# log_u   = theta contains the values u1,...,uk (FALSE; default) or log(u1),...,log(uk) (TRUE)
freqest <- function(tau=NULL,times=NULL,delta=1,types=NULL,k=NULL,Wpred=NULL,Wobs=NULL){
  p <- ncol(Wpred) #number of predictors
  times_est <- seq(0,tau,by=delta)
  L <- length(times_est)
  # for safety, checks that Wobs and Wpred have the same number of colums
  if(ncol(Wobs)!=p) stop("Wpred and Wobs must have the same number of columns.") 
  # selects only the observations with regressors equal to Wpred
  I <- which(apply(Wobs,1,function(r){ all(r==Wpred) }))
  times <- times[I]
  types <- types[I]
  # computes the posterior expectation
  A      <- matrix(0,nrow=L,ncol=k)
  S      <- matrix(0,nrow=L,ncol=1)
  F.hat  <- matrix(0,nrow=L,ncol=k)
  l      <- rep(0,times=L)
  m      <- matrix(0,ncol=(k+1),nrow=L)
  den    <- matrix(0,nrow=L,ncol=1)
  S_ <- 1
  for(tt in 1:L){
    m[tt,] <- sapply(0:k,function(d) sum((times==times_est[tt])&(types==d)) )
    l[tt] <- sum(times>times_est[tt])
    numS <- sum(m[tt,-1])
    den[tt]  <- l[tt] + sum(m[tt,]) 
    S_ <- S_ * (1 - numS/den[tt])
    S[tt] <- S_
  }
  Slag   <- c(1,S[1:length(S)-1])
  for(C in 1:k){
    A_ <- 0
    for(tt in 1:L){
      numA <- m[tt,C+1]
      A_ <- A_ + numA/den[tt]
      A[tt,C] <- A_
    }
    Adelta <- c(A[1,C],diff(A[,C]))
    F.hat[,C] <- cumsum(Slag*Adelta)
  }
  
  # output
  return(list(A=A,S=S,F.hat=F.hat))
}

#####################################################################################
# Functions to work with random subdistribution functions generated according
# to the subdistribution beta-Stacy process distribution
#####################################################################################

# generates a random discrete-time sub-distribution function
# from a sub-beta-stacy with centering subdistribution function (F0)
# and vector of precision parameters
# INPUTS:
# F0    = tau x k matrix. the (t,c) element of this matrix if F0(t,c), i.e.
#         the value of F0 at time t (t=1,...,tau) and event type c (c=1,...,k)
# omega = vector of length tau whose t-th element is the precision parameter
#         omega_t corresponding to time t (t=1,...,tau)
sbstacy <- function(omega,F0){
  # accounts for possible infinite omegas
  omega[is.infinite(omega)] <- 10^6
  alpha <- matrix(0,ncol=ncol(F0)+1,nrow=nrow(F0))
  deltaF0 <- F0
  deltaF0[2:nrow(F0),] <- F0[2:nrow(F0),]-F0[1:(nrow(F0)-1),]
  alpha[,1] <- omega*(1-rowSums(F0))
  # accounts for the fact that 1-rowSums(F.post) may be negative by numerical error if 
  # rowSums(F.post) is very close to 1
  alpha[alpha<0] <- 0
  alpha[,2:ncol(alpha)] <- omega*deltaF0
  alpha[1,] <- rep(1,times=ncol(alpha))
  W <- apply(alpha[-1,],1,FUN=function(a) myrdirichlet(1,a)) # does not consider the first row of alpha since at that point the subdistribution function is fixed at zero
  w0 <- cumprod(W[1,])
  w0[2:length(w0)] <- w0[1:(length(w0)-1)]
  w0[1] <- 1
  w0 <- matrix(w0,nrow=nrow(W)-1,ncol=ncol(W),byrow=TRUE)
  P <- W[2:nrow(W),]*w0
  P[is.na(P)] <- 0
  subF <- apply(P,1,FUN=cumsum)
  subF <- rbind(rep(0,times=ncol(subF)),subF) # adjust so that subF(0,c)=0 for all c
  return(subF)
}

### evaluates a sub-distribution function represented by the matrix subF 
eval.subdf <- function(subF,time,cause){
  ifelse(time<=0,0,subF[time,cause])
}

### Function to simulate data from a subdistribution function 
# Inputs:
# n = number of samples to generate
# tau = grid of discrete time points 
# subF = matrix representing the subdistribution function
sample.subdf <- function(n,tau,subF){
  P <- rep(0,times=nrow(subF)+1)
  P[1:nrow(subF)] <- rowSums(subF)
  P[nrow(subF)+1] <- 1
  w <- pmax(0,diff(P))
  T <- sample(1:length(w),size=n,prob = w,replace=TRUE)
  times <- tau[T+1]
  times[T==length(w)] <- tau[length(w)]
  C <- vector("numeric",length=n)
  for(i in 1:n){
    if(T[i]==length(w)){
      C[i] <- 0
    } else {
      C[i] <- sample(1:ncol(subF),size=1,prob=(subF[T[i]+1,]-subF[T[i],])/w[T[i]])
    }
  }
  sim <- cbind(times,C)
  return(sim)
}


#####################################################################################
# Compiles all functions (gives a small speed gain)
#####################################################################################

rgamss              <- compiler::cmpfun(rgamss)
myrdirichlet        <- compiler::cmpfun(myrdirichlet)
getF0omega          <- compiler::cmpfun(getF0omega)
logposterior_param  <- compiler::cmpfun(logposterior_param)
postest             <- compiler::cmpfun(postest)
loglik              <- compiler::cmpfun(loglik)
logposterior        <- compiler::cmpfun(logposterior)
sbstacy             <- compiler::cmpfun(sbstacy)
freqest             <- compiler::cmpfun(freqest)

