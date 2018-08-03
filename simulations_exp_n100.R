
start <- Sys.time()

library(nnet)
library(timereg)
library(compiler)
library(parallel)
set.seed(12345)
enableJIT(3)
#source("C:/Users/andre/Dropbox/sub-beta-stacy/programs/functions.R")
source("functions.R")

# sample from the posterior distributions of the models
samp.post.par <- function(trace=NULL,W=NULL,times=NULL,
                          delta=NULL,k=NULL,dist="weibull",
                          log_u=FALSE){
  samp <- lapply(1:nrow(trace), function(i){
    l <- getF0omega(theta=trace[i,],
                    times=times,
                    delta=delta,
                    k=k,
                    dist=dist,
                    log_u=log_u,
                    W=W) 
    return(exp(l$F0))
  })
  return(samp)
}

samp.post.sbs <- function(trace=NULL,
                          tau=NULL,
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
  samp <- lapply(1:nrow(trace), function(i){
    l <- postmean(theta=trace[i,],
                  tau=tau,
                  times=times,
                  delta=delta,
                  types=types,
                  k=k,
                  Wpred=Wpred,
                  Wobs=Wobs,
                  omega_m=omega_m,
                  log_u=log_u,
                  dist=dist) 
    return(l)
  })
  return(samp)
}

###############################################################################

#### Simulation parameters
nsamp <- 100       # Sample size
nsims <- 100       # Number of simulated samples

# MCMC parameters
nburnin = 1000     # num of burned iterations
niter   = 1000     # num of iterations after burn in
thin    = 10       # thinning

#### Estimates data-generating parameters for the Weibull model from the 
#### melanoma data

# load data
data("melanoma")

# build useful quantities
W <- matrix(1,nrow=nrow(melanoma))    # cov matrix (sex=0 for females)
times = melanoma$days                # event times
types = melanoma$status              # event types:      
types[types==2] = 0                  # 0. censored
types[types==3] = 2                  # 1. dead from melanoma
                                     # 2. dead from other causes
delta = 1                            
tau <- 7000

# Preliminary parametric estimates, used as initialization values 
I <- which(types!=0)
tt <- seq(0,tau,delta)
int <- findInterval(times,tt,all.inside = TRUE)
f <- relevel(as.factor(types[I]),ref=2)
ww <- 1/predict(glm(types!=0 ~ -1+W ))
b <- coef(multinom(f~-1+W[I,],weights = ww[I]))
uv1 <- survreg(Surv(tt[int]-delta,tt[int],event=3*(types==1),type = "interval")~-1+W,dist="weibull")
uv2 <- survreg(Surv(tt[int]-delta,tt[int],event=3*(types==2),type = "interval")~-1+W,dist="weibull")
u1 <- 1/uv1$scale
u2 <- 1/uv2$scale
v1 <- -coef(uv1)*u1
v2 <- -coef(uv2)*u2
theta.init <- c(b,v1,v2,log(u1),log(u2))

# Maximum likelihood estimates from the lognormal model
mle <- optim(fn         = function(theta,...) -loglik_param(theta,...),
             par        = theta.init,
             times      = times,
             delta      = 1,
             types      = types,
             k          = 2,
             W          = W,
             dist       = "weibull",
             log_u      = TRUE,
             method     = "BFGS",
             hessian    = TRUE,
             control=list(trace=1))

#### simulation parameters
dist <- "weibull"
tau <- 7000
theta.true <- mle$par
delta <- 1
k <- 2
times <- seq(0,tau,delta)
W <- matrix(1,nrow=length(times),ncol=1)
deltaF0 <- exp(getF0omega(theta=theta.true,
                  times=times,
                  W=W,
                  k=2,
                  dist=dist,
                  log_u=TRUE)$F0delta) # "true" subdistribution function
deltaS <- rowSums(deltaF0)
condevent <- deltaF0 / deltaS

# Run the simulation
sim <- function(i){
  # Simulates time to events
  tobs <- sample(c(times,Inf),
                 size=nsamp,
                 replace = TRUE,
                 prob = c(deltaS,max(0,1-sum(deltaS))))
  
  # Simulates event type/ censoring indicators anc censoring times
  event <- rep(0,times=nsamp)
  I <- which(is.finite(tobs))
  for(i in I){
    J <- which(times==tobs[i])
    event[i] <- sample(c(1,2),1,prob=condevent[J,])
  }
  tobs <- pmin(tobs,tau)
  tryCatch(
    {
    # Fit the Bayesian model using the same settings and priors from the main analysis
    # but forces the parametric functional form to be exponential
    W <- matrix(1,nrow=nsamp,ncol=1)
    I <- which(event!=0)
    tt <- seq(0,tau,delta)
    int <- findInterval(tobs,tt,all.inside = TRUE,left.open = TRUE)
    f <- relevel(as.factor(event[I]),ref=2)
    ww <- 1/predict(glm(event!=0 ~ -1+W ))
    b <- coef(multinom(f~-1+W[I,],weights = ww[I]))
    uv1 <- survreg(Surv(pmax(0.1,tt[int]),tt[int]+delta,event=3*(event==1),type = "interval")~-1+W,dist="exponential")
    uv2 <- survreg(Surv(pmax(0.1,tt[int]),tt[int]+delta,event=3*(event==2),type = "interval")~-1+W,dist="exponential")
    u1 <- 1
    u2 <- 1
    v1 <- -coef(uv1)*u1
    v2 <- -coef(uv2)*u2
    theta.init <- c(b,v1,v2,log(u1),log(u2))
    k  <- length(unique(types))-1              # number of event types
    p  <- ncol(W)                              # number of regressors  
    b0 <- matrix(0,p,k-1)                      # b prior mean
    Sb <- array(diag(1,p),c(p,p,k-1))          # b prior variance
    v0 <- matrix(c(log(log(2)/3650),0),p,k)    # v prior mean
    Sv <- array(diag(1,p),c(p,p,k))            # v prior variance
    pu <- rep(1+10,k)                          # u alpha (ignored)
    qu <- rep(10,k)                            # u beta (ignored)
    opt <- optim(fn         = function(theta,...){ -logposterior_param(c(theta[1:(length(theta)-k)],rep(0,k)),...)},
                 par        = theta.init,
                 times      = tobs,
                 delta      = delta,
                 types      = event,
                 k          = k,
                 W          = W,
                 b0         = b0,
                 Sb         = Sb,
                 v0         = v0,
                 Sv         = Sv,
                 pu         = pu,
                 qu         = qu,
                 dist       ="weibull",
                 verbose    = FALSE,
                 log_u      = TRUE,
                 method     = "BFGS",
                 hessian    = TRUE,
                 control=list(trace=1))
    tune <- 2.4/sqrt(p*(2*k-1)+k)
    H <- opt$hessian
    d <- diag(H)
    d[(length(theta.init)-k+1):length(theta.init)] <- 1
    diag(H) <- d
    V    <- solve(opt$hessian+diag(1e-4,ncol=ncol(H),nrow=nrow(H)))
    # Lognormal parametric model
    out <- MCMCpack::MCMCmetrop1R(fun        = function(theta,...){ logposterior_param(c(theta[1:(length(theta)-k)],rep(0,k)),...)},
                                  theta.init = opt$par,
                                  burnin     = nburnin,
                                  mcmc       = niter,
                                  thin       = thin,
                                  tune       = tune,
                                  times      = tobs,
                                  delta      = delta,
                                  types      = event,
                                  k          = k,
                                  W          = W,
                                  b0         = b0,
                                  Sb         = Sb,
                                  v0         = v0,
                                  Sv         = Sv,
                                  pu         = pu,
                                  qu         = qu,
                                  dist       = "weibull",
                                  verbose    = 100,
                                  log_u      = TRUE,
                                  V          = V)  
    # Nonparametric subdistribution beta-Stacy model 
    out.sbs.1 <- MCMCpack::MCMCmetrop1R(fun      = function(theta,...){ logposterior(c(theta[1:(length(theta)-k)],rep(0,k)),...)},
                                      theta.init = opt$par,
                                      burnin     = nburnin,
                                      mcmc       = niter,
                                      thin       = thin,
                                      tune       = tune,
                                      times      = tobs,
                                      delta      = delta,
                                      tau        = tau,
                                      types      = event,
                                      k          = k,
                                      W          = W,
                                      b0         = b0,
                                      Sb         = Sb,
                                      v0         = v0,
                                      Sv         = Sv,
                                      pu         = pu,
                                      qu         = qu,
                                      omega_m    = 1,
                                      dist       = "weibull",
                                      verbose    = 100,
                                      log_u      = TRUE,
                                      V          = V)
    out.sbs.1000 <- MCMCpack::MCMCmetrop1R(fun     = function(theta,...){ logposterior(c(theta[1:(length(theta)-k)],rep(0,k)),...)},
                                        theta.init = opt$par,
                                        burnin     = nburnin,
                                        mcmc       = niter,
                                        thin       = thin,
                                        tune       = tune,
                                        times      = tobs,
                                        delta      = delta,
                                        tau        = tau,
                                        types      = event,
                                        k          = k,
                                        W          = W,
                                        b0         = b0,
                                        Sb         = Sb,
                                        v0         = v0,
                                        Sv         = Sv,
                                        pu         = pu,
                                        qu         = qu,
                                        omega_m    = 1000,
                                        dist       = "weibull",
                                        verbose    = 100,
                                        log_u      = TRUE,
                                        V          = V)
    out.sbs.1e6 <- MCMCpack::MCMCmetrop1R(fun         = function(theta,...){ logposterior(c(theta[1:(length(theta)-k)],rep(0,k)),...)},
                                           theta.init = opt$par,
                                           burnin     = nburnin,
                                           mcmc       = niter,
                                           thin       = thin,
                                           tune       = tune,
                                           times      = tobs,
                                           delta      = delta,
                                           tau        = tau,
                                           types      = event,
                                           k          = k,
                                           W          = W,
                                           b0         = b0,
                                           Sb         = Sb,
                                           v0         = v0,
                                           Sv         = Sv,
                                           pu         = pu,
                                           qu         = qu,
                                           omega_m    = 1e6,
                                           dist       = "weibull",
                                           verbose    = 100,
                                           log_u      = TRUE,
                                           V          = V)
    # Gets the posterior means
    out[,4:5] <- 0
    post.par <- samp.post.par(trace=out,W=matrix(1,ncol=1,nrow=length(tt),byrow=TRUE),
                              times=tt,delta=delta,k=2,dist="weibull",log_u=TRUE)
    post.par.1 <- sapply(post.par,function(l) l[,1])
    pred.par <- rowMeans(post.par.1)
    
    out.sbs.1[,4:5] <- 0
    post.sbs.1 <- samp.post.sbs(trace=out.sbs.1,
                              tau=tau,
                              Wpred=matrix(1,ncol=1,nrow=length(tt),byrow=TRUE),
                              Wobs=W,
                              times=tobs,
                              types=event,
                              delta=delta,
                              omega_m=1,
                              k=2,
                              dist="weibull",
                              log_u=TRUE)
    post.sbs.1.1 <- sapply(post.sbs.1,function(l) l[,1])
    pred.sbs.1 <- rowMeans(post.sbs.1.1)
    
    out.sbs.1000[,4:5] <- 1
    post.sbs.1000 <- samp.post.sbs(trace=out.sbs.1000,
                                tau=tau,
                                Wpred=matrix(1,ncol=1,nrow=length(tt),byrow=TRUE),
                                Wobs=W,
                                times=tobs,
                                types=event,
                                delta=delta,
                                omega_m=1000,
                                k=2,
                                dist="weibull",
                                log_u=TRUE)
    post.sbs.1.1000 <- sapply(post.sbs.1000,function(l) l[,1])
    pred.sbs.1000 <- rowMeans(post.sbs.1.1000)
    
    out.sbs.1e6[,4:5] <- 0
    post.sbs.1e6 <- samp.post.sbs(trace=out.sbs.1e6,
                                tau=tau,
                                Wpred=matrix(1,ncol=1,nrow=length(tt),byrow=TRUE),
                                Wobs=W,
                                times=tobs,
                                types=event,
                                delta=delta,
                                omega_m=1e6,
                                k=2,
                                dist="weibull",
                                log_u=TRUE)
    post.sbs.1.1e6 <- sapply(post.sbs.1e6,function(l) l[,1])
    pred.sbs.1e6 <- rowMeans(post.sbs.1.1e6)
    
    # Computes the Kolmogorov-Smirnov statistics
    F0.1         <- cumsum(deltaF0[,1])
    KS.par       <- max(abs(pred.par      - F0.1))
    KS.sbs.1     <- max(abs(pred.sbs.1    - F0.1))
    KS.sbs.1000  <- max(abs(pred.sbs.1000 - F0.1))
    KS.sbs.1e6   <- max(abs(pred.sbs.1e6  - F0.1))
    
    # Compute the Cramer-von-Mises statistic
    CVM.par      <- sum( (pred.par      - F0.1)^2*delta )
    CVM.sbs.1    <- sum( (pred.sbs.1    - F0.1)^2*delta )
    CVM.sbs.1000 <- sum( (pred.sbs.1000 - F0.1)^2*delta )
    CVM.sbs.1e6  <- sum( (pred.sbs.1e6  - F0.1)^2*delta )
    
    return(list(KS.par=KS.par,
                KS.sbs.1=KS.sbs.1,
                KS.sbs.1000=KS.sbs.1000,
                KS.sbs.1e6=KS.sbs.1e6,
                CVM.par=CVM.par,
                CVM.sbs.1=CVM.sbs.1,
                CVM.sbs.1000=CVM.sbs.1000,
                CVM.sbs.1e6=CVM.sbs.1e6,
                tobs=tobs,
                event=event))
  },
  error=function(cond) return(list(KS.par=NA,
                                   KS.sbs.1=NA,
                                   KS.sbs.1000=NA,
                                   KS.sbs.1e6=NA,
                                   CVM.par=NA,
                                   CVM.sbs.1=NA,
                                   CVM.sbs.1000=NA,
                                   CVM.sbs.1e6=NA,
                                   tobs=tobs,
                                   event=event))
  # Outputs the results
)}
cl1<-makeCluster(2)
clusterExport(cl1, ls())
clusterEvalQ(cl1, library(nnet))
clusterEvalQ(cl1, library(timereg))
clusterEvalQ(cl1, library(MCMCpack))
out <- parLapply(cl1,1:nsims,sim)
stopCluster(cl1)

save(list=c("out"),file="out_exp_n100.Rdata")

end <- Sys.time()
end-start





