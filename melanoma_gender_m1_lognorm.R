# preliminaries
set.seed(123)
source("functions.R")

# load data
library(nnet)
library(timereg)
data("melanoma")

# build useful quantities
W <- cbind(1,melanoma[,c("sex")])    # cov matrix (sex=0 for females); includes colum for intercept term
times = melanoma$days                # event times
types = melanoma$status              # event types:      
types[types==2] = 0                  # 0. censored
types[types==3] = 2                  # 1. dead from melanoma
                                     # 2. dead from other causes
delta <- 1                           # Length of time intervals  
tau <- 7000                          # Maximum prediction time

# hyper-parameters specifications
k  <- 2                                            # number of event types
p  <- 2                                            # number of regressors (including the intercept) 
b0 <- matrix(0,p,k-1)                              # b prior mean
Sb <- array(diag(1,p),c(p,p,k-1))                  # b prior variance
v0 <- matrix(c(log(3650),0),p,k)                   # v prior mean
Sv <- array(diag(1,p),c(p,p,k))                    # v prior variance
pu <- rep(1+1,k)                                   # u alpha
qu <- rep(1,k)                                     # u beta
omega_m <- 1                                       # m (reinforcement parameter)

# standardizes the predictors to improve mixing
Wstd <- W
Wstd[,2:ncol(W)] <- scale(Wstd[,2:ncol(W)])
M <- diag(apply(W,2,sd))
M[1,] <- colMeans(W)
b0_std <- M%*%b0
Sb_std <- sapply(1:(dim(Sb)[3]),function(i) M%*%Sb[,,i]%*%t(M),simplify="array")
v0_std <- M%*%v0
Sv_std <- sapply(1:(dim(Sv)[3]),function(i) M%*%Sv[,,i]%*%t(M),simplify="array")

# Preliminary parametric estimates, used as initialization values for the MCMC algorithms
I <- which(types!=0)
tt <- seq(0,tau,delta)
int <- findInterval(times,tt,all.inside = TRUE)
f <- relevel(as.factor(types[I]),ref=2)
ww <- 1/predict(glm(types!=0 ~ -1+Wstd ))
b <- coef(multinom(f~-1+Wstd[I,],weights = ww[I]))
uv1 <- survreg(Surv(tt[int]-delta,tt[int],event=3*(types==1),type = "interval")~-1+Wstd,dist="lognormal")
uv2 <- survreg(Surv(tt[int]-delta,tt[int],event=3*(types==2),type = "interval")~-1+Wstd,dist="lognormal")
u1 <- uv1$scale
u2 <- uv2$scale
v1 <- coef(uv1)
v2 <- coef(uv2)
theta.init <- c(b,v1,v2,log(u1),log(u2))

# Find MAP to initialize MCMC sampler
opt <- optim(fn         = function(theta,...) -logposterior(theta,...),
             par        = theta.init,
             times      = times,
             delta      = delta,
             types      = types,
             k          = k,
             W          = Wstd,
             b0         = b0_std,
             Sb         = Sb_std,
             v0         = v0_std,
             Sv         = Sv_std,
             pu         = pu,
             qu         = qu,
             tau        = tau,
             omega_m    = omega_m,
             dist       = "lognormal",
             verbose    = FALSE,
             log_u      = TRUE,
             method     = "BFGS",
             hessian    = TRUE,
             control=list(trace=1))

# Covariance matrix of the Gaussian random walk proposal
tune <- 2.4/sqrt(p*(2*k-1)+k)
V    <- solve(opt$hessian)

### MCMC sampling from the marginal posterior distribution of the hyperparameters

# Tuning parameters for the MCMC sampler
nburnin = 1000      # num of burned iterations
niter   = 25000     # num of iterations after burn in
thin    = 25         # thinning

# run MCMC 
out <- MCMCpack::MCMCmetrop1R(fun        = logposterior,
                              theta.init = opt$par,
                              burnin     = nburnin,
                              mcmc       = niter,
                              thin       = thin,
                              tune       = tune,
                              times      = times,
                              delta      = delta,
                              tau        = tau,
                              types      = types,
                              k          = k,
                              W          = Wstd,
                              b0         = b0_std,
                              Sb         = Sb_std,
                              v0         = v0_std,
                              Sv         = Sv_std,
                              pu         = pu,
                              qu         = qu,
                              dist       = "lognormal",
                              omega_m    = omega_m,
                              verbose    = 100,
                              log_u      = TRUE,
                              V          = V)

# transform the parameters to account for scaling and reparametrization
lls <- c("b01","b11",
         "v01","v11",
         "v02","v12",
         "log_u1","log_u2")
colnames(out) <- parse(text=lls)
out_rescaled  <- out
# rescale b
out_rescaled[,"b11"]    <- out_rescaled[,"b11"]/sd(W[,2])
out_rescaled[,"b01"]    <- out_rescaled[,"b01"]-out_rescaled[,"b11"]*mean(W[,2])
# rescale v
out_rescaled[,"v11"]    <- out_rescaled[,"v11"]/sd(W[,2])
out_rescaled[,"v01"]    <- out_rescaled[,"v01"]-out_rescaled[,"v11"]*mean(W[,2])
out_rescaled[,"v12"]    <- out_rescaled[,"v12"]/sd(W[,2])
out_rescaled[,"v02"]    <- out_rescaled[,"v02"]-out_rescaled[,"v12"]*mean(W[,2])
# Rescale u
out_rescaled[,"log_u1"] <- exp(out_rescaled[,"log_u1"])
out_rescaled[,"log_u2"] <- exp(out_rescaled[,"log_u2"])
lls_rescaled <- c("b01","b11",
                  "v01","v11",
                  "v02","v12",
                  "u1","u2")
colnames(out_rescaled) <- parse(text=lls_rescaled)

# save output
save.image(file="melanoma_gender_m1_lognormal.RData")

