
source("C:/Users/andre/Dropbox/sub-beta-stacy/programs/functions.R")
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/plot_pred_dist.Rdata")
set.seed(314271)
library(timereg)
data("melanoma")

times <-  seq(0,7000,by=1)

############ WEIBULL - FEMALES ############################

# Simulates from the prior
k <- 2
p <- 2
W <- matrix(c(1,0),ncol = 2,nrow=length(times),byrow = TRUE) # (sex=0 for females)
m <- c(1e0,1e3,1e6)
niter <- 50
fp <- array(0,dim = c(niter,length(times),length(m)))
fs <- array(0,dim = c(niter,length(times),length(m)))
for(i in 1:niter){
  print(i)
  # Prior
  b <- rnorm(2,0,1)
  v1 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
  v2 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
  u <- rgamma(2,11,10)
  theta <- c(b,v1,v2,u)
  for(j in 1:length(m)){
    lf0 <- getF0omega(theta=theta, delta=1, times = times, W=W, 
                      k=2,log_u = FALSE, omega_m=m[j])
    f0 <- sbstacy(omega=lf0$omega,F0 = exp(lf0$F0))
    fs[i,,j] <- f0[,1]
    fp[i,,j] <- exp(lf0$F0)[,1]
  }
}

# Extracts a sample of curves from the posterior samples
I.1 <- sample(1:length(post.female.sbs.m1),niter,replace = FALSE)
I.1000 <- sample(1:length(post.female.sbs.m1000),niter,replace = FALSE)
I.1e6 <- sample(1:length(post.female.sbs.m1e6),niter,replace = FALSE)

pdf("d:/prior_sim_F.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,fs[1,,1],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,main=main,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,1],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.female.sbs.m1[[I.1[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 

main <- expression(bold(log[10](m)==3))
plot(tt,fs[1,,2],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,main=main,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,2],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.female.sbs.m1000[[I.1000[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 

main <- expression(bold(log[10](m)==6))
plot(tt,fs[1,,3],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,main=main,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,3],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.female.sbs.m1e6[[I.1e6[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 
dev.off()

############ WEIBULL - MALES ############################

# Simulates from the prior
k <- 2
p <- 2
W <- matrix(c(1,1),ncol = 2,nrow=length(times),byrow = TRUE) # (sex=1 for males)
m <- c(1e0,1e3,1e6)
niter <- 50
fp <- array(0,dim = c(niter,length(times),length(m)))
fs <- array(0,dim = c(niter,length(times),length(m)))
for(i in 1:niter){
  print(i)
  # Prior
  b <- rnorm(2,0,1)
  v1 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
  v2 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
  u <- rgamma(2,11,10)
  theta <- c(b,v1,v2,u)
  for(j in 1:length(m)){
    lf0 <- getF0omega(theta=theta, delta=1, times = times, W=W, 
                      k=2,log_u = FALSE, omega_m=m[j])
    f0 <- sbstacy(omega=lf0$omega,F0 = exp(lf0$F0))
    fs[i,,j] <- f0[,1]
    fp[i,,j] <- exp(lf0$F0)[,1]
  }
}

# Extracts a sample of curves from the posterior samples
I.1 <- sample(1:length(post.male.sbs.m1),niter,replace = FALSE)
I.1000 <- sample(1:length(post.male.sbs.m1000),niter,replace = FALSE)
I.1e6 <- sample(1:length(post.male.sbs.m1e6),niter,replace = FALSE)

pdf("d:/prior_sim_M.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,fs[1,,1],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,1],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.male.sbs.m1[[I.1[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 

main <- expression(bold(log[10](m)==3))
plot(tt,fs[1,,2],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,2],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.male.sbs.m1000[[I.1000[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 

main <- expression(bold(log[10](m)==6))
plot(tt,fs[1,,3],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2,
     xlab="Days since surgery",ylab="Cumulative incidence")
for(i in 2:niter){
  points(tt,fs[i,,3],ylim=c(0,1),type="s",col=rgb(0.5,0.5,0.5,alpha=0.4),lwd=2)
  points(tt,post.male.sbs.m1e6[[I.1e6[i]]][,1],lwd=2,type="s",col=rgb(0,0,0,alpha=0.6))
} 
dev.off()



# ############ LOG-NORMAL ############################
# 
# k <- 2
# p <- 2
# 
# # build useful quantities
# W <- matrix(c(1,0),ncol = 2,nrow=length(times),byrow = TRUE) # (sex=0 for females)
# m <- c(1e3)
# niter <- 100
# fp <- array(0,dim = c(niter,length(times),length(m)))
# fs <- array(0,dim = c(niter,length(times),length(m)))
# for(i in 1:niter){
#   print(i)
#   # Prior
#   b <- rnorm(2,0,1)
#   v1 <- c(rnorm(1,log(3650),1),rnorm(1,0,1))
#   v2 <- c(rnorm(1,log(3650),1),rnorm(1,0,1))
#   u <- rgamma(2,2,1)
#   theta <- c(b,v1,v2,u)
#   for(j in 1:length(m)){
#     lf0 <- getF0omega(theta=theta, delta=1, times = times, W=W, 
#                       k=2,log_u = FALSE, omega_m=m[j], dist="lognormal")
#     f0 <- sbstacy(omega=lf0$omega,F0 = exp(lf0$F0))
#     fs[i,,j] <- f0[,1]
#     fp[i,,j] <- exp(lf0$F0)[,1]
#   }
# }
# 
# plot(times,0*times,col="white",ylim=c(0,1),xlab="Time since surgery",ylab="Standard deviation")
# for(j in 1:niter) points(times,fs[j,,1],type="s")
# points(times,apply(fp, c(2,3), mean),type="s",col="red")
