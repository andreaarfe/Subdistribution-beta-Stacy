
library(timereg)
source("C:/Users/andre/Dropbox/sub-beta-stacy/programs/functions.R")
set.seed(314271)
data("melanoma")

# Simulate from the prior distribution to estimate:
# 1) the prior mean
# 2) the prior variance Var(F(t,c)-F0(t,c,|theta))=E[Var(F(t,c)|theta)] 
# Var(F(t,c)|theta) is computed using the formulas in Lemma 2.1 of the main text
# Focuses on women as the reference groups.
times <- seq(1,7000,by=1)
tt <- c(0,times)
k <- 2
p <- 2
W <- matrix(c(1,0),ncol = 2,nrow=length(times),byrow = TRUE) # (sex=0 for females)
m <- c(1e0,1e1,1e2,1e3,1e4,1e5)
niter <- 1000
meanF <- array(0,dim=c(length(m),niter,k,length(times)))
varF <- array(0,dim=c(length(m),niter,k,length(times)))
prob <- array(0,dim=c(length(m),niter))
for(j in 1:length(m)){
  for(i in 1:niter){
    # Prior
    b <- rnorm(2,0,1)
    v1 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
    v2 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
    u <- rgamma(2,11,10)
    theta <- c(b,v1,v2,u)
    logit <- W[1,]%*%b
    prob[j,i] <- exp(logit)/(1+exp(logit))
    lf0 <- getF0omega(theta=theta, delta=1, times = times, W=W, 
                      k=2,log_u = FALSE, omega_m=m[j])
    F0delta    <- exp(lf0$F0delta)
    F0         <- exp(lf0$F0)
    cumF0      <- rowSums(F0)
    omega      <- lf0$omega; 
    A          <- t(cbind(1-cumF0,F0delta)*omega) 
    sumA <- colSums(A)
    #divA <- t(t(A)/sumA)
    #mean_deltaF <- t(divA[2:(k+1),])*c(1,cumprod(divA[1,1:(ncol(divA)-1)]))
    divA_plus1 <- t((1+t(A))/(1+sumA))
    mean_deltaF2_plus1 <- t(divA_plus1[2:(k+1),])*c(1,cumprod(divA_plus1[1,1:(ncol(divA_plus1)-1)]))
    meanF[j,i,,] <- t(F0)
    varF[j,i,,] <- t(F0delta*(mean_deltaF2_plus1 - F0delta))
    #varF[j,i,,] <- mean_deltaF*(mean_deltaF2_plus1 - mean_deltaF)
  }
}
mm <- apply(meanF, c(3,4), mean)
vv <- apply(varF, c(1,3,4), mean)

# Plot of the prior mean
plot(mm[1,],type="s",xlab="Time since surgery",
     ylab=expression(paste(E,"[",F[0],"(t,1| ",theta,")]")))
plot(mm[2,],type="s",xlab="Time since surgery",
     ylab=expression(paste(E,"[",F[0],"(t,2| ",theta,")]")))

# Plot of prior conditional variance E(var(F|theta))=Var(F-F_0(theta)) 
# for given values of m
pdf(file="d:/prior_stddev.pdf")
par(mar=c(5,5,3,3))
tt <- c(0,times)
plot(tt,c(0,sqrt(vv[1,1,])),ylim=c(0,0.026),type="s",lty=2,lwd=2,
     xlab="Days since surgery",
     ylab=expression(sqrt(Var(paste(Delta,F)(t,1)-paste(Delta,F)[0](t,paste("1| ", theta))))))
points(tt,c(0,sqrt(vv[2,1,])),type="s",lty=3,lwd=2)
points(tt,c(0,sqrt(vv[3,1,])),type="s",lty=4,lwd=2)
points(tt,c(0,sqrt(vv[4,1,])),type="s",lty=5,lwd=2)
points(tt,c(0,sqrt(vv[5,1,])),type="s",lty=6,lwd=2)
points(tt,c(0,sqrt(vv[6,1,])),type="s",lty=8,lwd=2)
points(tt,c(0,sqrt(diff(c(0,mm[1,]))*(1-diff(c(0,mm[1,]))))),type="s",lwd=3)
legend("topright",lty=c(2,3,4,5,6,8,1),lwd=c(2,2,2,2,2,2,3),bty="n",
       legend=c(expression(log[10](m)==0),
                expression(log[10](m)==1),
                expression(log[10](m)==2),
                expression(log[10](m)==3),
                expression(log[10](m)==4),
                expression(log[10](m)==5),
                expression(log[10](m)==+infinity)))
dev.off()

# Same plot but for men.
times <- seq(1,7000,by=1)
tt <- c(0,times)
k <- 2
p <- 2
W <- matrix(c(1,1),ncol = 2,nrow=length(times),byrow = TRUE) # (sex=1 for men)
m <- c(1e0,1e1,1e2,1e3,1e4,1e5)
niter <- 1000
meanF <- array(0,dim=c(length(m),niter,k,length(times)))
varF <- array(0,dim=c(length(m),niter,k,length(times)))
prob <- array(0,dim=c(length(m),niter))
for(j in 1:length(m)){
  for(i in 1:niter){
    # Prior
    b <- rnorm(2,0,1)
    v1 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
    v2 <- c(rnorm(1,log(log(2)/3650),1),rnorm(1,0,1))
    u <- rgamma(2,11,10)
    theta <- c(b,v1,v2,u)
    logit <- W[1,]%*%b
    prob[j,i] <- exp(logit)/(1+exp(logit))
    lf0 <- getF0omega(theta=theta, delta=1, times = times, W=W, 
                      k=2,log_u = FALSE, omega_m=m[j])
    F0delta    <- exp(lf0$F0delta)
    F0         <- exp(lf0$F0)
    cumF0      <- rowSums(F0)
    omega      <- lf0$omega; 
    A          <- t(cbind(1-cumF0,F0delta)*omega) 
    sumA <- colSums(A)
    divA_plus1 <- t((1+t(A))/(1+sumA))
    mean_deltaF2_plus1 <- t(divA_plus1[2:(k+1),])*c(1,cumprod(divA_plus1[1,1:(ncol(divA_plus1)-1)]))
    meanF[j,i,,] <- t(F0)
    varF[j,i,,] <- t(F0delta*(mean_deltaF2_plus1 - F0delta))
  }
}
mm <- apply(meanF, c(3,4), mean)
vv <- apply(varF, c(1,3,4), mean)

# Plot of the prior mean
plot(mm[1,],type="s",xlab="Time since surgery",
     ylab=expression(paste(E,"[",F[0],"(t,1| ",theta,")]")))
plot(mm[2,],type="s",xlab="Time since surgery",
     ylab=expression(paste(E,"[",F[0],"(t,2| ",theta,")]")))

# Plot of prior conditional variance E(var(F|theta))=Var(F-F_0(theta)) 
# for given values of m
pdf(file="d:/prior_stddev_M.pdf")
par(mar=c(5,5,3,3))
plot(tt,c(0,sqrt(vv[1,1,])),ylim=c(0,0.026),type="s",lty=2,lwd=2,
     xlab="Days since surgery",
     ylab=expression(sqrt(Var(paste(Delta,F)(t,1)-paste(Delta,F)[0](t,paste("1| ", theta))))))
points(tt,c(0,sqrt(vv[2,1,])),type="s",lty=3,lwd=2)
points(tt,c(0,sqrt(vv[3,1,])),type="s",lty=4,lwd=2)
points(tt,c(0,sqrt(vv[4,1,])),type="s",lty=5,lwd=2)
points(tt,c(0,sqrt(vv[5,1,])),type="s",lty=6,lwd=2)
points(tt,c(0,sqrt(vv[6,1,])),type="s",lty=8,lwd=2)
points(tt,c(0,sqrt(diff(c(0,mm[1,]))*(1-diff(c(0,mm[1,]))))),type="s",lwd=3)
legend("topright",lty=c(2,3,4,5,6,8,1),lwd=c(2,2,2,2,2,2,3),bty="n",
       legend=c(expression(log[10](m)==0),
                expression(log[10](m)==1),
                expression(log[10](m)==2),
                expression(log[10](m)==3),
                expression(log[10](m)==4),
                expression(log[10](m)==5),
                expression(log[10](m)==+infinity)))
dev.off()



#save.image("d:/prior_mean_variance.Rdata")
#load("d:/prior_mean_variance.Rdata")
