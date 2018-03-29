
source("functions.R")
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/plot_pred_dist.Rdata")
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


# Frequentist estimates
# Women
Wpred_W <- matrix(c(1,0),ncol=2)
freq_W <- freqest(tau=tau,times=times,delta=delta,types=types,k=k,Wobs=W,Wpred=Wpred_W)$F.hat
# Men
Wpred_M <- matrix(c(1,1),ncol=2)
freq_M <- freqest(tau=tau,times=times,delta=delta,types=types,k=k,Wobs=W,Wpred=Wpred_M)$F.hat
# Censoring distribution
cens <- freqest(tau=tau,
                times=times,
                delta=delta,
                types=(types==0),
                k=1,
                Wobs= matrix(c(1),ncol=1,nrow=length(times)),
                Wpred=matrix(c(1),ncol=1))$F.hat 

# Computes the predictive distributions
# Females
post.female.sbs.m1.1 <- sapply(post.female.sbs.m1,function(l) l[,1])
post.female.sbs.m1000.1 <- sapply(post.female.sbs.m1000,function(l) l[,1])
post.female.sbs.m1e6.1 <- sapply(post.female.sbs.m1e6,function(l) l[,1])
pred.female.sbs.m1.1 <- rowMeans(post.female.sbs.m1.1)
pred.female.sbs.m1000.1 <- rowMeans(post.female.sbs.m1000.1)
pred.female.sbs.m1e6.1 <- rowMeans(post.female.sbs.m1e6.1)
post.female.sbs.m1.2 <- sapply(post.female.sbs.m1,function(l) l[,2])
post.female.sbs.m1000.2 <- sapply(post.female.sbs.m1000,function(l) l[,2])
post.female.sbs.m1e6.2 <- sapply(post.female.sbs.m1e6,function(l) l[,2])
pred.female.sbs.m1.2 <- rowMeans(post.female.sbs.m1.2)
pred.female.sbs.m1000.2 <- rowMeans(post.female.sbs.m1000.2)
pred.female.sbs.m1e6.2 <- rowMeans(post.female.sbs.m1e6.2)
pred.female.sbs.m1 <- cbind(pred.female.sbs.m1.1,pred.female.sbs.m1.2)
pred.female.sbs.m1000 <- cbind(pred.female.sbs.m1000.1,pred.female.sbs.m1000.2)
pred.female.sbs.m1e6 <- cbind(pred.female.sbs.m1e6.1,pred.female.sbs.m1e6.2)
# Males
post.male.sbs.m1.1 <- sapply(post.male.sbs.m1,function(l) l[,1])
post.male.sbs.m1000.1 <- sapply(post.male.sbs.m1000,function(l) l[,1])
post.male.sbs.m1e6.1 <- sapply(post.male.sbs.m1e6,function(l) l[,1])
pred.male.sbs.m1.1 <- rowMeans(post.male.sbs.m1.1)
pred.male.sbs.m1000.1 <- rowMeans(post.male.sbs.m1000.1)
pred.male.sbs.m1e6.1 <- rowMeans(post.male.sbs.m1e6.1)
post.male.sbs.m1.2 <- sapply(post.male.sbs.m1,function(l) l[,2])
post.male.sbs.m1000.2 <- sapply(post.male.sbs.m1000,function(l) l[,2])
post.male.sbs.m1e6.2 <- sapply(post.male.sbs.m1e6,function(l) l[,2])
pred.male.sbs.m1.2 <- rowMeans(post.male.sbs.m1.2)
pred.male.sbs.m1000.2 <- rowMeans(post.male.sbs.m1000.2)
pred.male.sbs.m1e6.2 <- rowMeans(post.male.sbs.m1e6.2)
pred.male.sbs.m1 <- cbind(pred.male.sbs.m1.1,pred.male.sbs.m1.2)
pred.male.sbs.m1000 <- cbind(pred.male.sbs.m1000.1,pred.male.sbs.m1000.2)
pred.male.sbs.m1e6 <- cbind(pred.male.sbs.m1e6.1,pred.male.sbs.m1e6.2)

# Simulate several predicted replicates of the dataset and computes the Kalbflesih-Prentice estimatore
# Censoring times are simulated from the empirical distribution of the observed censoring times
tt <- seq(0,7000,1)
nW <- sum(melanoma$sex==0)
nM <- sum(melanoma$sex==1)
NREPS <- 100
# Women
rep_W.1 <- vector("list",length=NREPS)
rep_W.1000 <- vector("list",length=NREPS)
rep_W.1e6 <- vector("list",length=NREPS)
for(i in 1:NREPS){
  rW.1 <- sample.subdf(nW,tt,pred.female.sbs.m1)
  rW.1000 <- sample.subdf(nW,tt,pred.female.sbs.m1000)
  rW.1e6 <- sample.subdf(nW,tt,pred.female.sbs.m1e6)
  rcens.1 <- sample.subdf(nW,tt,matrix(cens[complete.cases(cens),],ncol=1))
  rcens.1000 <- sample.subdf(nW,tt,matrix(cens[complete.cases(cens),],ncol=1))
  rcens.1e6 <- sample.subdf(nW,tt,matrix(cens[complete.cases(cens),],ncol=1))
  I.cens.1 <- which(rW.1[,1]>rcens.1[,1])
  rW.1[I.cens.1,1] <- rcens.1[I.cens.1,1]
  rW.1[I.cens.1,2] <- 0
  I.cens.1000 <- which(rW.1000[,1]>rcens.1000[,1])
  rW.1000[I.cens.1000,1] <- rcens.1000[I.cens.1000,1]
  rW.1000[I.cens.1000,2] <- 0
  I.cens.1e6 <- which(rW.1e6[,1]>rcens.1e6[,1])
  rW.1e6[I.cens.1e6,1] <- rcens.1e6[I.cens.1e6,1]
  rW.1e6[I.cens.1e6,2] <- 0
  rep_W.1[[i]] <- freqest(tau=tau,times=rW.1[,1],delta=delta,types=rW.1[,2],k=k,
                        Wobs=matrix(c(1,0),ncol=2,nrow=nW,byrow = TRUE),
                        Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow = TRUE))$F.hat
  rep_W.1000[[i]] <- freqest(tau=tau,times=rW.1000[,1],delta=delta,types=rW.1000[,2],k=k,
                          Wobs=matrix(c(1,0),ncol=2,nrow=nW,byrow = TRUE),
                          Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow = TRUE))$F.hat
  rep_W.1e6[[i]] <- freqest(tau=tau,times=rW.1e6[,1],delta=delta,types=rW.1e6[,2],k=k,
                          Wobs=matrix(c(1,0),ncol=2,nrow=nW,byrow = TRUE),
                          Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow = TRUE))$F.hat
}
# Men
rep_M.1 <- vector("list",length=NREPS)
rep_M.1000 <- vector("list",length=NREPS)
rep_M.1e6 <- vector("list",length=NREPS)
for(i in 1:NREPS){
  rM.1 <- sample.subdf(nM,tt,pred.male.sbs.m1)
  rM.1000 <- sample.subdf(nM,tt,pred.male.sbs.m1000)
  rM.1e6 <- sample.subdf(nM,tt,pred.male.sbs.m1e6)
  rcens.1 <- sample.subdf(nM,tt,matrix(cens[complete.cases(cens),],ncol=1))
  rcens.1000 <- sample.subdf(nM,tt,matrix(cens[complete.cases(cens),],ncol=1))
  rcens.1e6 <- sample.subdf(nM,tt,matrix(cens[complete.cases(cens),],ncol=1))
  I.cens.1 <- which(rM.1[,1]>rcens.1[,1])
  rM.1[I.cens.1,1] <- rcens.1[I.cens.1,1]
  rM.1[I.cens.1,2] <- 0
  I.cens.1000 <- which(rM.1000[,1]>rcens.1000[,1])
  rM.1000[I.cens.1000,1] <- rcens.1000[I.cens.1000,1]
  rM.1000[I.cens.1000,2] <- 0
  I.cens.1e6 <- which(rM.1e6[,1]>rcens.1e6[,1])
  rM.1e6[I.cens.1e6,1] <- rcens.1e6[I.cens.1e6,1]
  rM.1e6[I.cens.1e6,2] <- 0
  rep_M.1[[i]] <- freqest(tau=tau,times=rM.1[,1],delta=delta,types=rM.1[,2],k=k,
                          Wobs=matrix(c(1,1),ncol=2,nrow=nM,byrow = TRUE),
                          Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow = TRUE))$F.hat
  rep_M.1000[[i]] <- freqest(tau=tau,times=rM.1000[,1],delta=delta,types=rM.1000[,2],k=k,
                             Wobs=matrix(c(1,1),ncol=2,nrow=nM,byrow = TRUE),
                             Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow = TRUE))$F.hat
  rep_M.1e6[[i]] <- freqest(tau=tau,times=rM.1e6[,1],delta=delta,types=rM.1e6[,2],k=k,
                            Wobs=matrix(c(1,1),ncol=2,nrow=nM,byrow = TRUE),
                            Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow = TRUE))$F.hat
}

# plots the results

pdf("d:/pred_checks_melanoma_M.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(5,5,3,3))
main <- expression(bold(log[10](m)==0))
plot(tt,freq_M[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_M.1[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_M[,1],type="s",col="red",lwd=2)

main <- expression(bold(log[10](m)==3))
plot(tt,freq_M[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_M.1000[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_M[,1],type="s",col="red",lwd=2)

main <- expression(bold(log[10](m)==6))
plot(tt,freq_M[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_M.1e6[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_M[,1],type="s",col="red",lwd=2)
par(mfrow=c(1,1))
dev.off()


pdf("d:/pred_checks_melanoma_W.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(5,5,3,3))
main <- expression(bold(log[10](m)==0))
plot(tt,freq_W[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_W.1[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_W[,1],type="s",col="red",lwd=2)

main <- expression(bold(log[10](m)==3))
plot(tt,freq_W[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_W.1000[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_W[,1],type="s",col="red",lwd=2)

main <- expression(bold(log[10](m)==6))
plot(tt,freq_W[,1],type="s",lwd=3,ylim=c(0,0.8),
     xlab="Days since surgery",ylab="Cumulative incidence",main=main,
     cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="red")
for(i in 1:NREPS) points(tt,rep_W.1e6[[i]][,1],type="s",col=rgb(0,0,0,alpha=0.15))
points(tt,freq_W[,1],type="s",col="red",lwd=2)
par(mfrow=c(1,1))
dev.off()

# 
# pdf("d:/pred_checks_other_M.pdf",height = 7,width = 21)
# par(mfcol=c(1,3))
# par(mar=c(5,5,3,3))
# main <- expression(bold(log[10](m)==0))
# plot(tt,freq_M[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_M.1[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.male.sbs.m1.2,type="s",col="red",lwd=2)
# 
# main <- expression(bold(log[10](m)==3))
# plot(tt,freq_M[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_M.1000[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.male.sbs.m1000.2,type="s",col="red",lwd=2)
# 
# main <- expression(bold(log[10](m)==6))
# plot(tt,freq_M[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_M.1e6[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.male.sbs.m1e6.2,type="s",col="red",lwd=2)
# par(mfrow=c(1,1))
# dev.off()
# 
# pdf("d:/pred_checks_other_W.pdf",height = 7,width = 21)
# par(mfcol=c(1,3))
# par(mar=c(5,5,3,3))
# main <- expression(bold(log[10](m)==0))
# plot(tt,freq_W[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_W.1[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.female.sbs.m1.2,type="s",col="red",lwd=2)
# 
# main <- expression(bold(log[10](m)==3))
# plot(tt,freq_W[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_W.1000[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.female.sbs.m1000.2,type="s",col="red",lwd=2)
# 
# main <- expression(bold(log[10](m)==6))
# plot(tt,freq_W[,2],type="s",lwd=3,ylim=c(0,0.8),
#      xlab="Days since surgery",ylab="Cumulative incidence",main=main,
#      cex.main=2.5,cex.lab=2,cex.axis=1.5)
# for(i in 1:NREPS) points(tt,rep_W.1e6[[i]][,2],type="s",col=rgb(0,0,0,alpha=0.15))
# points(tt,pred.female.sbs.m1e6.2,type="s",col="red",lwd=2)
# par(mfrow=c(1,1))
# dev.off()
