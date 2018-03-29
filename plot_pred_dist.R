
library(MCMCpack)
library(timereg)
data("melanoma")
source("C:/Users/andre/Dropbox/sub-beta-stacy/programs/functions.R")
set.seed(314271)

# load posterior samples 
in1 <- "C:/Users/andre/Dropbox/sub-beta-stacy/results/nonparam - weibull/melanoma_gender_m1.RData"
in2 <- "C:/Users/andre/Dropbox/sub-beta-stacy/results/nonparam - weibull/melanoma_gender_m1000.RData"
in3 <- "C:/Users/andre/Dropbox/sub-beta-stacy/results/nonparam - weibull/melanoma_gender_m1e6.RData"

load(in1); sbsw_m1    <- out_rescaled; 
load(in2); sbsw_m1000 <- out_rescaled; 
load(in3); sbsw_m1e6  <- out_rescaled;


# time discretization
delta <- 1                            
tau <- 7000
tt <- seq(0,tau,by=delta)

# build useful quantities
W <- cbind(1,melanoma[,c("sex")])    # cov matrix (sex=0 for females); includes column for the intercept term
times = melanoma$days                # event times
types = melanoma$status              # event types:      
types[types==2] = 0                  # 0. censored
types[types==3] = 2                  # 1. dead from melanoma
                                     # 2. dead from other causes

# Kalbfleish-Prentice estimators
kp_F <- freqest(tau=tau,times=times,delta=delta,types=types,Wobs=W,k=2,Wpred=matrix(c(1,0),ncol=2))
kp_M <- freqest(tau=tau,times=times,delta=delta,types=types,Wobs=W,k=2,Wpred=matrix(c(1,1),ncol=2))

# sample from the posterior distributions of the subdistribution beta-stacy model
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
    l <- postest(theta=trace[i,],
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
    s <- sbstacy(l$omega.post,l$F.post)
    return(s)
  })
  return(samp)
}
# Males
post.male.sbs.m1 <- samp.post.sbs(trace=sbsw_m1,
                               tau=tau,
                               Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow=TRUE),
                               Wobs=W,
                               times=times,
                               types=types,
                               delta=delta,
                               omega_m=1,
                               k=2)
post.male.sbs.m1000 <- samp.post.sbs(trace=sbsw_m1000,
                                     tau=tau,
                                     Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow=TRUE),
                                     Wobs=W,
                                     times=times,
                                     types=types,
                                     delta=delta,
                                     omega_m=1000,
                                     k=2)
post.male.sbs.m1e6 <- samp.post.sbs(trace=sbsw_m1e6,
                                     tau=tau,
                                     Wpred=matrix(c(1,1),ncol=2,nrow=1,byrow=TRUE),
                                     Wobs=W,
                                     times=times,
                                     types=types,
                                     delta=delta,
                                     omega_m=1e6,
                                     k=2)

# Females
post.female.sbs.m1 <- samp.post.sbs(trace=sbsw_m1,
                                 tau=tau,
                                 Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE),
                                 Wobs=W,
                                 times=times,
                                 types=types,
                                 delta=delta,
                                 omega_m=1,
                                 k=2)
post.female.sbs.m1000 <- samp.post.sbs(trace=sbsw_m1000,
                                       tau=tau,
                                       Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE),
                                       Wobs=W,
                                       times=times,
                                       types=types,
                                       delta=delta,
                                       omega_m=1000,
                                       k=2)
post.female.sbs.m1e6 <- samp.post.sbs(trace=sbsw_m1e6,
                                       tau=tau,
                                       Wpred=matrix(c(1,0),ncol=2,nrow=1,byrow=TRUE),
                                       Wobs=W,
                                       times=times,
                                       types=types,
                                       delta=delta,
                                       omega_m=1e6,
                                       k=2)


# Predictive distributions and posterior summaries
# Males
post.male.sbs.m1.1 <- sapply(post.male.sbs.m1,function(l) l[,1])
post.male.sbs.m1000.1 <- sapply(post.male.sbs.m1000,function(l) l[,1])
post.male.sbs.m1e6.1 <- sapply(post.male.sbs.m1e6,function(l) l[,1])
pred.male.sbs.m1.1 <- rowMeans(post.male.sbs.m1.1)
pred.male.sbs.m1000.1 <- rowMeans(post.male.sbs.m1000.1)
pred.male.sbs.m1e6.1 <- rowMeans(post.male.sbs.m1e6.1)
lcl95.male.sbs.m1.1 <- apply(post.male.sbs.m1.1,1,function(x) quantile(x,probs = 0.025))
lcl95.male.sbs.m1000.1 <- apply(post.male.sbs.m1000.1,1,function(x) quantile(x,probs = 0.025))
lcl95.male.sbs.m1e6.1 <- apply(post.male.sbs.m1e6.1,1,function(x) quantile(x,probs = 0.025))
ucl95.male.sbs.m1.1 <- apply(post.male.sbs.m1.1,1,function(x) quantile(x,probs = 0.975))
ucl95.male.sbs.m1000.1 <- apply(post.male.sbs.m1000.1,1,function(x) quantile(x,probs = 0.975))
ucl95.male.sbs.m1e6.1 <- apply(post.male.sbs.m1e6.1,1,function(x) quantile(x,probs = 0.975))
# Females
post.female.sbs.m1.1 <- sapply(post.female.sbs.m1,function(l) l[,1])
post.female.sbs.m1000.1 <- sapply(post.female.sbs.m1000,function(l) l[,1])
post.female.sbs.m1e6.1 <- sapply(post.female.sbs.m1e6,function(l) l[,1])
pred.female.sbs.m1.1 <- rowMeans(post.female.sbs.m1.1)
pred.female.sbs.m1000.1 <- rowMeans(post.female.sbs.m1000.1)
pred.female.sbs.m1e6.1 <- rowMeans(post.female.sbs.m1e6.1)
lcl95.female.sbs.m1.1 <- apply(post.female.sbs.m1.1,1,function(x) quantile(x,probs = 0.025))
lcl95.female.sbs.m1000.1 <- apply(post.female.sbs.m1000.1,1,function(x) quantile(x,probs = 0.025))
lcl95.female.sbs.m1e6.1 <- apply(post.female.sbs.m1e6.1,1,function(x) quantile(x,probs = 0.025))
ucl95.female.sbs.m1.1 <- apply(post.female.sbs.m1.1,1,function(x) quantile(x,probs = 0.975))
ucl95.female.sbs.m1000.1 <- apply(post.female.sbs.m1000.1,1,function(x) quantile(x,probs = 0.975))
ucl95.female.sbs.m1e6.1 <- apply(post.female.sbs.m1e6.1,1,function(x) quantile(x,probs = 0.975))

# Plot of the posterior summaries and predictive distributions

pdf("d:/pred_dist.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,kp_M$F.hat[,1],type="s",col="white",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
     xlab="Time since surgery",ylab="Cumulative incidence")

main <- expression(bold(log[10](m)==3))
plot(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
     xlab="Time since surgery",ylab="Cumulative incidence")

main <- expression(bold(log[10](m)==6))
plot(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lty=1,lwd=2,
     xlab="Time since surgery",ylab="Cumulative incidence")
par(mfrow=c(1,1))
dev.off()


pdf("d:/pred_dist_F.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,kp_F$F.hat[,1],type="s",col="white",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")

main <- expression(bold(log[10](m)==3))
plot(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")

main <- expression(bold(log[10](m)==6))
plot(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,ucl95.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3,lty=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lty=1,lwd=2,
       xlab="Time since surgery",ylab="Cumulative incidence")
par(mfrow=c(1,1))
dev.off()

#save.image("C:/Users/andre/Dropbox/sub-beta-stacy/results/plot_pred_dist.Rdata")
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/plot_pred_dist.Rdata")

