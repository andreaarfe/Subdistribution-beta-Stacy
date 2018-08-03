
library(MCMCpack)
library(timereg)
data("melanoma")
source("C:/Users/andre/Dropbox/sub-beta-stacy/programs/functions.R")
set.seed(314271)

load("C:/Users/andre/Dropbox/sub-beta-stacy/results/param - weibull/melanoma_param_gender_03082018.RData")
out_rescaled_param <- out_rescaled
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/plot_pred_dist.Rdata")


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

# sample from the posterior distributions of the subdistribution beta-stacy model
samp.post.param <- function(trace=NULL,
                          times=NULL,
                          delta=1,
                          k=NULL,
                          W=NULL,
                          theta=NULL,
                          omega_m=1,
                          log_u=FALSE,
                          dist="weibull"){
  samp <- lapply(1:nrow(trace), function(i){
    l <- getF0omega(theta=trace[i,],
                 times=times,
                 delta=delta,
                 k=k,
                 W=W,
                 omega_m=omega_m,
                 log_u=log_u,
                 dist=dist) 
    return(exp(l$F0))
  })
  return(samp)
}
# Males
post.male.param <- samp.post.param(trace=out_rescaled_param,
                                  W=matrix(c(1,1),ncol=2,nrow=length(tt),byrow=TRUE),
                                  times=tt,
                                  delta=delta,
                                  k=2)
# Females
post.female.param <- samp.post.param(trace=out_rescaled_param,
                                   W=matrix(c(1,0),ncol=2,nrow=length(tt),byrow=TRUE),
                                   times=tt,
                                   delta=delta,
                                   k=2)

# Predictive distributions and posterior summaries
# Males
post.male <- sapply(post.male.param,function(l) l[,1])
pred.male <- rowMeans(post.male)
lcl95.male <- apply(post.male,1,function(x) quantile(x,probs = 0.025))
ucl95.male <- apply(post.male,1,function(x) quantile(x,probs = 0.975))
# Females
post.female <- sapply(post.female.param,function(l) l[,1])
pred.female <- rowMeans(post.female)
lcl95.female <- apply(post.female,1,function(x) quantile(x,probs = 0.025))
ucl95.female <- apply(post.female,1,function(x) quantile(x,probs = 0.975))

pdf("C:/Users/andre/Dropbox/sub-beta-stacy/Figures/pred_dist_with_weibull_M.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,kp_M$F.hat[,1],type="s",col="white",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.male.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.male,lty=2,type="l",lwd=1)
points(tt,lcl95.male,lty=2,type="l",lwd=1)
points(tt,ucl95.male,lty=2,type="l",lwd=1)

main <- expression(bold(log[10](m)==3))
plot(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.male.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.male,lty=2,type="l",lwd=1)
points(tt,lcl95.male,lty=2,type="l",lwd=1)
points(tt,ucl95.male,lty=2,type="l",lwd=1)

main <- expression(bold(log[10](m)==6))
plot(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.male.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_M$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lty=1,lwd=2,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.male,lty=2,type="l",lwd=1)
points(tt,lcl95.male,lty=2,type="l",lwd=1)
points(tt,ucl95.male,lty=2,type="l",lwd=1)
par(mfrow=c(1,1))
dev.off()

pdf("C:/Users/andre/Dropbox/sub-beta-stacy/Figures/pred_dist_with_weibull_F.pdf",height = 3,width = 7)
par(mfcol=c(1,3))
par(mar=c(4,4,3,1))
main <- expression(bold(log[10](m)==0))
plot(tt,kp_F$F.hat[,1],type="s",col="white",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.female.sbs.m1.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.female,lty=2,type="l",lwd=1)
points(tt,lcl95.female,lty=2,type="l",lwd=1)
points(tt,ucl95.female,lty=2,type="l",lwd=1)
main <- expression(bold(log[10](m)==3))
plot(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.female.sbs.m1000.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=1,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.female,lty=2,type="l",lwd=1)
points(tt,lcl95.female,lty=2,type="l",lwd=1)
points(tt,ucl95.female,lty=2,type="l",lwd=1)
main <- expression(bold(log[10](m)==6))
plot(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lwd=2,lty=2,
     xlab="Days since surgery",ylab="Cumulative incidence",cex.main=1.5,cex.lab=1.2,cex.axis=1.2)
points(tt,pred.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,lcl95.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,ucl95.female.sbs.m1e6.1,type="s",ylim=c(0,1),col="gray60",lwd=3)
points(tt,kp_F$F.hat[,1],type="s",col="black",ylim=c(0,1),main=main,lty=1,lwd=2,
       xlab="Time since surgery",ylab="Cumulative incidence")
points(tt,pred.female,lty=2,type="l",lwd=1)
points(tt,lcl95.female,lty=2,type="l",lwd=1)
points(tt,ucl95.female,lty=2,type="l",lwd=1)
par(mfrow=c(1,1))
dev.off()


