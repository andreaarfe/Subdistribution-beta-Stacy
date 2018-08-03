
library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))


load("C:/Users/andre/Dropbox/sub-beta-stacy/results/simulations/out_exp_n100.Rdata")
out.100 <- out
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/simulations/out_exp_n500.Rdata")
out.500 <- out
load("C:/Users/andre/Dropbox/sub-beta-stacy/results/simulations/out_exp_n1000.Rdata")
out.1000 <- out
rm(list=c("out"))

KS.par <- unlist(lapply(out.100, function(l) l$KS.par))
KS.sbs.1 <- unlist(lapply(out.100, function(l) l$KS.sbs.1))
KS.sbs.1000 <- unlist(lapply(out.100, function(l) l$KS.sbs.1000))
KS.sbs.1e6 <- unlist(lapply(out.100, function(l) l$KS.sbs.1e6))
res <- data.frame(
  Model=rep(c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"),each=length(KS.par)),
  KS=c(KS.par,KS.sbs.1,KS.sbs.1000,KS.sbs.1e6)
)
res$Model <- factor(res$Model,levels=c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"))
res$n <- "n=100"
d <- res

KS.par <- unlist(lapply(out.500, function(l) l$KS.par))
KS.sbs.1 <- unlist(lapply(out.500, function(l) l$KS.sbs.1))
KS.sbs.1000 <- unlist(lapply(out.500, function(l) l$KS.sbs.1000))
KS.sbs.1e6 <- unlist(lapply(out.500, function(l) l$KS.sbs.1e6))
res <- data.frame(
  Model=rep(c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"),each=length(KS.par)),
  KS=c(KS.par,KS.sbs.1,KS.sbs.1000,KS.sbs.1e6)
)
res$Model <- factor(res$Model,levels=c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"))
res$n <- "n=500"
d <- rbind(d,res)

KS.par      <- unlist(lapply(out.1000, function(l) l$KS.par))
KS.sbs.1    <- unlist(lapply(out.1000, function(l) l$KS.sbs.1))
KS.sbs.1000 <- unlist(lapply(out.1000, function(l) l$KS.sbs.1000))
KS.sbs.1e6  <- unlist(lapply(out.1000, function(l) l$KS.sbs.1e6))
res <- data.frame(
  Model=rep(c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"),each=length(KS.par)),
  KS=c(KS.par,KS.sbs.1,KS.sbs.1000,KS.sbs.1e6)
)
res$Model <- factor(res$Model,levels=c("Parametric","SBS, m=1","SBS, m=1000","SBS, m=10^6"))
res$n <- "n=1000"
d <- rbind(d,res)
d$n <- factor(d$n,levels=c("n=100","n=500","n=1000"))


pdf("C:/Users/andre/Dropbox/sub-beta-stacy/Figures/sim_results.pdf",height = 4,width = 9)
d$Model <- factor(unclass(d$Model),labels=c("Parametric","m=10^0","m=10^3","m=10^6"))
bwplot(KS~Model|n,data=d,
       ylab="Kolmogorov-Smirnov distance", scales=list(cex=0.75),pch="|")
dev.off()


