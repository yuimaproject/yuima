### R code from vignette source '/Users/jago/Dropbox (VOICES)/yuima-book/chapter6.Rnw'

###################################################
### code chunk number 1: chapter6.Rnw:3-6
###################################################
options(width=60)
options(continue="  ")
require(yuima)


###################################################
### code chunk number 2: carma.brown
###################################################
carma.mod<-setCarma(p=3,q=1,loc.par="c0",Carma.var="y",Latent.var="X")
carma.mod


###################################################
### code chunk number 3: carma.brown.str
###################################################
str(carma.mod)


###################################################
### code chunk number 4: carma.brown.par
###################################################
par.carma<-list(a1=4,a2=4.75,a3=1.5,b0=1,b1=0.23,c0=0)
samp<-setSampling(Terminal=100, n=3000)  
set.seed(123)
carma <-simulate(carma.mod,
    true.parameter=par.carma, sampling=samp)


###################################################
### code chunk number 5: plot-carma
###################################################
plot(carma)


###################################################
### code chunk number 6: chapter6.Rnw:193-197
###################################################
pdf("figures/plot-carma.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(carma)
dev.off()


###################################################
### code chunk number 7: chapter6.Rnw:224-225 (eval = FALSE)
###################################################
## CarmaNoise(yuima, param, data=NULL)


###################################################
### code chunk number 8: carma.brown.qmle
###################################################
 fit <- qmle(carma, start=par.carma)
 fit


###################################################
### code chunk number 9: Carm21Comp0
###################################################
modCP<-setCarma(p=2,q=1,Carma.var="y",
 measure=list(intensity="Lamb",df=list("dnorm(z, mu, sig)")),
 measure.type="CP") 
true.parmCP <-list(a1=1.39631,a2=0.05029,b0=1,b1=2,
                  Lamb=1,mu=0,sig=1)


###################################################
### code chunk number 10: sim_Carm21Comp0
###################################################
samp.L<-setSampling(Terminal=200,n=4000)
set.seed(123)
simCP<-simulate(modCP,true.parameter=true.parmCP,sampling=samp.L)


###################################################
### code chunk number 11: plot-simCP
###################################################
plot(simCP,main="CP CARMA(2,1) model")


###################################################
### code chunk number 12: chapter6.Rnw:270-274
###################################################
pdf("figures/plot-simCP.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(simCP,main="CP CARMA(2,1) model")
dev.off()


###################################################
### code chunk number 13: plot_Carm21Comp0
###################################################
carmaoptCP <- qmle(simCP, start=true.parmCP, Est.Incr="Incr")
summary(carmaoptCP)


###################################################
### code chunk number 14: plot-carmaoptCP
###################################################
plot(carmaoptCP,ylab="Incr.",type="l",
 main="Compound Poisson with normal jump size")


###################################################
### code chunk number 15: chapter6.Rnw:291-295
###################################################
pdf("figures/plot-carmaoptCP.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(carmaoptCP,main="Compound Poisson with normal jump size",ylab="Incr.",type="l")
dev.off()


###################################################
### code chunk number 16: Carm21vg
###################################################
modVG<-setCarma(p=2,q=1,Carma.var="y",
     measure=list("rvgamma(z,lambda,alpha,beta,mu)"),
     measure.type="code") 
true.parmVG <-list(a1=1.39631, a2=0.05029, b0=1, b1=2,
                   lambda=1, alpha=1, beta=0, mu=0)


###################################################
### code chunk number 17: PlotCarm21vg
###################################################
set.seed(100)
simVG<-simulate(modVG, true.parameter=true.parmVG, 
 sampling=samp.L)
plot(simVG,main="VG CARMA(2,1) model")


###################################################
### code chunk number 18: chapter6.Rnw:321-325
###################################################
pdf("figures/plot-simVG.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(simVG,main="VG CARMA(2,1) model")
dev.off()


###################################################
### code chunk number 19: EstPlotCarm21vg
###################################################
carmaoptVG <- qmle(simVG, start=true.parmVG, Est.Incr="Incr")
summary(carmaoptVG)
plot(carmaoptVG,xlab="Time",
 main="Variance Gamma increments",ylab="Incr.",type="l")


###################################################
### code chunk number 20: chapter6.Rnw:339-343
###################################################
pdf("figures/plot-carmaoptVG.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(carmaoptVG,main="Variance Gamma increments",ylab="Incr.",xlab="Time",type="l")
dev.off()


###################################################
### code chunk number 21: Carm21Comp3
###################################################
modNIG<-setCarma(p=2,q=1,Carma.var="y",
   measure=list("rNIG(z,alpha,beta,delta1,mu)"),
   measure.type="code") 
IncMod<-setModel(drift="0",diffusion="0",jump.coeff="1",
  measure=list("rNIG(z,1,0,1,0)"),measure.type="code")
set.seed(100)
simLev<-simulate(IncMod,sampling=samp.L)
incrLevy<-diff(as.numeric(get.zoo.data(simLev)[[1]]))


###################################################
### code chunk number 22: plot-incrLevy
###################################################
plot(incrLevy,main="simulated noise increments",type="l")


###################################################
### code chunk number 23: chapter6.Rnw:366-370
###################################################
pdf("figures/plot-incrLevy.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(incrLevy,main="simulated noise increments",type="l")
dev.off()


###################################################
### code chunk number 24: Carm21sim3
###################################################
true.parmNIG <-list(a1=1.39631,a2=0.05029,b0=1,b1=2,
                  alpha=1,beta=0,delta1=1,mu=0)
simNIG<-simulate(modNIG,true.parameter=true.parmNIG,sampling=samp.L)


###################################################
### code chunk number 25: plot-simNIG
###################################################
plot(simNIG,main="NIG CARMA(2,1) model")


###################################################
### code chunk number 26: chapter6.Rnw:387-391
###################################################
pdf("figures/plot-simNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(simNIG,main="NIG CARMA(2,1) model")
dev.off()


###################################################
### code chunk number 27: plot_Carm21Comp3
###################################################
carmaoptNIG <- qmle(simNIG, start=true.parmNIG, Est.Incr="Incr")
summary(carmaoptNIG)


###################################################
### code chunk number 28: plot-carmaoptNIG
###################################################
plot(carmaoptNIG,main="Normal Inverse Gaussian",ylab="Incr.",type="l")


###################################################
### code chunk number 29: chapter6.Rnw:407-411
###################################################
pdf("figures/plot-carmaoptNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(carmaoptNIG,main="Normal Inverse Gaussian",ylab="Incr.",type="l")
dev.off()


###################################################
### code chunk number 30: Incrtime1Levy3a
###################################################
NIG.Inc<-as.numeric(coredata(carmaoptNIG@Incr.Lev))
NIG.freq<-frequency(carmaoptNIG@Incr.Lev)


###################################################
### code chunk number 31: Incrtime1Levy3b
###################################################
t.idx <- seq(from=1, to=length(NIG.Inc), by=NIG.freq)
Unitary.NIG.Inc<-diff(cumsum(NIG.Inc)[t.idx])


###################################################
### code chunk number 32: Incrtime1Levy3est
###################################################
library(GeneralizedHyperbolic)
FitInc.NIG.Lev<-nigFit(Unitary.NIG.Inc)
summary(FitInc.NIG.Lev, hessian = TRUE, hessianMethod = "tsHessian")


###################################################
### code chunk number 33: plot-fitNIG
###################################################
par(mfrow = c(1, 2))
plot(FitInc.NIG.Lev, which = 2:3,
       plotTitles = paste(c("Histogram of NIG ",
        "Log-Histogram of NIG ",
        "Q-Q Plot of NIG "), "Incr.", sep = ""))


###################################################
### code chunk number 34: chapter6.Rnw:447-456
###################################################
pdf("figures/plot-fitNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
par(mfrow = c(1, 2))
plot(FitInc.NIG.Lev, which = 2:3,
       plotTitles = paste(c("Histogram of NIG ",
                            "Log-Histogram of NIG ",
                            "Q-Q Plot of NIG "), "Incr.",
                          sep = ""))
dev.off()


###################################################
### code chunk number 35: chapter6.Rnw:468-473
###################################################
library(quantmod)
getSymbols("^VIX",  to="2016-12-31")
X <- VIX$VIX.Close
VIX.returns <- log(X)
plot(VIX.returns, main="VIX daily log-Returns")


###################################################
### code chunk number 36: chapter6.Rnw:476-480
###################################################
pdf("figures/plot-VIXret.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(VIX.returns, main="VIX daily log-Returns")
dev.off()


###################################################
### code chunk number 37: chapter6.Rnw:487-488
###################################################
acf(VIX.returns)


###################################################
### code chunk number 38: chapter6.Rnw:491-495
###################################################
pdf("figures/plot-acfVIX.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(acf(VIX.returns))
dev.off()


###################################################
### code chunk number 39: chapter6.Rnw:502-513
###################################################
library(TSA)
eacf(VIX.returns,ar.max = 3, ma.max = 4)

Delta <- 1/252
VIX.Data<-setData(VIX.returns,delta=Delta)
Normal.model<-setCarma(p=2, q=1,loc.par="mu")
Normal.CARMA<-setYuima(data=VIX.Data, model=Normal.model)
Normal.start <- list(a1=36,a2=56,b0=21,b1=1,mu=0)
Normal.est <- qmle(yuima=Normal.CARMA, start=Normal.start, 
 Est.Incr="Incr")
summary(Normal.est )


###################################################
### code chunk number 40: chapter6.Rnw:516-518
###################################################
inc <-Normal.est@Incr.Lev
shapiro.test(as.numeric(inc))


###################################################
### code chunk number 41: chapter6.Rnw:522-523
###################################################
plot(acf(as.numeric(inc)))


###################################################
### code chunk number 42: chapter6.Rnw:526-530
###################################################
pdf("figures/plot-acf2VIX.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(acf(as.numeric(inc)))
dev.off()


###################################################
### code chunk number 43: chapter6.Rnw:537-539
###################################################
Box.test(x=as.numeric(inc), lag = 10, type ="Ljung-Box")
Box.test(x=as.numeric(inc), lag = 10, type ="Box-Pierce")


###################################################
### code chunk number 44: chapter6.Rnw:543-567
###################################################
VG.model <- setCarma(p=2, q=1,loc.par="mu", 
 measure=list("rvgamma(z,lambda,alpha,beta,mu0)"), 
 measure.type="code")
NIG.model <- setCarma(p=2, q=1,loc.par="mu",
 measure=list(df=list("rNIG(z, alpha, beta, delta1, mu0)")), 
 measure.type="code")

VG.CARMA<-setYuima(data=VIX.Data, model=VG.model)
NIG.CARMA<-setYuima(data=VIX.Data, model=NIG.model)

VG.start <- list(a1=36,a2=56,b0=21,b1=1,mu=0,
 lambda=1,alpha=1,beta=0,mu0=0) 
NIG.start <- list(a1=36,a2=56,b0=21,b1=1,mu=0,
 alpha=2,beta=1,delta1=1,mu0=0)

fit.VG <- qmle(yuima=VG.CARMA,start=VG.start,
 Est.Incr="IncrPar",aggregation=FALSE)
fit.NIG <- qmle(yuima=NIG.CARMA,start=NIG.start,
  Est.Incr="IncrPar",aggregation=FALSE)
cf.VG <- coef(fit.VG )
cf.NIG <- coef(fit.NIG )

summary(fit.VG)
summary(fit.NIG)


###################################################
### code chunk number 45: chapter6.Rnw:571-586
###################################################
d.N <- function(u) log( 1+dnorm(u, mean=mean(inc), sd=sd(inc)) ) 
d.VG <- function(u) { 
 log(1+dvgamma(u, lambda=cf.VG["lambda"]*Delta, 
  alpha=cf.VG["alpha"], beta=cf.VG["beta"], mu=cf.VG["mu0"]*Delta)) 
}
d.NIG <- function(u) {
 log(1+dNIG(u,alpha=cf.NIG["alpha"], beta=cf.NIG["beta"],
   delta=cf.NIG["delta1"]*Delta, mu=cf.NIG["mu0"]*Delta))
}
d.Emp <- density(inc)
plot(d.Emp$x, log(1+d.Emp$y),type="l", 
 main="Rescaled log-densities")
curve(d.N, min(d.Emp$x), max(d.Emp$x), col="blue",add=TRUE, lty=3)
curve(d.VG, min(d.Emp$x), max(d.Emp$x), col="red",add=TRUE,lty=4)
curve(d.NIG, min(d.Emp$x), max(d.Emp$x), col="green",add=TRUE,lty=2)


###################################################
### code chunk number 46: chapter6.Rnw:589-596
###################################################
pdf("figures/plot-densVIX.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(d.Emp$x, log(1+d.Emp$y),type="l", main="Rescaled log-densities")
curve(d.N, min(d.Emp$x), max(d.Emp$x), col="blue",add=TRUE,lty=3, n=500)
curve(d.VG, min(d.Emp$x), max(d.Emp$x), col="red",add=TRUE,lty=4, n=500)
curve(d.NIG, min(d.Emp$x), max(d.Emp$x), col="green",add=TRUE,lty=2, n=500)
dev.off()


