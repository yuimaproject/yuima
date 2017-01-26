## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
tidy=FALSE,
width.cutoff = 60,
strip.white=TRUE,
warning=FALSE
)

## ----include=FALSE-------------------------------------------------------
options(width=55)
options(continue="  ")
require(yuima)

## ----carma.brown,echo=TRUE,eval=TRUE-------------------------------------
carma.mod<-setCarma(p=3,q=1,loc.par="c0",Carma.var="y",Latent.var="X")
carma.mod

## ----carma.brown.str,results='hide'--------------------------------------
str(carma.mod)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(carma.mod)),width=60))

## ----carma.brown.par,echo=TRUE,eval=TRUE---------------------------------
par.carma<-list(a1=4,a2=4.75,a3=1.5,b0=1,b1=0.23,c0=0)
samp<-setSampling(Terminal=100, n=3000)  
set.seed(123)
carma <-simulate(carma.mod,
    true.parameter=par.carma, sampling=samp)

## ----plot-carma,echo=TRUE,fig.keep='none',results='hide'-----------------
plot(carma)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-carma.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(carma)
dev.off()

## ----eval=FALSE----------------------------------------------------------
## CarmaNoise(yuima, param, data=NULL)

## ----carma.brown.qmle,echo=TRUE,eval=TRUE--------------------------------
 fit <- qmle(carma, start=par.carma)
 fit

## ----Carm21Comp0, echo=TRUE----------------------------------------------
modCP<-setCarma(p=2,q=1,Carma.var="y",
 measure=list(intensity="Lamb",df=list("dnorm(z, mu, sig)")),
 measure.type="CP") 
true.parmCP <-list(a1=1.39631,a2=0.05029,b0=1,b1=2,
                  Lamb=1,mu=0,sig=1)

## ----sim_Carm21Comp0, echo=TRUE------------------------------------------
samp.L<-setSampling(Terminal=200,n=4000)
set.seed(123)
simCP<-simulate(modCP,true.parameter=true.parmCP,sampling=samp.L)

## ----plot-simCP,echo=TRUE,fig.keep='none',results='hide'-----------------
plot(simCP,main="CP CARMA(2,1) model")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-simCP.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(simCP,main="CP CARMA(2,1) model")
dev.off()

## ----plot_Carm21Comp0,echo=TRUE,fig.width=14,fig.height=7----------------
carmaoptCP <- qmle(simCP, start=true.parmCP, Est.Incr="Incr")
summary(carmaoptCP)

## ----plot-carmaoptCP,echo=TRUE,fig.keep='none',results='hide'------------
plot(carmaoptCP,ylab="Incr.",type="l",
 main="Compound Poisson with normal jump size")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-carmaoptCP.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(carmaoptCP,main="Compound Poisson with normal jump size",ylab="Incr.",type="l")
dev.off()

## ----Carm21vg, echo=TRUE-------------------------------------------------
modVG<-setCarma(p=2,q=1,Carma.var="y",
     measure=list("rvgamma(z,lambda,alpha,beta,mu)"),
     measure.type="code") 
true.parmVG <-list(a1=1.39631, a2=0.05029, b0=1, b1=2,
                   lambda=1, alpha=1, beta=0, mu=0)

## ----PlotCarm21vg, echo=TRUE,fig.width=14,fig.height=7,fig.keep='none',results='hide'----
set.seed(100)
simVG<-simulate(modVG, true.parameter=true.parmVG, 
 sampling=samp.L)
plot(simVG,main="VG CARMA(2,1) model")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-simVG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(simVG,main="VG CARMA(2,1) model")
dev.off()

## ----EstPlotCarm21vg,echo=TRUE, fig.keep='none',results='hide'-----------
carmaoptVG <- qmle(simVG, start=true.parmVG, Est.Incr="Incr")
summary(carmaoptVG)
plot(carmaoptVG,xlab="Time",
 main="Variance Gamma increments",ylab="Incr.",type="l")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-carmaoptVG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(carmaoptVG,main="Variance Gamma increments",ylab="Incr.",xlab="Time",type="l")
dev.off()

## ----Carm21Comp3, echo=TRUE----------------------------------------------
modNIG<-setCarma(p=2,q=1,Carma.var="y",
   measure=list("rNIG(z,alpha,beta,delta1,mu)"),
   measure.type="code") 
IncMod<-setModel(drift="0",diffusion="0",jump.coeff="1",
  measure=list("rNIG(z,1,0,1,0)"),measure.type="code")
set.seed(100)
simLev<-simulate(IncMod,sampling=samp.L)
incrLevy<-diff(as.numeric(get.zoo.data(simLev)[[1]]))

## ----plot-incrLevy,echo=TRUE,fig.keep='none',results='hide'--------------
plot(incrLevy,main="simulated noise increments",type="l")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-incrLevy.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(incrLevy,main="simulated noise increments",type="l")
dev.off()

## ----Carm21sim3, echo=TRUE-----------------------------------------------
true.parmNIG <-list(a1=1.39631,a2=0.05029,b0=1,b1=2,
                  alpha=1,beta=0,delta1=1,mu=0)
simNIG<-simulate(modNIG,true.parameter=true.parmNIG,sampling=samp.L)

## ----plot-simNIG,echo=TRUE,fig.keep='none',results='hide'----------------
plot(simNIG,main="NIG CARMA(2,1) model")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-simNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(simNIG,main="NIG CARMA(2,1) model")
dev.off()

## ----plot_Carm21Comp3 ,echo=TRUE,fig.width=14,fig.height=7---------------
carmaoptNIG <- qmle(simNIG, start=true.parmNIG, Est.Incr="Incr")
summary(carmaoptNIG)

## ----plot-carmaoptNIG,echo=TRUE,fig.keep='none',results='hide'-----------
plot(carmaoptNIG,main="Normal Inverse Gaussian",ylab="Incr.",type="l")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-carmaoptNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(carmaoptNIG,main="Normal Inverse Gaussian",ylab="Incr.",type="l")
dev.off()

## ----Incrtime1Levy3a, echo=TRUE------------------------------------------
NIG.Inc<-as.numeric(coredata(carmaoptNIG@Incr.Lev))
NIG.freq<-frequency(carmaoptNIG@Incr.Lev)

## ----Incrtime1Levy3b, echo=TRUE------------------------------------------
t.idx <- seq(from=1, to=length(NIG.Inc), by=NIG.freq)
Unitary.NIG.Inc<-diff(cumsum(NIG.Inc)[t.idx])

## ----Incrtime1Levy3est, echo=TRUE----------------------------------------
library(GeneralizedHyperbolic)
FitInc.NIG.Lev<-nigFit(Unitary.NIG.Inc)
summary(FitInc.NIG.Lev, hessian = TRUE, hessianMethod = "tsHessian")

## ----plot-fitNIG,echo=TRUE,fig.keep='none',results='hide'----------------
par(mfrow = c(1, 2))
plot(FitInc.NIG.Lev, which = 2:3,
       plotTitles = paste(c("Histogram of NIG ",
        "Log-Histogram of NIG ",
        "Q-Q Plot of NIG "), "Incr.", sep = ""))

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-fitNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
par(mfrow = c(1, 2))
plot(FitInc.NIG.Lev, which = 2:3,
       plotTitles = paste(c("Histogram of NIG ",
                            "Log-Histogram of NIG ",
                            "Q-Q Plot of NIG "), "Incr.",
                          sep = ""))
dev.off()

## ----message=FALSE,fig.keep='none'---------------------------------------
library(quantmod)
getSymbols("^VIX",  to="2016-12-31")
X <- VIX$VIX.Close
VIX.returns <- log(X)
plot(VIX.returns, main="VIX daily log-Returns")

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-VIXret.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(VIX.returns, main="VIX daily log-Returns")
dev.off()

## ----fig.keep='none'-----------------------------------------------------
acf(VIX.returns)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-acfVIX.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(acf(VIX.returns))
dev.off()

## ----message=FALSE,fig.keep='none'---------------------------------------
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

## ------------------------------------------------------------------------
inc <-Normal.est@Incr.Lev
shapiro.test(as.numeric(inc))

## ----fig.keep='none'-----------------------------------------------------
plot(acf(as.numeric(inc)))

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-acf2VIX.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(acf(as.numeric(inc)))
dev.off()

## ------------------------------------------------------------------------
Box.test(x=as.numeric(inc), lag = 10, type ="Ljung-Box")
Box.test(x=as.numeric(inc), lag = 10, type ="Box-Pierce")

## ------------------------------------------------------------------------
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

## ----fig.keep='none'-----------------------------------------------------
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

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-densVIX.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(d.Emp$x, log(1+d.Emp$y),type="l", main="Rescaled log-densities")
curve(d.N, min(d.Emp$x), max(d.Emp$x), col="blue",add=TRUE,lty=3, n=500)
curve(d.VG, min(d.Emp$x), max(d.Emp$x), col="red",add=TRUE,lty=4, n=500)
curve(d.NIG, min(d.Emp$x), max(d.Emp$x), col="green",add=TRUE,lty=2, n=500)
dev.off()

