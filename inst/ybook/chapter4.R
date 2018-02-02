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

## ------------------------------------------------------------------------
set.seed(123)
mu <- 0
sigma <- 1
lambda <- 10
samp <- setSampling(Terminal=10, n=1000)
mod10b <- setPoisson(intensity="lambda", df=list("dnorm(z,mu,sigma)"))
y10b <- simulate(mod10b,sampling=samp,
  true.par=list(lambda=lambda,mu=0.1, sigma=2))
y10b

## ----fig.keep='none'-----------------------------------------------------
BGmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
 measure.type="code", measure=list(df="rbgamma(z, delta.plus=1.4, 
 gamma.plus=0.3, delta.minus=2,
 gamma.minus=0.6)"))
n <- 1000
samp <- setSampling(Terminal=1, n=n)
BGyuima <- setYuima(model=BGmodel, sampling=samp)
set.seed(127)
for (i in 1:5) {
 result <- simulate(BGyuima)
 plot(result,xlim=c(0,1),ylim=c(-6,6),
  main="Paths of bilateral gamma process",col=i,par(new=T))
}

## ----plot-BGprocess,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-BGprocess.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(127)
for (i in 1:5) {
 result <- simulate(BGyuima)
 plot(result,xlim=c(0,1),ylim=c(-6,6),
  main="Paths of bilateral gamma process",col=i,par(new=T))
}
dev.off()

## ----fig.keep='none'-----------------------------------------------------
VGmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
 measure.type="code", measure=list(df="rbgamma(z, delta.minus=2,
 gamma.minus=0.6, delta.plus=2, gamma.plus=0.3)"))
VGyuima <- setYuima(model=VGmodel, sampling=samp)
set.seed(127)
for (i in 1:5) {
 result <- simulate(VGyuima)
 plot(result,xlim=c(0,1),ylim=c(-4,12),
  main="Paths of variance gamma process",col=i,par(new=T))
}

## ----plot-VGprocess,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-VGprocess.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(127)
for (i in 1:5) {
result <- simulate(VGyuima)
 plot(result,xlim=c(0,1),ylim=c(-4,12),
  main="Paths of variance gamma process",col=i,par(new=T))
}
dev.off()

## ----eval=FALSE----------------------------------------------------------
## Gmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
##  measure.type="code", measure=list(df="rgamma(z,
##  shape=0.7, scale=1)"))
## n <- 10000
## samp <- setSampling(Terminal=1, n=n)
## Gyuima <- setYuima(model=Gmodel, sampling=samp)
## set.seed(129)
## for (i in 1:5){
## result <- simulate(Gyuima)
##  plot(result,xlim=c(0,1),ylim=c(-0.1,1.2),
##   main="Paths of gamma process",col=i,par(new=T))
## }

## ----plot-Gprocess,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-Gprocess.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
Gmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
measure.type="code", measure=list(df="rgamma(z,
shape=0.7, scale=1)"))
n <- 10000
samp <- setSampling(Terminal=1, n=n)
Gyuima <- setYuima(model=Gmodel, sampling=samp)
set.seed(129)
for (i in 1:5){
 result <- simulate(Gyuima)
 plot(result,xlim=c(0,1),ylim=c(-0.1,1.2),
  main="Paths of gamma process",col=i,par(new=T))
}
dev.off()

## ----eval=FALSE----------------------------------------------------------
## n <- 5
## sampling <- setSampling(Terminal=1, n=n)
## Gmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
## measure.type="code", measure=list(df="rgamma(z,
## shape=0.7, scale=1)"))
## Gyuima <- setYuima(model=Gmodel, sampling=samp)
## simdata <- NULL
## set.seed(127)
## for (i in 1:3000){
## result <- simulate(Gyuima)
##  x1 <- result@data@original.data[n+1,1]
##  simdata <- c(simdata,as.numeric(x1))
## }
## hist(simdata, xlim=c(0,2), ylim=c(0,3), breaks=100, freq=FALSE,
##  main=expression(paste("Distribution of ", X[1],
##  " and Density of Gamma(0.7,1)")))
## curve(dgamma(x,0.7,1),add=TRUE,col="red")

## ----plot-dgamma,echo=FALSE,results='hide'-------------------------------
pdf("figures/plot-dgamma.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
n <- 5
samp <- setSampling(Terminal=1, n=n)
Gmodel <- setModel(drift="0", xinit="0", jump.coeff="1",
measure.type="code", measure=list(df="rgamma(z,
shape=0.7, scale=1)"))
Gyuima <- setYuima(model=Gmodel, sampling=samp)
simdata <- NULL
set.seed(127)
for (i in 1:3000){
 result <- simulate(Gyuima)
 x1 <- result@data@original.data[n+1,1]
 simdata <- c(simdata,as.numeric(x1))
}
hist(simdata, xlim=c(0,2), ylim=c(0,3), breaks=100, freq=FALSE,
 main=expression(paste("Distribution of ", X[1],
 " and Density of Gamma(0.7,1)")))
curve(dgamma(x,0.7,1),add=TRUE,col="red")
dev.off()

## ----fig.keep='none'-----------------------------------------------------
delta <- 1
gamma <- 2
set.seed(127)
x <- rIG(100000,delta,gamma)
hist(x,xlim=c(0,2),ylim=c(0,2),breaks=100,freq=FALSE)
curve(dIG(x,delta,gamma),add=TRUE,col="red", 
 from=min(x), to=max(x), n=500)
mean(x)
var(x)

## ----plot-dIG,echo=FALSE,results='hide'----------------------------------
pdf("figures/plot-dIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
hist(x,xlim=c(0,2),ylim=c(0,2),breaks=100,freq=FALSE)
curve(dIG(x,delta,gamma),add=TRUE,col="red", 
 from=min(x), to=max(x), n=500)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
IGmodel <- setModel(drift=0, xinit=0, jump.coeff=1, 
 measure.type="code", measure=list(df="rIG(z, delta=1, gamma=2)"))
n <- 1000
samp <- setSampling(Terminal=1, n=n)
IGyuima <- setYuima(model=IGmodel, sampling=samp)
set.seed(127)
for (i in 1:5){
 result <- simulate(IGyuima,xinit=0)
 plot(result, xlim=c(0,1), ylim=c(0,1),
  main="Paths of IG process (delta=1, gamma=2)",par(new=T),col=i)
}

## ----plot-dIGproc,echo=FALSE,results='hide'------------------------------
pdf("figures/plot-dIGproc.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(127)
for (i in 1:5){
 result <- simulate(IGyuima,xinit=0)
 plot(result, xlim=c(0,1), ylim=c(0,1),
  main="Paths of IG process (delta=1, gamma=2)",par(new=T),col=i)
}
dev.off()

## ----eval=FALSE----------------------------------------------------------
## n <- 5
## samp <- setSampling(Terminal=1, n=n)
## IGyuima <- setYuima(model=IGmodel, sampling=samp)
## IGsimdata <- NULL
## for (i in 1:3000){
##  result <- simulate(IGyuima)
##  x1 <- result@data@original.data[n+1,1]
##  IGsimdata <- c(IGsimdata,as.numeric(x1))
## }
## hist(IGsimdata,xlim=c(0,2), ylim=c(0,2), breaks=100, freq=FALSE,
##  main=expression(paste("Distribution of ",X[1],
##  " and Density of IG(1,2)")))
## curve(dIG(x,delta,gamma),add=TRUE,col="red",
##  from = 0.001, to = 5, n=500)

## ----plot-IGprocd,echo=FALSE,results='hide'------------------------------
pdf("figures/plot-IGprocd.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
n <- 5
samp <- setSampling(Terminal=1, n=n)
IGyuima <- setYuima(model=IGmodel, sampling=samp)
IGsimdata <- NULL
for (i in 1:3000){
 result <- simulate(IGyuima)
 x1 <- result@data@original.data[n+1,1]
 IGsimdata <- c(IGsimdata,as.numeric(x1))
}
hist(IGsimdata,xlim=c(0,2), ylim=c(0,2), breaks=100, freq=FALSE,
main=expression(paste("Distribution of ",X[1]," and Density of IG(1,2)")))
curve(dIG(x,delta,gamma),add=TRUE,col="red",
 from = 0.001, to = 5, n=500)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
rep <- 3000000
set.seed(129)
X1 <- rpts(rep,0.5,0.2,1)
hist(X1,xlim=c(0,3),ylim=c(0,3),breaks=100,
 main=expression(X[1]),probability=TRUE)
X05 <- rpts(rep,0.5,0.1,1)
X05.prime <- rpts(rep,0.5,0.1,1)
Xsum <- X05+X05.prime
summary(X1)
summary(Xsum)
ks.test(X1,Xsum)

## ----plot-X1pts,echo=FALSE,results='hide'--------------------------------
pdf("figures/plot-X1pts.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
hist(X1,xlim=c(0,3),ylim=c(0,3),breaks=100,main=expression(paste(X[1]," positive tempered stable distribution")),probability=TRUE)
dev.off()
rm(X1)
rm(Xsum)
rm(X05)
rm(X05.prime)

## ----fig.keep='none'-----------------------------------------------------
lambda <- 2
alpha <- 1.5
beta <- -0.7
mu <- 3
xinit <- 0
gamma <- sqrt(alpha^2-beta^2)
n <- 1000
T <- 1.8
VGPmodel <- setModel(drift=0, jump.coeff=1, measure.type="code",
 measure=list(df="rvgamma(z,lambda,alpha,beta,mu)"))
samp <- setSampling(Terminal=T, n=n)
VGPyuima <- setYuima(model=VGPmodel, sampling=samp)
# simulation
set.seed(127)
for (i in 1:7) {
 result <- simulate(VGPyuima, xinit=xinit,
 true.par=list(lambda=lambda,alpha=alpha,beta=beta,mu=mu))
plot(result,xlim=c(0,T),ylim=c(-5,6),col=i,
main="Paths of variance gamma process",par(new=T))
}

## ----plot-VGprocess2,echo=FALSE,results='hide'---------------------------
pdf("figures/plot-VGprocess2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(127)
for (i in 1:7) {
 result <- simulate(VGPyuima, xinit=xinit,
  true.par=list(lambda=lambda,alpha=alpha,beta=beta,mu=mu))
 plot(result,xlim=c(0,T),ylim=c(-5,6),col=i, main="Paths of variance gamma process",par(new=T))
}
dev.off()

## ----eval=FALSE----------------------------------------------------------
## n <- 5
## samp <- setSampling(Terminal=T, n=n)
## VGPyuima <- setYuima(model=VGPmodel, sampling=samp)
## VGPsimdata <- NULL
## for (i in 1:5000){
##  result <- simulate(VGPyuima, xinit=xinit,
##   true.par=list(lambda=lambda,alpha=alpha,beta=beta,mu=mu))
##   x1 <- result@data@original.data[n+1,1]
##   VGPsimdata <- c(VGPsimdata,as.numeric(x1[1]))
## }
## hist(VGPsimdata,xlim=c(-7,10),ylim=c(0,0.22),breaks=100,freq=FALSE,
##  main=expression(paste("Distribution of ",X[1.8],
##  " and Density of VG")))
## curve(dvgamma(x,lambda*T,alpha,beta,mu*T),add=TRUE,col="red")

## ----plot-VGPproc2,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-VGPproc2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
n <- 5
samp <- setSampling(Terminal=T, n=n)
VGPyuima <- setYuima(model=VGPmodel, sampling=samp)
VGPsimdata <- NULL
for (i in 1:5000){
 result <- simulate(VGPyuima, xinit=xinit,
  true.par=list(lambda=lambda,alpha=alpha,beta=beta,mu=mu))
  x1 <- result@data@original.data[n+1,1]
  VGPsimdata <- c(VGPsimdata,as.numeric(x1[1]))
}
hist(VGPsimdata,xlim=c(-7,10),ylim=c(0,0.22),breaks=100,freq=FALSE,
 main=expression(paste("Distribution of ",X[1.8], " and Density of VG")))
curve(dvgamma(x,lambda*T,alpha,beta,mu*T),add=TRUE,col="red")
dev.off()

## ----fig.keep='none'-----------------------------------------------------
delta <- 0.5
alpha <- 1.5
beta <- -0.7
mu <- 3
gamma <- sqrt(alpha^2-beta^2)
n <- 10000
T <- 1.8
set.seed(127)
normal.rn <- rnorm(n,0,1)
iv.rn <- rIG(n,delta*T,gamma)
z <- mu*T+beta*iv.rn+sqrt(iv.rn)*normal.rn
title <- expression(paste(NIGP[1.8],
 " built by subordination (green) and rNIG (white)"))
nig.rn <- rNIG(n,alpha,beta,delta*T,mu*T)
hist(z,xlim=c(-1,10),ylim=c(0,0.61),breaks=100, freq=FALSE,
 col="green", main=title, xlab=expression(X[1.8]) )
curve(dNIG(x,alpha,beta,delta*T,mu*T),add=TRUE,col="red")
par(new=T)
hist(nig.rn,xlim=c(-1,10),ylim=c(0,0.61),breaks=100,
 freq=FALSE, main="", xlab="")

## ----plot-NIGproc2,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-NIGproc2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
hist(z,xlim=c(-1,10),ylim=c(0,0.61),breaks=100, freq=FALSE,
col="green", main=title, xlab=expression(X[1.8]) )
curve(dNIG(x,alpha,beta,delta*T,mu*T),add=TRUE,col="red")
par(new=T)
hist(nig.rn,xlim=c(-1,10),ylim=c(0,0.61),breaks=100, freq=FALSE,
main="", xlab="")
dev.off()

## ----eval=FALSE----------------------------------------------------------
## delta1 <- 0.5
## alpha <- 1.5
## beta <- -0.7
## mu <- 3
## xinit <- 0
## gamma <- sqrt(alpha^2-beta^2)
## n <- 1000
## T <- 1.8
## NIG2model <- setModel(drift=0, jump.coeff=1, measure.type="code",
##  measure=list(df="rNIG(z,alpha,beta,delta1,mu)"))
## samp <- setSampling(Terminal=T, n=n)
## NIG2yuima <- setYuima(model=NIG2model, sampling=samp)
## set.seed(127)
## for (i in 1:10) {
##  result <- simulate(NIG2yuima, xinit=xinit,
##   true.par=list(delta1=delta1, alpha=alpha, beta=beta,
##   mu=mu, gamma=gamma))
## plot(result,xlim=c(0,T),ylim=c(-1,10),col=i,
##  main="Paths of NIG process",par(new=T))
## }

## ----plot-NIGproc3,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-NIGproc3.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
delta1 <- 0.5
alpha <- 1.5
beta <- -0.7
mu <- 3
xinit <- 0
gamma <- sqrt(alpha^2-beta^2)
n <- 1000
T <- 1.8
NIG2model <- setModel(drift=0, jump.coeff=1, measure.type="code",
measure=list(df="rNIG(z,alpha,beta,delta1,mu)"))
samp <- setSampling(Terminal=T, n=n)
NIG2yuima <- setYuima(model=NIG2model, sampling=samp)
set.seed(127)
for (i in 1:10) {
result <- simulate(NIG2yuima, xinit=xinit,
true.par=list(delta1=delta1, alpha=alpha, beta=beta,
mu=mu, gamma=gamma))
plot(result,xlim=c(0,T),ylim=c(-1,10),col=i,
main="Paths of NIG process",par(new=T))
}
dev.off()

## ----fig.keep='none'-----------------------------------------------------
n <- 5
samp <- setSampling(Terminal=T, n=n)
NIG2yuima <- setYuima(model=NIG2model, sampling=samp)
NIG2data <- NULL
for (i in 1:3000){
 result <- simulate(NIG2yuima, xinit=xinit,
  true.par=list(delta1=delta1, alpha=alpha, beta=beta,
  mu=mu, gamma=gamma))
 x1 <- result@data@original.data[n+1,1]
 NIG2data <- c(NIG2data,as.numeric(x1[1]))
}
hist(NIG2data,xlim=c(2,8),ylim=c(0,0.8),breaks=100, freq=FALSE,
 main=expression(paste("Distribution of ",X[1.8],
 " and Density of NIG")))
curve(dNIG(x,alpha,beta,delta*T,mu*T),add=TRUE,col="red")

## ----plot-NIGproc4,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-NIGproc4.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
hist(NIG2data,xlim=c(2,8),ylim=c(0,0.8),breaks=100, freq=FALSE,
main=expression(paste("Distribution of ",X[1.8], " and Density of NIG")))
curve(dNIG(x,alpha,beta,delta*T,mu*T),add=TRUE,col="red")
dev.off()

## ----fig.keep='none'-----------------------------------------------------
nrep <- 100000
alpha <- 0.5
delta <- 0.2
gamma <- 1
beta <- 1
mu <- -0.7
Lambda <- matrix(1,1,1)
t <- 1.5
par(mfrow=c(2,2))
set.seed(127)
x <- rnts(nrep,alpha,delta*t,gamma,beta,mu*t,Lambda)
s <- rpts(nrep,alpha,delta*t,gamma)
w <- rnorm(nrep,0,1)
y <- rep(mu*t,nrep) + beta*s + sqrt(s)*w
hist(x,xlim=c(-3,3),ylim=c(0,1.2),breaks=200,
 main=expression(X[t]),probability=TRUE)
hist(y,xlim=c(-3,3),ylim=c(0,1.2),breaks=200,
 main=expression(Y[t]),probability=TRUE,col="red")
## experiment by convolution
nrep <- 3000000
Xt <- rnts(nrep,alpha,delta*t,gamma,beta,mu*t,Lambda)
X05 <- rnts(nrep,alpha,delta*t/2,gamma,beta,mu*t/2,Lambda)
X05.prime <- rnts(nrep,alpha,delta*t/2,gamma,beta,mu*t/2,Lambda)
Xsum <- X05+X05.prime
hist(Xt,xlim=c(-3,3),ylim=c(0,1.2),breaks=300,
 main=expression(X[t]),probability=TRUE)
hist(Xsum,xlim=c(-3,3),ylim=c(0,1.2),breaks=300,
 main=expression(paste(X[t/2]+X[t/2],"'")),
 probability=TRUE,col="red")
ks.test(Xt,Xsum)

## ----plot-NTSPproc,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-NTSPproc.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
par(mfrow=c(2,2))
hist(x,xlim=c(-3,3),ylim=c(0,1.2),breaks=200,
main=expression(X[t]),probability=TRUE)
hist(y,xlim=c(-3,3),ylim=c(0,1.2),breaks=200,
main=expression(Y[t]),probability=TRUE,col="red")
hist(Xt,xlim=c(-3,3),ylim=c(0,1.2),breaks=300,
main=expression(X[t]),probability=TRUE)
hist(Xsum,xlim=c(-3,3),ylim=c(0,1.2),breaks=300,
main=expression(paste(X[t/2]+X[t/2],"'")),
probability=TRUE,col="red")
dev.off()
rm(Xt)
rm(X05)
rm(X05.prime)
rm(Xsum)

## ----eval=FALSE----------------------------------------------------------
## alpha <- 0.5
## beta <- -0.4
## sigma <- 0.7
## gamma <- 0.5
## n <- 1000
## T <- 1.8
## ASmodel <- setModel(drift=0, jump.coeff=1, measure.type="code",
##  measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
## samp <- setSampling(Terminal=T, n=n)
## ASyuima <- setYuima(model=ASmodel, sampling=samp)
## set.seed(129)
## for (i in 1:10) {
##  result <- simulate(ASyuima, true.par=list(alpha=alpha,
##   beta=beta,sigma=sigma,gamma=gamma))
##  plot(result,xlim=c(0,T),ylim=c(-40,10),col=i,
##   main=expression(paste("Paths of stable process (",
##   alpha==0.5,",",beta==-0.4,")")),par(new=T))
##  }
## 
## #param2
## alpha <- 1
## beta <- -0.4
## sigma <- 0.7
## gamma <- 0.5
## AS2model <- setModel(drift=0, jump.coeff=1, measure.type="code",
##  measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
## AS2yuima <- setYuima(model=AS2model, sampling=samp)
## for (i in 1:10) {
##  result <- simulate(AS2yuima, true.par=list(alpha=alpha,
##   beta=beta,sigma=sigma,gamma=gamma))
##  plot(result,xlim=c(0,T),ylim=c(-5,5),col=i,
##  main=expression(paste("Paths of stable process (",
##  alpha==1,",",beta==-0.4,")")),par(new=T))
## }
## 
## #param3
## alpha <- 1
## beta <- 0.4
## sigma <- 0.7
## gamma <- 0.5
## AS3model <- setModel(drift=0, jump.coeff=1, measure.type="code",
##  measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
## AS3yuima <- setYuima(model=AS3model, sampling=samp)
## for (i in 1:10) {
##  result <- simulate(AS3yuima, true.par=list(alpha=alpha,
##   beta=beta,sigma=sigma,gamma=gamma))
## plot(result,xlim=c(0,T),ylim=c(-5,5),col=i,
##  main=expression(paste("Paths of stable process (",
##  alpha==1,",",beta==0.4,")")),par(new=T))
## }
## 
## #param4
## alpha <- 1.5
## beta <- 0.4
## sigma <- 0.7
## gamma <- 0.5
## AS4model <- setModel(drift=0, jump.coeff=1, measure.type="code",
##  measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
## AS4yuima <- setYuima(model=AS4model, sampling=samp)
## for (i in 1:10) {
##  result <- simulate(AS4yuima, true.par=list(alpha=alpha,
##   beta=beta, sigma=sigma,gamma=gamma))
##  plot(result,xlim=c(0,T),ylim=c(-3,5),col=i,
##   main=expression(paste("Paths of stable process (",
##   alpha==1.5,",",beta==0.4,")")),par(new=T))
## }

## ----plot-ASproc,echo=FALSE,results='hide'-------------------------------
pdf("figures/plot-ASproc1.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
alpha <- 0.5
beta <- -0.4
sigma <- 0.7
gamma <- 0.5
n <- 1000
T <- 1.8
ASmodel <- setModel(drift=0, jump.coeff=1, measure.type="code",
measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
samp <- setSampling(Terminal=T, n=n)
ASyuima <- setYuima(model=ASmodel, sampling=samp)
set.seed(129)
for (i in 1:10) {
result <- simulate(ASyuima, true.par=list(alpha=alpha,
beta=beta,sigma=sigma,gamma=gamma))
plot(result,xlim=c(0,T),ylim=c(-40,10),col=i,
main=expression(paste("Paths of stable process (",
alpha==0.5,",",beta==-0.4,")")),par(new=T))
}
dev.off()
pdf("figures/plot-ASproc2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
#param2
alpha <- 1
beta <- -0.4
sigma <- 0.7
gamma <- 0.5
AS2model <- setModel(drift=0, jump.coeff=1, measure.type="code",
measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
AS2yuima <- setYuima(model=AS2model, sampling=samp)
for (i in 1:10) {
result <- simulate(AS2yuima, true.par=list(alpha=alpha,
beta=beta,sigma=sigma,gamma=gamma))
plot(result,xlim=c(0,T),ylim=c(-5,5),col=i,
main=expression(paste("Paths of stable process (",
alpha==1,",",beta==-0.4,")")),par(new=T))
}
dev.off()
pdf("figures/plot-ASproc3.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
#param3
alpha <- 1
beta <- 0.4
sigma <- 0.7
gamma <- 0.5
AS3model <- setModel(drift=0, jump.coeff=1, measure.type="code",
measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
AS3yuima <- setYuima(model=AS3model, sampling=samp)
for (i in 1:10) {
result <- simulate(AS3yuima, true.par=list(alpha=alpha,
beta=beta,sigma=sigma,gamma=gamma))
plot(result,xlim=c(0,T),ylim=c(-5,5),col=i,
main=expression(paste("Paths of stable process (",
alpha==1,",",beta==0.4,")")),par(new=T))
}
dev.off()
pdf("figures/plot-ASproc4.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
#param4
alpha <- 1.5
beta <- 0.4
sigma <- 0.7
gamma <- 0.5
AS4model <- setModel(drift=0, jump.coeff=1, measure.type="code",
measure=list(df="rstable(z,alpha,beta,sigma,gamma)"))
AS4yuima <- setYuima(model=AS4model, sampling=samp)
for (i in 1:10) {
result <- simulate(AS4yuima, true.par=list(alpha=alpha,beta=beta,
sigma=sigma,gamma=gamma))
plot(result,xlim=c(0,T),ylim=c(-3,5),col=i,
main=expression(paste("Paths of stable process (",
alpha==1.5,",",beta==0.4,")")),par(new=T))
}
dev.off()

## ----fig.keep='none'-----------------------------------------------------
modJump <- setModel(drift = c("-theta*x"), diffusion = "sigma",
 jump.coeff=c("gamma+x/sqrt(1+x^2)"),
 measure = list(intensity="lambda",df=list("dnorm(z, -3, 1)")),
 measure.type="CP", solve.variable="x")
modJump
samp <- setSampling(n=10000,Terminal=10)
set.seed(125)
X <- simulate(modJump, xinit=2, sampling=samp,
 true.par= list(theta=2, sigma=0.5,gamma=0.3,lambda=0.5))
plot(X)

## ----plot-modelJump,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-modelJump.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(X)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
x0 <- 2
a <- 0.1
c <- -1
model.ig <- setModel(drift="a*x", xinit=x0, jump.coeff=c, 
 measure.type="code", measure=list(df="rIG(z, delta0, gamma)"))
model.ig
sampling.ig <- setSampling(Terminal=10, n=10000)
yuima.ig <- setYuima(model=model.ig, sampling=sampling.ig)
set.seed(128)
result.ig <- simulate(yuima.ig,true.par=list(delta0=0.55,gamma=2))
plot(result.ig)

## ----plot-modelIG,echo=FALSE,results='hide'------------------------------
pdf("figures/plot-modelIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(128)
result.ig <- simulate(yuima.ig,true.par=list(delta0=0.55,gamma=2))
plot(result.ig)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
x0 <- 2
a <- 0.1
c <- -1
model.nig <- setModel(drift="a*x", xinit=x0, jump.coeff=c,
 measure.type="code",measure=list(df="rNIG(z, alpha,
 beta, delta0, mu)"))
sampling.nig <- setSampling(Terminal=10, n=10000)
yuima.nig <- setYuima(model=model.nig, sampling=sampling.ig)
set.seed(128)
result.nig <- simulate(yuima.nig,true.par=list(alpha=2, beta=0,
 delta0=0.55, mu=0))
plot(result.nig)

## ----plot-modelNIG,echo=FALSE,results='hide'-----------------------------
pdf("figures/plot-modelNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(128)
result.nig <- simulate(yuima.nig,true.par=list(alpha=2, beta=0,
delta0=0.55, mu=0))
plot(result.nig)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
x0 <- 2
a <- 0.1
c <- -1
Lambda <- matrix(1,1,1)
model.nig <- setModel(drift="a*x", xinit=x0, jump.coeff=c,
 measure.type="code",measure=list(df="rNIG(z, alpha,
 beta, delta0, mu, Lambda)"))
sampling.nig <- setSampling(Terminal=10, n=10000)
yuima.nig <- setYuima(model=model.nig, sampling=sampling.ig)
set.seed(128)
result.nig <- simulate(yuima.nig,true.par=list(alpha=2,
 beta=0, delta0=0.55, mu=0, Lambda=Lambda))
plot(result.nig)

## ----plot-modelNIG2,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-modelNIG2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(128)
result.nig <- simulate(yuima.nig,true.par=list(alpha=2,
beta=0, delta0=0.55, mu=0, Lambda=Lambda))
plot(result.nig)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
x0 <- c(2,3)
a1 <- function(t,x1,x2){ x1*cos(2*pi*t)-x2*sin(2*pi*t) }
a2 <- function(t,x1,x2){ x1*sin(2*pi*t)+x2*cos(2*pi*t) }
a <- c("a1(t,x1,x2)","a2(t,x1,x2)")
b <- matrix(c("t*x2","1","0","x1"),2,2)
c <- matrix(c("cos(2*pi*t)", "(5-t)*x1","sin(2*pi*t)",1),2,2)
alpha <- 2
beta <- c(0,0)
delta0 <- 0.55
mu <- c(0,0)
Lambda <- matrix(c(1,0,0,1),2,2)
model.mnig <- setModel(drift=a, xinit=x0, diffusion=b,
  jump.coeff=c, measure.type="code",
  measure=list(df="rNIG(z, alpha, beta, delta0, mu, Lambda)"),
  state.variable=c("x1","x2"),solve.variable=c("x1","x2") )
model.mnig
sampling.mnig <- setSampling(Terminal=1, n=10000)
yuima.mnig <- setYuima(model=model.mnig, sampling=sampling.mnig)
set.seed(128)
result.mnig <- simulate(yuima.mnig,true.par=list(alpha=alpha,
 beta=beta, delta0=delta0, mu=mu, Lambda=Lambda))
plot(result.mnig)

## ----plot-modelMNIG,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-modelMNIG.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
set.seed(128)
result.mnig <- simulate(yuima.mnig,true.par=list(alpha=alpha,
beta=beta, delta0=delta0, mu=mu, Lambda=Lambda))
plot(result.mnig)
dev.off()

## ----fig.keep='none', results='hide'-------------------------------------
mod5 <- setModel(drift = c("-theta*x"), diffusion = "sigma",
jump.coeff=c("gamma+x/sqrt(1+x^2)"),
measure = list(intensity="lambda",df=list("dnorm(z, 2, 0.1)")),
measure.type="CP", solve.variable="x")
theta <- 2
sigma <- 0.5
gamma <- 0.3
lambda <- 2.5
T <- 10
N <- 10000
delta <- T/N
h <- T/N
true <- list(theta=theta, sigma=sigma,gamma=gamma,lambda=lambda)
set.seed(125)
X <- simulate(mod5, true.p=true,xinit=2,
sampling=setSampling(n=N,Terminal=T))
plot(X)
r <- h^0.4
est.qmle <- qmle(yuima=X, start=true,
 lower=list(theta=1,sigma=0,gamma=0.1,lambda=0.1), 
 upper=list(theta=3,sigma=2,gamma=0.8,lambda=20), method="L-BFGS-B",
 threshold=r)
unlist(true)
summary(est.qmle)

## ----echo=FALSE----------------------------------------------------------
unlist(true)
writeLines(strwrap(capture.output(summary(est.qmle)),width=60))

## ----plot-modelSDEJ,echo=FALSE,results='hide'----------------------------
pdf("figures/plot-modelSDEJ.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(X)
dev.off()

## ------------------------------------------------------------------------
est.qmle1 <- qmle(yuima=X, start=true,
 lower=list(theta=1,sigma=0,gamma=0.1,lambda=0.1),
 upper=list(theta=3,sigma=2,gamma=0.8,lambda=20), method="L-BFGS-B",
 threshold=2) # too large
coef(est.qmle1)
est.qmle2 <- qmle(yuima=X, start=true,
 lower=list(theta=1,sigma=0,gamma=0.1,lambda=0.1),
 upper=list(theta=3,sigma=2,gamma=10,lambda=1000), method="L-BFGS-B",
 threshold=0.03) ## too low
coef(est.qmle2)

## ----fig.keep='none',message=FALSE---------------------------------------
require(quantmod)
getSymbols("ENI.MI",to="2016-12-31")
S <- ENI.MI$ENI.MI.Adjusted
Z <- na.omit(diff(log(S)))
Dt <- 1/252
# geometric Brownian motion estimation
model1 <- setModel(drift="mu*x", diff="sigma*x")
gBm <- setYuima(model=model1, data=setData(S,delta=Dt))
gBm.fit <- qmle(gBm, start=list(mu=0,sigma=1),method="BFGS")
gBm.cf <- coef(gBm.fit)
zMin <- min(Z)
zMax <- max(Z)
# Gaussian-Levy estimation
model3 <- setPoisson( df="dnorm(z,mu,sigma)")
Norm <- setYuima(model=model3, data=setData(cumsum(Z),delta=Dt))
Norm.fit <- qmle(Norm,start=list(mu=1, sigma=1),
 lower=list(mu=1e-7,sigma=0.01),method="L-BFGS-B")
Norm.cf <- coef(Norm.fit)
# NIG-Levy estimation
model2 <- setPoisson( df="dNIG(z,alpha,beta,delta1,mu)")
NIG <- setYuima(model=model2, data=setData(cumsum(Z),delta=Dt))
NIG.fit <- qmle(NIG,start=list(alpha=10, beta=1, delta1=1,mu=1),
 lower=list(alpha=1,beta=-2, delta1=0.001,mu=0.0001),
  method="L-BFGS-B")
NIG.cf <- coef(NIG.fit)
myfgBm <- function(u) 
 dnorm(u, mean=gBm.cf["mu"], sd=gBm.cf["sigma"])
myfNorm <- function(u) 
 dnorm(u, mean=Norm.cf["mu"],sd=Norm.cf["sigma"])
myfNIG <- function(u) 
 dNIG(u, alpha=NIG.cf["alpha"],beta=NIG.cf["beta"],
 delta=NIG.cf["delta1"], mu=NIG.cf["mu"])
plot(density(Z,na.rm=TRUE),main="Gaussian versus NIG")
curve(myfgBm, zMin, zMax, add=TRUE, lty=2)
curve(myfNorm, zMin, zMax, col="red", add=TRUE, lty=4)
curve(myfNIG, zMin, zMax, col="blue", add=TRUE,lty=3)

## ----plot-modelExpLevy,echo=FALSE,results='hide'-------------------------
pdf("figures/plot-modelExpLevy.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(density(Z,na.rm=TRUE),main="Gaussian versus NIG")
curve(myfgBm, zMin, zMax, add=TRUE, lty=2)
curve(myfNorm, zMin, zMax, col="red", add=TRUE, lty=4)
curve(myfNIG, zMin, zMax, col="blue", add=TRUE,lty=3)
dev.off()

## ------------------------------------------------------------------------
AIC(gBm.fit)
AIC(Norm.fit)
AIC(NIG.fit)

