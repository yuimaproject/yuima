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

## ----results='hide',tidy=FALSE,width.cutoff = 60-------------------------
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)", 
 xinit="rnorm(1)")
str(mod1)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(mod1)),width=60))

## ----plot-mod1diff,echo=TRUE,fig.keep='none',results='hide'--------------
set.seed(123)
x1 <- simulate(mod1)
x2 <- simulate(mod1)
par(mfrow=c(1,2))
plot(x1)
plot(x2)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1diff.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
set.seed(123)
x1 <- simulate(mod1)
x2 <- simulate(mod1)
par(mfrow=c(1,2))
plot(x1)
plot(x2)
dev.off()

## ----results='hide'------------------------------------------------------
mod2 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)", 
 xinit="rnorm(1, mean=mu)")
mod2
str(mod2)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(mod2)),width=60))

## ----eval=FALSE----------------------------------------------------------
## x <- simulate(mod2, true.par=list(mu=1))

## ----plot-mod1diff2,echo=TRUE,fig.keep='none',results='hide'-------------
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)")
set.seed(123)
x1 <- simulate(mod1, xinit=1)
x2 <- simulate(mod1, xinit=expression(rnorm(1)))
x3 <- simulate(mod2, xinit=3)
par(mfrow=c(1,3))
plot(x1, main="mod1, xinit=1")
plot(x2, main="mod1, xinit=expression(rnorm(1))")
plot(x3, main="mod2, xinit=3")
par(mfrow=c(1,1))

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1diff2.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)")
set.seed(123)
x1 <- simulate(mod1, xinit=1)
x2 <- simulate(mod1, xinit=expression(rnorm(1)))
x3 <- simulate(mod2, xinit=3)
par(mfrow=c(1,3))
plot(x1, main="mod1, xinit=1")
plot(x2, main="mod1, xinit=expression(rnorm(1))")
plot(x3, main="mod2, xinit=3")
dev.off()

## ------------------------------------------------------------------------
ou <- setModel(drift="-theta*x", diffusion=1)

## ------------------------------------------------------------------------
gBm <- setModel(drift="mu*x", diffusion="sigma*x")

## ------------------------------------------------------------------------
vasicek <- setModel(drift="theta1-theta2*x", diffusion="theta3")

## ------------------------------------------------------------------------
cev <- setModel(drift="mu*x", diffusion="sigma*x^gamma")

## ------------------------------------------------------------------------
cir <- setModel(drift="theta1-theta2*x", diffusion="theta3*sqrt(x)")

## ------------------------------------------------------------------------
ckls <- setModel(drift="theta1-theta2*x", diffusion="theta3*x^theta4")

## ------------------------------------------------------------------------
hyper1 <- setModel( diff="sigma", 
 drift="(sigma^2/2)*(beta-alpha*((x-mu)/(sqrt(delta^2+(x-mu)^2))))")

## ------------------------------------------------------------------------
hyper1
str(hyper1@parameter)

## ------------------------------------------------------------------------
hyper2 <- setModel(drift="0", 
 diffusion = "sigma*exp(0.5*alpha*sqrt(delta^2+(x-mu)^2)-
  beta*(x-mu))")

## ------------------------------------------------------------------------
hyper2
str(hyper2@parameter)

## ------------------------------------------------------------------------
set.seed(123)
modA <- setModel(drift="-0.3*x", diffusion=1)
modB <- setModel(drift="0.3*x", diffusion=1)

## Set the model in an `yuima' object with a sampling scheme. 
Terminal <- 1
n <- 500
mod.sampling <- setSampling(Terminal=Terminal, n=n)
yuima1 <- setYuima(model=modA, sampling=mod.sampling)
yuima2 <- setYuima(model=modB, sampling=mod.sampling)

##use original increment
delta <- Terminal/n
my.dW <- matrix( rnorm(n , 0, sqrt(delta)), nrow=1, ncol=n)

## Solve SDEs using Euler-Maruyama method.
y1 <- simulate(yuima1, xinit=1, increment.W=my.dW)
y2 <- simulate(yuima2, xinit=1, increment.W=my.dW)

## ----plot-modAB,echo=TRUE,fig.keep='none',results='hide'-----------------
plot(y1)
lines(get.zoo.data(y2)[[1]], col="red",lty=3)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-modAB.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(y1)
lines(get.zoo.data(y2)[[1]], col="red",lty=3)
dev.off()

## ----echo=TRUE,results='hide'--------------------------------------------
sol <- c("x1","x2") # variable for numerical solution
a <- c("-3*x1","-x1-2*x2")   # drift vector 
b <- matrix(c("1","x1","0","3","x2","0"),2,3)  #  diffusion matrix
mod3 <- setModel(drift = a, diffusion = b, solve.variable = sol)

## ----sim-mod3,echo=TRUE,fig.keep='none',results='hide'-------------------
set.seed(123)
X <- simulate(mod3)
plot(X, plot.type="single",lty=1:2)

## ----plot-mod3,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-mod3.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X, plot.type="single",lty=1:2)
dev.off()

## ----echo=TRUE,results='hide'--------------------------------------------
mu <- 0.1
sig <- 0.2
rho <- -0.7

g <- function(t) {0.4 + (0.1 + 0.2*t)* exp(-2*t)}


f1 <- function(t, x1, x2, x3) {
    ret <- 0
    if(x1 > 0 && x2 > 0) ret <- x2*exp(log(x1)*2/3)
    return(ret)
}

f2 <- function(t, x1, x2, x3) {
     ret <- 0
     if(x3 > 0) ret <- rho*sig*x3
     return(ret)
}

f3 <- function(t, x1, x2, x3) {
     ret <- 0
     if(x3 > 0) ret <- sqrt(1-rho^2)*sig*x3
     return(ret)
}

diff.coef.matrix <- matrix(c("f1(t,x1,x2,x3)",
 "f2(t,x1,x2,x3) * g(t)", "f2(t,x1,x2,x3)", "0", 
     "f3(t,x1,x2,x3)*g(t)", "f3(t,x1,x2,x3)"),  3, 2)

sabr.mod <- setModel(drift = c("0", "mu*g(t)*x3", "mu*x3"), 
diffusion = diff.coef.matrix, state.variable = c("x1", "x2", "x3"),
    solve.variable = c("x1", "x2", "x3"))
str(sabr.mod@parameter)

## ----echo=TRUE,results='hide'--------------------------------------------
 f2 <- function(t, x1, x2, x3, rho, sig) {
     ret <- 0
     if(x3 > 0) ret <- rho*sig*x3
     return(ret)
 }

 f3 <- function(t, x1, x2, x3, rho, sig) {
     ret <- 0
     if(x3 > 0) ret <- sqrt(1-rho^2)*sig*x3
     return(ret)
 }

 diff.coef.matrix <- matrix(c("f1(t,x1,x2,x3)", 
 "f2(t,x1,x2,x3,rho, sig) * g(t)", "f2(t,x1,x2,x3,rho,sig)", 
 "0", "f3(t,x1,x2,x3,rho,sig)*g(t)", "f3(t,x1,x2,x3,rho,sig)"),  3, 2)

 sabr.mod <- setModel(drift = c("0", "mu*g(t)*x3", "mu*x3"), 
 diffusion = diff.coef.matrix, state.variable = c("x1", "x2", "x3"), 
 solve.variable = c("x1", "x2", "x3"))
str(sabr.mod@parameter)

## ------------------------------------------------------------------------
Sigma <- matrix(c(0.5, 0.7, 0.7, 2), 2, 2)
C <- chol(Sigma)
C
crossprod(C)
Sigma

## ------------------------------------------------------------------------
set.seed(123)
drift <- c("mu*x1", "kappa*(theta-x2)")
diffusion <- matrix(c("c11*sqrt(x2)*x1", "0", 
  "c12*sqrt(x2)*x1", "c22*epsilon*sqrt(x2)"),2,2)
heston <- setModel(drift=drift, diffusion=diffusion,  
 state.var=c("x1","x2"))
X <- simulate(heston, true.par=list(theta=0.5, mu=1.2, kappa=2, 
 epsilon=0.2, c11=C[1,1], c12=C[1,2], c22=C[2,2]), 
 xinit=c(100,0.5))

## ----plot-heston,echo=FALSE, fig.keep='none',results='hide'--------------
pdf("figures/plot-heston.pdf",width=9,height=6)
set.seed(123)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ----echo=TRUE,results='hide'--------------------------------------------
ymodel <- setModel(drift="(2-theta2*x)", diffusion="(1+x^2)^theta1")
n <- 750
ysamp <- setSampling(Terminal = n^(1/3), n = n)
yuima <- setYuima(model = ymodel, sampling = ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit = 1, 
true.parameter = list(theta1 = 0.2, theta2 = 0.3))

## ----echo=TRUE,results='hide'--------------------------------------------
param.init <- list(theta2=0.5,theta1=0.5)
low.par <-  list(theta1=0, theta2=0)
upp.par <-  list(theta1=1, theta2=1)
mle1 <- qmle(yuima, start = param.init,  
  lower = low.par, upper = upp.par)

## ----results='hide'------------------------------------------------------
summary(mle1)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(mle1)),width=60))

## ----echo=TRUE-----------------------------------------------------------
prior <- list(theta2=list(measure.type="code",df="dunif(theta2,0,1)"),
 theta1=list(measure.type="code",df="dunif(theta1,0,1)"))

## ----echo=TRUE, results='hide'-------------------------------------------
lower <- list(theta1=0,theta2=0)
upper <- list(theta1=1,theta2=1)
bayes1 <- adaBayes(yuima, start=param.init, prior=prior,
 lower=lower,upper=upper, method="nomcmc")

## ----echo=TRUE-----------------------------------------------------------
coef(summary(bayes1))
coef(summary(mle1))

## ------------------------------------------------------------------------
n <- 500
ysamp <- setSampling(Terminal = n^(1/3), n = n)
yuima <- setYuima(model = ymodel, sampling = ysamp)
set.seed(123)
yuima <- simulate(yuima, xinit = 1, 
true.parameter = list(theta1 = 0.2, theta2 = 0.3))
param.init <- list(theta2=0.5,theta1=0.5)
lower <- list(theta1=0, theta2=0) 
upper <- list(theta1=1, theta2=1)
mle2 <- qmle(yuima, start =param.init , 
lower = lower, upper = upper)
bayes2 <- adaBayes(yuima, start=param.init, prior=prior,
 lower=lower,upper=upper)

## ------------------------------------------------------------------------
coef(summary(bayes2))
coef(summary(mle2))

## ------------------------------------------------------------------------
ymodel <- setModel(drift="(2-theta2*x)", diffusion="(1+x^2)^theta1")
n <- 100000
ysamp <- setSampling(delta=0.001, n = n)
mod <- setYuima(model = ymodel, sampling=ysamp)
set.seed(123)
yuima <- simulate(mod, xinit = 1, 
true.parameter = list(theta1 = 0.2, theta2 = 0.3))
param.init <- list(theta2=0.5,theta1=0.5)
yuima0.01 <- subsampling(yuima, 
 sampling=setSampling(delta=0.01,n=NULL,Terminal=100))
yuima0.1 <- subsampling(yuima, 
 sampling=setSampling(delta=0.1,n=NULL,Terminal=100))
yuima1.0 <- subsampling(yuima, 
 sampling=setSampling(delta=1,n=NULL,Terminal=100))

## ----echo=TRUE, fig.keep='none',results='hide'---------------------------
par(mfrow=c(2,2))
plot(yuima,main="delta=0.001, n=100000")
plot(yuima0.01,main="delta=0.01, n=10000")
plot(yuima0.1,main="delta=0.1, n=1000")
plot(yuima1.0,main="delta=1.0, n=100")

## ----plot-delta,echo=FALSE, fig.keep='none',results='hide'---------------
pdf("figures/plot-delta.pdf",width=9,height=6)
par(mar=c(4,4,3,1))
par(mfrow=c(2,2))
plot(yuima,main="delta=0.001, n=100000")
plot(yuima0.01,main="delta=0.01, n=10000")
plot(yuima0.1,main="delta=0.1, n=1000")
plot(yuima1.0,main="delta=1.0, n=100")
dev.off()

## ------------------------------------------------------------------------
low <- list(theta1=0, theta2=0)
up <-  list(theta1=1, theta2=1)
mle0.001 <- qmle(yuima, start = param.init, lower = low, upper = up)
summary(mle0.001)@coef

mle0.01 <- qmle(yuima0.01, start = param.init, lower = low, 
 upper = up)
summary(mle0.01)@coef

mle0.1 <- qmle(yuima0.1, start = param.init, lower = low, upper = up)
summary(mle0.1)@coef

mle1.0 <- qmle(yuima1.0, start = param.init, lower = low, upper = up)
summary(mle1.0)@coef

## ----echo=FALSE----------------------------------------------------------
est <- rbind( t(summary(mle0.001)@coef),  t(summary(mle0.01)@coef),  
 t(summary(mle0.1)@coef),   t(summary(mle1.0)@coef))

## ----message=FALSE-------------------------------------------------------
library(quantmod)
getSymbols("AAPL",to="2016-12-31")
head(AAPL)
S <- AAPL$AAPL.Adjusted

## ------------------------------------------------------------------------
Delta <- 1/252
gBm <- setModel(drift="mu*x", diffusion="sigma*x")
mod <- setYuima(model=gBm, data=setData(S, delta=Delta))

## ----appl,echo=TRUE,fig.keep='none',results='hide'-----------------------
set.seed(123)
plot(S)

## ----plot-aapl,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-aapl.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(S)
dev.off()

## ----results='hide'------------------------------------------------------
fit <- qmle(mod, start=list(mu=1, sigma=1), 
  lower=list(mu=0.1, sigma=0.1), 
  upper=list(mu=100, sigma=10))
summary(fit)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit)),width=60))

## ------------------------------------------------------------------------
X <- diff(log(S))
X <- as.numeric(na.omit(diff(log(S))))
alpha <- mean(X)/Delta
sigma <- sqrt(var(X)/Delta)
mu <- alpha +0.5*sigma^2
mu
sigma
coef(fit)

## ------------------------------------------------------------------------
getSymbols("DEXUSEU", src="FRED")
DEXUSEU <- DEXUSEU["/2016"]
head(DEXUSEU)
meanCIR <- mean(DEXUSEU, na.rm=TRUE)
meanCIR

## ----dexuseu,echo=TRUE,fig.keep='none',results='hide'--------------------
set.seed(123)
plot(DEXUSEU)

## ----plot-dexuseu,echo=FALSE, fig.keep='none',results='hide'-------------
pdf("figures/plot-dexuseu.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(DEXUSEU)
dev.off()

## ------------------------------------------------------------------------
cir1 <- setModel(drift="theta1-theta2*x", diffusion="sigma*sqrt(x)")
cir2 <- setModel(drift="kappa*(mu-x)", diffusion="sigma*sqrt(x)")
mod1 <- setYuima(model=cir1, data=setData(na.omit(DEXUSEU), 
 delta=Delta))
mod2 <- setYuima(model=cir2, data=setData(na.omit(DEXUSEU), 
 delta=Delta))

## ----results='hide'------------------------------------------------------
fit1 <- qmle(mod1, start=list(theta1=1, theta2=1, sigma=0.5),  
  lower=list(theta1=0.1, theta2=0.1, sigma=0.1),
  upper=list(theta1=10, theta2=10, sigma=100),
   method="L-BFGS-B")
summary(fit1)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit1)),width=60))

## ----results='hide'------------------------------------------------------
fit2 <- qmle(mod2, start=list(kappa=1, mu=meanCIR, sigma=0.5),  
  lower=list(kappa=0.1, mu=0.1, sigma=0.1),
  upper=list(kappa=10, mu=10, sigma=100),
   method="L-BFGS-B")
summary(fit2)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit2)),width=60))

## ------------------------------------------------------------------------
theta1 <- as.numeric( coef(fit2)["kappa"] * coef(fit2)["mu"] )
theta1
coef(fit1)["theta1"]

## ------------------------------------------------------------------------
mu <- as.numeric( coef(fit1)["theta1"] / coef(fit1)["theta2"] )
mu
coef(fit2)["mu"]

## ------------------------------------------------------------------------
model<- setModel(drift="t1*(t2-x)",diffusion="t3")

## ------------------------------------------------------------------------
T<-300
n<-3000
sampling <- setSampling(Terminal=T, n=n)
yuima<-setYuima(model=model, sampling=sampling)
h00 <- list(t1=0.3, t2=1, t3=0.25)
h01 <- list(t1=0.3, t2=0.2, t3=0.1)
set.seed(123)
X <- simulate(yuima, xinit=1, true.par=h00)

## ------------------------------------------------------------------------
phi1 <- function(x) 1-x+x*log(x)

## ------------------------------------------------------------------------
phi.test(X, H0=h00, phi=phi1, start=h00,
   lower=list(t1=0.1, t2=0.1, t3=0.1), 
   upper=list(t1=2,t2=2,t3=2),method="L-BFGS-B")

## ----echo=FALSE,results='hide'-------------------------------------------
pval <- phi.test(X, H0=h00, phi=phi1, start=h00,
   lower=list(t1=0.1, t2=0.1, t3=0.1), 
   upper=list(t1=2,t2=2,t3=2),method="L-BFGS-B")$pvalue

## ------------------------------------------------------------------------
phi.test(X, H0=h01, phi=phi1, start=h00, 
  lower=list(t1=0.1, t2=0.1, t3=0.1), 
  upper=list(t1=2,t2=2,t3=2),method="L-BFGS-B")

## ------------------------------------------------------------------------
library(quantmod)
Delta <- 1/252
getSymbols("DEXUSEU", src="FRED")
DEXUSEU <- DEXUSEU["/2016"]
USEU <- setData(na.omit(DEXUSEU),  delta=Delta)
meanCIR <- mean(get.zoo.data(USEU)[[1]])
gBm <- setModel(drift="mu*x", diffusion="sigma*x")
mod <- setYuima(model=gBm, data=USEU)
cir1 <- setModel(drift="theta1-theta2*x", diffusion="sigma*sqrt(x)")
cir2 <- setModel(drift="kappa*(mu-x)", diffusion="sigma*sqrt(x)")
ckls <- setModel(drift="theta1-theta2*x", diffusion="sigma*x^gamma")
mod1 <- setYuima(model=cir1, data=USEU)
mod2 <- setYuima(model=cir2, data=USEU)
mod3 <- setYuima(model=ckls, data=USEU)
gBm.fit <- qmle(mod, start=list(mu=1, sigma=1), 
  lower=list(mu=0.1, sigma=0.1), 
  upper=list(mu=100, sigma=10))
cir1.fit <- qmle(mod1, start=list(theta1=1, theta2=1, sigma=0.5),  
  lower=list(theta1=0.1, theta2=0.1, sigma=0.1),
  upper=list(theta1=10, theta2=10, sigma=100),
   method="L-BFGS-B")
cir2.fit <- qmle(mod2, start=list(kappa=1, mu=meanCIR, sigma=0.5),  
  lower=list(kappa=0.1, mu=0.1, sigma=0.1),
  upper=list(kappa=10, mu=10, sigma=100),
   method="L-BFGS-B")
ckls.fit <- qmle(mod3, start=list(theta1=1, theta2=1, sigma=0.5, 
  gamma=0.5),  lower=list(theta1=0.1, theta2=0.1, sigma=0.1, 
  gamma=0.1), upper=list(theta1=10, theta2=10, sigma=10, 
  gamma=2), method="L-BFGS-B")

## ------------------------------------------------------------------------
AIC(gBm.fit,cir1.fit,cir2.fit,ckls.fit)

## ----echo=FALSE----------------------------------------------------------
tmp <- AIC(gBm.fit,cir1.fit,cir2.fit,ckls.fit)

## ------------------------------------------------------------------------
set.seed(123)
S <- simulate(gBm, true.par=list(mu=1, sigma=0.25), 
  sampling=setSampling(T=1, n=1000), xinit=100)
mod <- setYuima(model=gBm, data=S@data)
mod1 <- setYuima(model=cir1, data=S@data)
mod2 <- setYuima(model=cir2, data=S@data)
mod3 <- setYuima(model=ckls, data=S@data)
 gBm.fit <- qmle(mod, start=list(mu=1, sigma=1), 
  lower=list(mu=0.1, sigma=0.1), 
  upper=list(mu=100, sigma=10))
cir1.fit <- qmle(mod1, start=list(theta1=1, theta2=1, sigma=0.5),  
  lower=list(theta1=0.1, theta2=0.1, sigma=0.1),
  upper=list(theta1=10, theta2=10, sigma=100),
   method="L-BFGS-B")
cir2.fit <- qmle(mod2, start=list(kappa=1, mu=meanCIR, sigma=0.5),  
  lower=list(kappa=0.1, mu=0.1, sigma=0.1),
  upper=list(kappa=10, mu=10, sigma=100),
   method="L-BFGS-B")
ckls.fit <- qmle(mod3, 
  start=list(theta1=1, theta2=1, sigma=0.5, gamma=0.5),  
  lower=list(theta1=0.1, theta2=0.1, sigma=0.1, gamma=0.1),
  upper=list(theta1=10, theta2=10, sigma=10, gamma=2),
   method="L-BFGS-B")

## ------------------------------------------------------------------------
AIC(gBm.fit,cir1.fit,cir2.fit,ckls.fit)

## ----echo=FALSE----------------------------------------------------------
tmp <- AIC(gBm.fit,cir1.fit,cir2.fit,ckls.fit)

## ------------------------------------------------------------------------
a <- c("1-mu11*X1+mu12*X2","2+mu21*X1-mu22*X2")
b <- matrix(c("s1*X1","s2*X1", "-s3*X2","s4*X2"),2,2)
mod.est <- setModel(drift=a, diffusion=b,
 solve.var=c("X1","X2"),state.variable=c("X1","X2"))
truep <- list(mu11=.9, mu12=0, mu21=0, mu22=0.7, 
 s1=.3, s2=0,s3=0,s4=.2)
low <- list(mu11=1e-5, mu12=1e-5, mu21=1e-5, mu22=1e-5, 
 s1=1e-5, s2=1e-5, s3=1e-5,s4=1e-5)
upp <- list(mu11=2, mu12=2, mu21=1, mu22=1, 
 s1=1, s2=1, s3=1,s4=1)
set.seed(123)
n <- 1000
X <- simulate(mod.est, T=n^(1/3), n=n, xinit=c(1,1), 
 true.parameter=truep)

## ----results='hide'------------------------------------------------------
myest <- lasso(X,  delta=2, start=truep, lower=low, upper=upp, 
 method="L-BFGS-B")
myest

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(myest),width=60))

## ------------------------------------------------------------------------
fit1 <- qmle(X, start=truep, lower=low, upper=upp, method="L-BFGS-B")

## ------------------------------------------------------------------------
a <- c("1-mu11*X1","2-mu22*X2")
b <- matrix(c("s1*X1","0", "0","s4*X2"),2,2)
mod.est2 <- setModel(drift=a, diffusion=b,
 solve.var=c("X1","X2"),state.variable=c("X1","X2"))
truep <- list(mu11=.9, mu22=0.7, s1=.3,s4=.2)
low <- list(mu11=1e-5,   mu22=1e-5, s1=1e-5,  s4=1e-5)
upp <- list(mu11=2,   mu22=2, s1=1,  s4=1)
Y <- setYuima(model=mod.est2, data=X@data)
fit2 <- qmle(Y, start=truep, lower=low, upper=upp, method="L-BFGS-B")
summary(fit1)
summary(fit2)
AIC(fit1, fit2)

## ----results='hide',fig.keep='none',message=FALSE------------------------
library("Ecdat")
data("Irates")
rates <- Irates[,"r1"]
plot(rates)
X <- window(rates, start=1964.471, end=1989.333)
mod <- setModel(drift="alpha+beta*x", diffusion="sigma*x^gamma")
yuima <- setYuima(data=setData(X,delta=1/12), model=mod)
start <- list(alpha=1, beta =-.1, sigma =.1, gamma =1)
low <- list(alpha=-5, beta =-5, sigma =-5, gamma =-5)
upp <- list(alpha=8, beta =8, sigma =8, gamma =8)

## ----plot-irates,echo=FALSE, fig.keep='none',results='hide'--------------
pdf("figures/plot-irates.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(rates)
dev.off()

## ----echo=TRUE,results='hide'--------------------------------------------
lasso.est <- lasso(yuima, start=start, lower=low, upper=upp,
   method="L-BFGS-B", delta=2)
lasso.est

## ----results='hide'------------------------------------------------------
mod1 <- setModel(drift="alpha", diffusion="sigma*x^gamma")
yuima1 <- setYuima(data=setData(X,delta=1/12), model=mod1)
start1 <- list(alpha=1, sigma =.1, gamma =1)
low1 <- list(alpha=-5, sigma =-5, gamma =-5)
upp1 <- list(alpha=8, sigma =8, gamma =8)
fit <- qmle(yuima, start=start, lower=low, upper=upp,
   method="L-BFGS-B")
fit1 <- qmle(yuima1, start=start1, lower=low1, upper=upp1,
   method="L-BFGS-B")
summary(fit)
summary(fit1)
AIC(fit, fit1)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit)),width=60))
writeLines(strwrap(capture.output(summary(fit1)),width=60))
writeLines(strwrap(capture.output(AIC(fit, fit1)),width=60))

## ----cpoint1-------------------------------------------------------------
diff.matrix <- matrix(c("theta1.k*x1","0*x2","0*x1","theta2.k*x2"), 
 2, 2)
drift.c <- c("sin(x1)", "3-x2")
drift.matrix <- matrix(drift.c, 2, 1)
ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, 
 time.variable="t", state.variable=c("x1", "x2"), 
 solve.variable=c("x1", "x2"))
ymodel

## ----cpoint3,results='hide'----------------------------------------------
n <- 1000

set.seed(123)

t0 <- list(theta1.k=0.5, theta2.k=0.3)
T <- 10
tau <- 4
pobs <- tau/T
ysamp1 <- setSampling(n=n*pobs, Initial=0, delta=0.01)
yuima1 <- setYuima(model=ymodel, sampling=ysamp1)
yuima1 <- simulate(yuima1, xinit=c(3, 3), true.parameter=t0)

v11 <- get.zoo.data(yuima1)[[1]] 
x1 <- as.numeric(v11[length(v11)]) # terminal value
v21 <- get.zoo.data(yuima1)[[2]]
x2 <- as.numeric(v21[length(v21)]) # terminal value

## ----cpoint3b,results='hide'---------------------------------------------
t1 <- list(theta1.k=0.2, theta2.k=0.4)
ysamp2 <- setSampling(Initial=n*pobs*0.01, n=n*(1-pobs), delta=0.01)
yuima2 <- setYuima(model=ymodel, sampling=ysamp2)
yuima2 <- simulate(yuima2, xinit=c(x1, x2), true.parameter=t1)

## ----cpoint3c,results='hide'---------------------------------------------
v12 <- get.zoo.data(yuima2)[[1]] 
v22 <- get.zoo.data(yuima2)[[2]]
v1 <- c(v11,v12[-1])
v2 <- c(v21,v22[-1])
new.data <- setData(zoo(cbind(v1,v2)),delta=0.01)
yuima <- setYuima(model=ymodel, data=new.data)

## ----cpoint4,fig.keep='none'---------------------------------------------
plot(yuima)

## ----plot-cpoint4,echo=FALSE, fig.keep='none',results='hide'-------------
pdf("figures/plot-cpoint4.pdf",width=9,height=5)
par(mar=c(4,4,1,1))
plot(yuima)
dev.off()

## ----cpoint4b------------------------------------------------------------
noDriftModel <- setModel(drift=c(0,0), diffusion=diff.matrix,
 time.variable="t", state.variable=c("x1", "x2"), 
 solve.variable=c("x1", "x2"))
noDriftModel <- setYuima(noDriftModel, data=new.data)
noDriftModel@model@drift
noDriftModel

## ----cpoint5-------------------------------------------------------------
t.est <- CPoint(yuima,param1=t0,param2=t1)
t.est$tau
t.est2 <- CPoint(noDriftModel,param1=t0,param2=t1)
t.est2$tau

## ------------------------------------------------------------------------
qmleL(noDriftModel, t=1.5, start=list(theta1.k=0.1, theta2.k=0.1),
 lower=list(theta1.k=0, theta2.k=0), 
 upper=list(theta1.k=1, theta2.k=1), 
 method="L-BFGS-B") -> estL
qmleR(noDriftModel, t=8.5, start=list(theta1.k=0.1, theta2.k=0.1),
 lower=list(theta1.k=0, theta2.k=0), 
 upper=list(theta1.k=1, theta2.k=1), 
 method="L-BFGS-B") -> estR
t0.est <- coef(estL)
t1.est <- coef(estR)

## ------------------------------------------------------------------------
t.est3 <- CPoint(noDriftModel,param1=t0.est,param2=t1.est)
t.est3

## ----eval=FALSE----------------------------------------------------------
## CPoint(noDriftModel,param1=t0.est,param2=t1.est, plot=TRUE)

## ----plot-cpoint-stat,echo=FALSE,results='hide'--------------------------
pdf("figures/plot-cpoint-stat.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
CPoint(noDriftModel,param1=t0.est,param2=t1.est, plot=TRUE)
dev.off()

## ------------------------------------------------------------------------
qmleL(noDriftModel, t=t.est3$tau, 
 start=list(theta1.k=0.1, theta2.k=0.1),
 lower=list(theta1.k=0, theta2.k=0), 
 upper=list(theta1.k=1, theta2.k=1), 
 method="L-BFGS-B") -> estL
qmleR(noDriftModel, t=t.est3$tau, 
 start=list(theta1.k=0.1, theta2.k=0.1),
 lower=list(theta1.k=0, theta2.k=0), 
 upper=list(theta1.k=1, theta2.k=1), 
 method="L-BFGS-B") -> estR
t02s.est <- coef(estL)
t12s.est <- coef(estR)
t2s.est3 <- CPoint(noDriftModel,param1=t02s.est,param2=t12s.est)
t2s.est3

## ------------------------------------------------------------------------
library(quantmod)
getSymbols("AAPL", to="2016-12-31")
S <- AAPL$AAPL.Adjusted
Delta <- 1/252
gBm <- setModel(drift="mu*x", diffusion="sigma*x")
mod <- setYuima(model=gBm, data=setData(S, delta=Delta))
lower <- list(mu=0.1, sigma=0.1)
upper <- list(mu=100, sigma=10)
start <- list(mu=1, sigma=1)
fit <- qmle(mod, start= start, upper=upper, lower=lower)
summary(fit)

## ------------------------------------------------------------------------
fit1 <- qmleL(mod, t=1, start= list(mu=1,sigma=1))
fit2 <- qmleR(mod, t=6, start= list(mu=1,sigma=1))
fit1
fit2

## ------------------------------------------------------------------------
cp <- CPoint(mod,param1=coef(fit1),param2=coef(fit2))
cp

## ----fig.keep='none'-----------------------------------------------------
X <- diff(log(get.zoo.data(mod)[[1]]))
plot(X)
abline(v=cp$tau, lty=3,lwd=2,col="red")

## ----plot-returns,echo=FALSE, fig.keep='none',results='hide'-------------
pdf("figures/plot-returns.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(X,main="log returns of AAPL")
abline(v=cp$tau, lty=3,lwd=2,col="red")
dev.off()

## ----plot-cpoint-aapl,echo=FALSE,results='hide'--------------------------
pdf("figures/plot-cpoint-aapl.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
CPoint(mod,param1=coef(fit1),param2=coef(fit2),plot=TRUE)
dev.off()

## ----echo=TRUE-----------------------------------------------------------
# diffusion coefficient for process 1
diff.coef.1 <- function(t,x1=0, x2=0) sqrt(1+t)
# diffusion coefficient for process 2
diff.coef.2 <- function(t,x1=0, x2=0) sqrt(1+t^2)
# correlation
cor.rho <- function(t,x1=0, x2=0) sqrt(1/2)
# coefficient matrix for diffusion term
diff.coef.matrix <- matrix( c( "diff.coef.1(t,x1,x2)", 
"diff.coef.2(t,x1,x2) * cor.rho(t,x1,x2)", "", 
"diff.coef.2(t,x1,x2) * sqrt(1-cor.rho(t,x1,x2)^2)"),2,2)
# Model SDE using yuima.model
cor.mod <- setModel(drift = c("",""), diffusion = diff.coef.matrix,
 solve.variable=c("x1","x2"))

## ----echo=TRUE-----------------------------------------------------------
CC.theta <- function( T, sigma1, sigma2, rho){
 tmp <- function(t) return( sigma1(t) * sigma2(t) * rho(t) )
 integrate(tmp,0,T)
}

## ------------------------------------------------------------------------
set.seed(123)
Terminal <- 1
n <- 1000
# Cumulative Covariance
theta <- CC.theta(T=Terminal, sigma1=diff.coef.1, 
sigma2=diff.coef.2, rho=cor.rho)$value
cat(sprintf("theta=%5.3f\n",theta))

## ----results='hide'------------------------------------------------------
yuima.samp <- setSampling(Terminal=Terminal,n=n)
yuima <- setYuima(model=cor.mod, sampling=yuima.samp)
X <- simulate(yuima)

## ------------------------------------------------------------------------
cce(X)

## ----cceplot1,fig.keep='none'--------------------------------------------
plot(X,main="complete data")

## ----plot-cceplot1,echo=FALSE, fig.keep='none',results='hide'------------
pdf("figures/plot-cceplot1.pdf",width=9,height=5)
par(mar=c(4,4,1,1))
plot(X,main="complete data")
dev.off()

## ------------------------------------------------------------------------
p1 <- 0.2
p2 <- 0.3
newsamp <- setSampling(random=list(rdist=c( 
  function(x) rexp(x, rate=p1*n/Terminal), 
  function(x) rexp(x, rate=p2*n/Terminal))) )

## ------------------------------------------------------------------------
Y <- subsampling(X, sampling=newsamp)

## ----cceplot2,fig.keep='none'--------------------------------------------
plot(Y,main="asynchronous data")

## ----plot-cceplot2,echo=FALSE, fig.keep='none',results='hide'------------
pdf("figures/plot-cceplot2.pdf",width=9,height=5)
par(mar=c(4,4,1,1))
plot(Y,main="asynchronous data")
dev.off()

## ------------------------------------------------------------------------
cce(Y)$covmat   # asynch data
cce(X)$covmat   # full data

## ------------------------------------------------------------------------
b1 <- function(x,y) y
b2 <- function(x,y) -x
s1 <- function(t,x,y) sqrt(abs(x)*(1+t))
s2 <- function(t,x,y) sqrt(abs(y))
cor.rho <- function(t,x,y) 1/(1+x^2)
diff.mat <- matrix(c("s1(t,x,y)", "s2(t,x,y) * cor.rho(t,x,y)","",
 "s2(t,x,y) * sqrt(1-cor.rho(t,x,y)^2)"), 2, 2) 
cor.mod <- setModel(drift = c("b1","b2"), diffusion = diff.mat,
solve.variable = c("x", "y"),state.var=c("x","y"))

## Generate a path of the process
set.seed(111) 
Terminal <- 1
n <- 10000
yuima.samp <- setSampling(Terminal = Terminal, n = n) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima, xinit=c(2,3)) 

## ------------------------------------------------------------------------
p1 <- 0.2
p2 <- 0.3
newsamp <- setSampling(random=list(rdist=c( 
 function(x) rexp(x, rate=p1*n/Terminal), 
 function(x) rexp(x, rate=p2*n/Terminal))) )
Y <- subsampling(yuima, sampling = newsamp)

## ----cceplot3,fig.keep='none'--------------------------------------------
plot(Y,main="asynchronous data (non linear case)")

## ----plot-cceplot3,echo=FALSE, fig.keep='none',results='hide'------------
pdf("figures/plot-cceplot3.pdf",width=9,height=5)
par(mar=c(4,4,1,1))
plot(Y,main="asynchronous data (non linear case)")
dev.off()

## ------------------------------------------------------------------------
cce(yuima)$covmat # full data
cce(Y)$covmat        # asynch data

## ------------------------------------------------------------------------
diff.coef.matrix <- matrix(c("sqrt(x1)", "3/5*sqrt(x2)",
 "1/3*sqrt(x3)", "", "4/5*sqrt(x2)","2/3*sqrt(x3)","","",
  "2/3*sqrt(x3)"), 3, 3) 
drift <- c("1-x1","2*(10-x2)","3*(4-x3)")
cor.mod <- setModel(drift = drift, diffusion = diff.coef.matrix,
  solve.variable = c("x1", "x2","x3")) 

set.seed(111) 
Terminal <- 1
yuima.samp <- setSampling(Terminal = Terminal, n = 1200) 
yuima <- setYuima(model = cor.mod, sampling = yuima.samp) 
yuima <- simulate(yuima, xinit=c(1,7,5)) 

# intentionally displace the second time series

data1 <- get.zoo.data(yuima)[[1]]
data2 <- get.zoo.data(yuima)[[2]]
time2 <- time( data2 )
theta2 <- 0.05   # the lag of x2 behind x1
stime2 <- time2 + theta2  
time(data2) <- stime2
data3 <- get.zoo.data(yuima)[[3]]
time3 <- time( data3 )
theta3 <- 0.12   # the lag of x3 behind x1
stime3 <- time3 + theta3 
time(data3) <- stime3
syuima <- setYuima(data=setData(merge(data1, data2, data3)))
yuima
syuima

## ----shifted,fig.keep='none'---------------------------------------------
plot(syuima,main="time shifted data")

## ----plot-shifted,echo=FALSE, fig.keep='none',results='hide'-------------
pdf("figures/plot-shifted.pdf",width=9,height=5)
par(mar=c(4,4,2,1))
plot(syuima,main="time shifted data")
dev.off()

## ------------------------------------------------------------------------
llag(yuima)
llag(syuima)

## ----plot-shifted-ci,echo=FALSE, fig.keep='none',results='hide'----------
pdf("figures/plot-shifted-ci.pdf",width=9,height=5)
par(mar=c(4,5,2,1))
par(mfrow=c(1,3))
llag(syuima,plot=TRUE,ci=TRUE)
dev.off()

## ------------------------------------------------------------------------
data2 <- get.zoo.data(yuima)[[2]]
time2 <- time( data2 )
theta2 <- 0.05   # the lag of x2 behind x1
stime2 <- time2 + theta2  
time(data2) <- stime2
data3 <- get.zoo.data(yuima)[[3]]
time3 <- time( data3 )
theta3 <- 0.12   # the lag of x3 behind x1
stime3 <- time3 + theta3 
time(data3) <- stime3
data1 <- data1[which(time(data1)>0.5 & time(data1)<1)]
data2 <- data2[which(time(data2)>0.5 & time(data2)<1)]
data3 <- data3[which(time(data3)>0.5 & time(data3)<1)]
syuima2 <- setYuima(data=setData(merge(data1, data2, data3)))
syuima2
llag(syuima2)

## ------------------------------------------------------------------------
p1 <- 0.2
p2 <- 0.3
p3 <- 0.4
n <- 1000
newsamp <- setSampling(
random=list(rdist=c( function(x) rexp(x, rate=p1*n/Terminal), 
function(x) rexp(x, rate=p2*n/Terminal),
function(x) rexp(x, rate=p3*n/Terminal))) )
psample <- subsampling(syuima, sampling = newsamp)
psample
llag(psample)

## ----results='hide'------------------------------------------------------
library(quantmod)
getSymbols("AAPL", from="2013-01-01", to="2013-12-31")
getSymbols("IBM", from="2013-01-01", to="2013-12-31")
getSymbols("AMZN", from="2013-01-01", to="2013-12-31")
getSymbols("EBAY", from="2013-01-01", to="2013-12-31")
getSymbols("FB", from="2013-01-01", to="2013-12-31")
getSymbols("MSFT", from="2013-01-01", to="2013-12-31")
data1 <- AAPL$AAPL.Close
data2 <- IBM$IBM.Close
data3 <- AMZN$AMZN.Close
data4 <- EBAY$EBAY.Close
data5 <- FB$FB.Close
data6 <- MSFT$MSFT.Close
market.data <- merge(data1, data2, data3, data4,data5,data6)
colnames(market.data) <- c("AAPL", "IBM", "AMZN", "EBAY", 
 "FB", "MSFT")
mkt <- setYuima(data=setData(market.data, delta=1/252))

## ------------------------------------------------------------------------
mkt
round(cce(mkt)$cormat,2) # correlation matrix

## ----market,fig.keep='none'----------------------------------------------
plot(mkt)

## ----plot-market,echo=FALSE, fig.keep='none',results='hide'--------------
pdf("figures/plot-market.pdf",width=9,height=5)
par(mar=c(4,4,1,1))
plot(mkt,main="")
dev.off()

## ------------------------------------------------------------------------
round(llag(mkt),4)

## ----corrplot,fig.keep='none',message=FALSE------------------------------
require(corrplot)
cols <- colorRampPalette(c("#7F0000", "red", "#FF7F00",
  "yellow", "white", "cyan", 
  "#007FFF", "blue", "#00007F"))
corrplot(cce(mkt)$cormat,method="ellipse", 
 cl.pos = "b", tl.pos = "d", tl.srt = 60, 
 col=cols(100), outline=TRUE)
corrplot(llag(mkt),method="ellipse",is.corr=FALSE,
 cl.pos = "b", tl.pos = "d", tl.srt = 60, 
 col=cols(100), outline=TRUE)

## ----plot-corrplot,echo=FALSE, fig.keep='none',results='hide'------------
pdf("figures/plot-corrplot1.pdf",width=6,height=6)
require(corrplot)
corrplot(cce(mkt)$cormat,method="ellipse", 
 cl.pos = "b", tl.pos = "d", tl.srt = 60, col=cols(100), outline=TRUE)
dev.off()
pdf("figures/plot-corrplot2.pdf",width=6,height=6)
corrplot(llag(mkt),method="ellipse",is.corr=FALSE,
 cl.pos = "b", tl.pos = "d", tl.srt = 60, col=cols(100), outline=TRUE)
dev.off()

## ----echo=TRUE, results='hide'-------------------------------------------
model <- setModel(drift = "x", diffusion = matrix( "x*e", 1,1))
T <- 1
xinit <- 150
K <- 100
f <- list( expression(x/T), expression(0))
F <- 0
e <- 0.5
yuima <- setYuima(model = model, 
  sampling = setSampling(Terminal=T, n=1000))
yuima <- setFunctional( yuima, f=f,F=F, xinit=xinit,e=e)

## ----echo=TRUE-----------------------------------------------------------
str(yuima@functional)

## ----echo=TRUE-----------------------------------------------------------
F0 <- F0(yuima)
F0

## ----echo=TRUE,results='hide'--------------------------------------------
rho <- expression(0)
epsilon <- e  # noise level
g <- function(x) {
  tmp <- (F0 - K) + (epsilon * x) 
 tmp[(epsilon * x) < (K-F0)] <- 0
 tmp
}

## ----echo=TRUE,results='hide'--------------------------------------------
asymp <- asymptotic_term(yuima, block=10, rho, g)
asymp

## ----echo=TRUE-----------------------------------------------------------
asy1 <- asymp$d0 + e * asymp$d1 
# 1st order asymp. exp. of asian call price
asy1
asy2 <- asymp$d0 + e * asymp$d1 +  e^2* asymp$d2 
# 2nd order asymp. exp. of asian call price
asy2

## ----message=FALSE-------------------------------------------------------
library("fExoticOptions")
levy <- LevyAsianApproxOption(TypeFlag = "c", S = xinit, SA = xinit, 
    X = K, Time = 1, time = 1, r = 0.0, b = 1, sigma = e)@price
levy

## ------------------------------------------------------------------------
a <- 0.9
e <- 0.4
Terminal <- 3

xinit <- 1
K <- 10

drift <- "a * x"
diffusion <- "e * sqrt(x)"

model <- setModel(drift=drift,diffusion=diffusion)

n <- 1000*Terminal
yuima <- setYuima(model = model,
  sampling = setSampling(Terminal=Terminal,n=n))

f <- list(c(expression(0)),c(expression(0)))
F <- expression(x)

yuima.ae <- setFunctional(yuima,f=f,F=F,xinit=xinit,e=e)
rho <- expression(0)
F1 <- F0(yuima.ae)


get_ge <- function(x,epsilon,K,F0){
	tmp <- (F0 - K) + (epsilon * x[1])
	tmp[(epsilon * x[1]) > (K - F0)] <- 0
	return( - tmp )
}


g <- function(x){
	return(get_ge(x,e,K,F1))
}

time1 <- proc.time()
asymp <- asymptotic_term(yuima.ae,block=100,rho,g)
time2 <- proc.time()

## ------------------------------------------------------------------------
ae.value0 <- asymp$d0
ae.value0
ae.value1 <- asymp$d0 + e * asymp$d1
ae.value1
ae.value2 <- as.numeric(asymp$d0 + e * asymp$d1 + e^2 * asymp$d2)
ae.value2
ae.time <- time2 - time1
ae.time

