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

## ----mod1----------------------------------------------------------------
mod1 <- setPoisson(intensity="lambda", df=list("dconst(z,1)"))
mod1

## ----poi1,fig.keep='none'------------------------------------------------
Terminal <- 30
samp <- setSampling(T=Terminal,n=3000)
set.seed(123)
poisson1 <- simulate(mod1, true.par=list(lambda=1),sampling=samp)
poisson1
plot(poisson1)

## ----plot-poi1,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi1.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson1,type="S")
dev.off()

## ----eval=FALSE----------------------------------------------------------
## setPoisson(intensity="lambda", df=list("dconst(z,1)"), scale=5)
## setPoisson(intensity="lambda", df=list("dconst(z,5)"))

## ----mod2,fig.keep='none'------------------------------------------------
mod2 <- setPoisson(intensity="lambda", df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson2 <- simulate(mod2, sampling=samp,
 true.par=list(lambda=1,mu=0, sigma=2))
poisson2
plot(poisson2)

## ----plot-poi2,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi2.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson2,type="S")
dev.off()

## ----fig.keep='none'-----------------------------------------------------
mod3 <- setPoisson(intensity="lambda", 
 df=list("dNIG(z,alpha,beta,gamma,mu)"))
poisson3 <- simulate(mod3, sampling=samp,
 true.par=list(lambda=10,alpha=2,beta=0.3,gamma=1,mu=0))
poisson3

## ----message=FALSE-------------------------------------------------------
require(fBasics)
mod4 <- setPoisson(intensity="lambda", 
 df=list("dnig(z,alpha,beta,gamma)"))
poisson4 <- simulate(mod4,  sampling=samp,
 true.par=list(lambda=10,alpha=2,beta=0.3,gamma=1))
poisson4

## ----mod5,fig.keep='none'------------------------------------------------
mod5 <- setPoisson(intensity="alpha+beta*t", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson5 <- simulate(mod5,  sampling=samp,
 true.par=list(alpha=2,beta=.5,mu=0, sigma=2))
plot(poisson5)
f <- function(t,alpha,beta) alpha + beta*t
curve(f(x,alpha=2,beta=0.5)-20,0,30,add=TRUE,col="red",lty=3,lwd=2)

## ----plot-poi5,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi5.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson5,type="S")
curve(f(x,alpha=2,beta=0.5)-20,0,30,add=TRUE,col="red",lty=3,lwd=2)
dev.off()

## ----mod6,fig.keep='none'------------------------------------------------
mod6 <- setPoisson(intensity="theta*t^(theta-1)", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson6 <- simulate(mod6,  sampling=samp,
 true.par=list(theta=1.5,mu=0, sigma=2))
plot(poisson6)
f <- function(t,theta) theta*t^(theta-1)
curve(f(x,theta=1.5),0,30,add=TRUE,col="red",lty=3,lwd=2)

## ----plot-poi6,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi6.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson6,type="S")
curve(f(x,theta=1.5),0,30,add=TRUE,col="red",lty=3,lwd=2)
dev.off()

## ----mod7,fig.keep='none'------------------------------------------------
mod7 <- setPoisson(intensity="beta*exp(-lambda*t)", 
 df=list("dexp(z,gamma)"))
set.seed(123)
poisson7 <- simulate(mod7, sampling=samp,
 true.par=list(lambda=.2,beta=10,gamma=1))
plot(poisson7)
f <- function(t,beta,lambda) beta*exp(-lambda*t)
curve(f(x,beta=10,lambda=0.2),0,30,add=TRUE,col="red",lty=3,lwd=2)

## ----plot-poi7,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi7.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson7,type="S")
curve(f(x,beta=10,lambda=0.2),0,30,add=TRUE,col="red",lty=3,lwd=2)
dev.off()

## ----mod8,fig.keep='none'------------------------------------------------
mod8 <- setPoisson(intensity="0.5*a*(1+cos(omega*t+phi))+lambda", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson8 <- simulate(mod8, sampling=samp,
 true.par=list(a=2,omega=0.5,phi=3.14,lambda=5,mu=0,sigma=1))
plot(poisson8)
f <- function(t,a,omega,phi,lambda) 0.5*a*(1+cos(omega*t+phi))+lambda
curve(f(x,a=2,omega=0.5,phi=3.14,lambda=5),0,30,add=TRUE,
 col="red",lty=3,lwd=2)

## ----plot-poi8,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi8.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson8,type="S")
curve(f(x,a=2,omega=0.5,phi=3.14,lambda=5),0,30,add=TRUE, col="red",lty=3,lwd=2)
dev.off()

## ----fig.keep='none'-----------------------------------------------------
mod9 <- setPoisson(intensity="a*cos(theta*t)+lambda", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson9 <- simulate(mod9, sampling=samp, 
true.par=list(a=1,theta=0.5,lambda=5,mu=0,sigma=1))
plot(poisson9)
f <- function(t,a,theta,lambda) a*cos(theta*t)+lambda
curve(f(x,a=1,theta=0.5,lambda=5),0,30,add=TRUE,col="red",lty=3,lwd=2)

## ----plot-poi9,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-poi9.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson9,type="S")
curve(f(x,a=1,theta=0.5,lambda=5),0,30,add=TRUE,col="red",lty=3,lwd=2)
dev.off()

## ----mod10,fig.keep='none'-----------------------------------------------
mod10 <- setPoisson(intensity="lambda*t", 
 df=list("dmvnorm(z,c(0.15,-0.1),matrix(c(2,-1.9,-1.9,4.3),2,2))"),
 dimension=2)
set.seed(123)
poisson10 <- simulate(mod10, true.par=list(lambda=5), sampling=samp)
poisson10
plot(poisson10)

## ----plot-poi10,echo=FALSE,results='hide'--------------------------------
pdf("figures/plot-poi10.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson10,type="S")
dev.off()

## ----mod111,fig.keep='none'----------------------------------------------
mod11 <- setPoisson(intensity="lambda*t", 
 df=list("dmvnorm(z,c(0.01,-0.01,.05),
 matrix(c(1,.5,0,.5,1,0,0,0,1),3,3))"),
 dimension=3)
set.seed(123)
poisson11 <- simulate(mod11, true.par=list(lambda=5), 
sampling=samp,xinit=c(-100,200,300))
plot(poisson11)

## ----plot-poi11,echo=FALSE,results='hide'--------------------------------
pdf("figures/plot-poi11.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(poisson11,type="S")
dev.off()

## ------------------------------------------------------------------------
r2DNIG <- function(n,alpha){
 alpha <- 2
 beta <- c(0,0)
 delta0 <- 0.55
 mu <- c(0,0)
 Lambda <- matrix(c(1,0,0,1),2,2)
 t(rNIG(n,alpha=alpha,beta=beta,delta=delta0,mu=mu,Lambda=Lambda))
}
# the next fake density plays no role in simulation
# but it is needed for model specification
d2DNIG <- function(n,alpha){ 
 rep(0,2)
}

## ------------------------------------------------------------------------
mod12 <- setPoisson(intensity="lambda", df=list("d2DNIG(z,)"),
 dim=2)
set.seed(123)
poisson12 <- simulate(mod12, true.par= list(lambda=1), 
 sampling=samp)
poisson12

## ------------------------------------------------------------------------
rMydis <- function(n,a=1){
 cbind(rnorm(n), rexp(n), rNIG(n,1,1,1,1))
}
dMydis <- function(n,a=1){
 rep(0,3)
}
mod13 <- setPoisson(intensity="lambda*t", 
 df=list("dMydis(z,1)"), dimension=3)
set.seed(123)
poisson13 <- simulate(mod13, true.par=list(lambda=5), 
sampling=samp)
poisson13

## ------------------------------------------------------------------------
mod14 <- setPoisson(intensity="alpha+lambda*t", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson14 <- simulate(mod14, sampling=samp,
true.par=list(alpha=1,lambda=.5,mu=0, sigma=2))
poisson14
fit14 <- qmle(poisson14, start=list(alpha=2,lambda=1,mu=0,sigma=1),
 lower=list(alpha=0.1, lambda=0.1,mu=-1,sigma=0.1), 
 upper=list(alpha=10,lambda=10,mu=3,sigma=4), 
 method="L-BFGS-B")
coef(fit14)

## ------------------------------------------------------------------------
summary(fit14)

## ----results='hide'------------------------------------------------------
mod15 <- setPoisson(intensity="lambda", 
 df=list("dNIG(z,alpha,beta,gamma,mu)"))
set.seed(123)
poisson15 <- simulate(mod15,sampling=samp, 
 true.par=list(lambda=10,alpha=2,beta=0.3,gamma=1,mu=0))
poisson15
fit15 <- qmle(poisson15, 
 start=list(beta=5,lambda=2,gamma=0.5,alpha=1,mu=0),
 lower=list(alpha=1,beta=0.1,lambda=0.1,gamma=0.1,mu=-1),
 upper=list(alpha=5,beta=0.99,lambda=20,gamma=2,mu=2), 
 method="L-BFGS-B")
summary(fit15)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit15)),width=60))

## ------------------------------------------------------------------------
mod16 <- setPoisson(intensity="beta*exp(-lambda*t)", 
 df=list("dexp(z,lambda)"))
set.seed(123)
poisson16 <- simulate(mod16, true.par=list(lambda=.2,beta=10),
 sampling=samp)
poisson16
fit16 <- qmle(poisson16, 
 start=list(beta=.5,lambda=2),
 lower=list(beta=0.1,lambda=0.1), 
 upper=list(beta=20,lambda=10), 
 method="L-BFGS-B")
summary(fit16)

## ----results='hide'------------------------------------------------------
mod17 <- setPoisson(intensity="lambda*t^(lambda-1)", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson17 <- simulate(mod17, true.par=list(lambda=2,mu=0, sigma=2),
 sampling=samp)
poisson17
fit17 <- qmle(poisson17, 
 start=list(lambda=5,mu=0,sigma=1),
 lower=list(lambda=0.1,mu=-1,sigma=0.1), 
 upper=list(lambda=10,mu=3,sigma=4), 
 method="L-BFGS-B")
summary(fit17)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit17)),width=60))

## ----results='hide'------------------------------------------------------
mod18 <- setPoisson(intensity="0.5*a*(1+cos(omega*t+phi))+lambda", 
 df=list("dnorm(z,mu,sigma)"))
set.seed(123)
poisson18 <- simulate(mod18, sampling=samp,
 true.par=list(a=2,omega=0.5,phi=3.14,lambda=5,mu=0,sigma=1))
 fit18 <- qmle(poisson18, 
 start=list(a=1, omega=0.2, phi=1, lambda=2, mu=1, sigma=2), 
 lower=list(a=0.1, omega=0.1, phi=0.1, lambda=0.1, mu=-2, sigma=0.1),
 upper=list(a=5, omega=1, phi=5, lambda=10, mu=2, sigma=3),
 method="L-BFGS-B")
summary(fit18)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(summary(fit18)),width=60))

