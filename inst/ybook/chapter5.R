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

## ----sim-mod4AB, echo=TRUE,fig.keep='none'-------------------------------
mod4A <- setModel(drift="3*y", diffusion=1, hurst=0.3, solve.var="y")
mod4A
mod4B <- setModel(drift="3*y", diffusion=1, hurst=0.7, solve.var="y")
mod4B
set.seed(123)
X1 <- simulate(mod4A,sampling=setSampling(n=1000))
X2 <- simulate(mod4B,sampling=setSampling(n=1000))
par(mfrow=c(2,1))
par(mar=c(2,3,1,1))
plot(X1,main="H=0.3")
plot(X2,main="H=0.7")

## ----results='hide'------------------------------------------------------
str(mod4A)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(mod4A)),width=60))

## ----plot-mod4AB,echo=FALSE,results='hide'-------------------------------
pdf("figures/plot-mod4AB.pdf",width=9,height=4)
par(mfrow=c(2,1))
par(mar=c(2,3,1,1))
plot(X1,main="H=0.3")
plot(X2,main="H=0.7")
dev.off()

## ------------------------------------------------------------------------
set.seed(123)
samp <- setSampling(Terminal=100, n=10000)
mod <- setModel(drift="-lambda*x", diffusion="sigma", hurst=NA)
ou <- setYuima(model=mod, sampling=samp)
fou <- simulate(ou, xinit=1, 
  true.param=list(lambda=2, sigma=1), hurst=0.7)
fou

## ------------------------------------------------------------------------
qgv(fou)

## ------------------------------------------------------------------------
mmfrac(fou)

## ----results='hide'------------------------------------------------------
data(MWK151)
str(MWK151)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(MWK151)),width=60))

## ----MWK151, echo=TRUE,fig.keep='none'-----------------------------------
par(mfrow=c(1,2))
plot(MWK151, main="Methuselah Walk ring widths", xlab="year")
plot(acf(MWK151))

## ----plot-MWK151,echo=FALSE,results='hide'-------------------------------
pdf("figures/plot-MWK151.pdf",width=9,height=6)
par(mar=c(4,4,2,1))
par(mfrow=c(1,2))
plot(MWK151, main="Methuselah Walk ring widths", xlab="year")
plot(acf(MWK151))
dev.off()

## ------------------------------------------------------------------------
mod <- setModel(drift="-lambda *x", diffusion="sigma", hurst=NA)
mwk  <- setYuima(model=mod, data=setData(MWK151))
mwk
mmfrac(mwk)

