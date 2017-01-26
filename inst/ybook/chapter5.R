### R code from vignette source '/Users/jago/Dropbox (VOICES)/yuima-book/chapter5.Rnw'

###################################################
### code chunk number 1: chapter5.Rnw:3-6
###################################################
options(width=60)
options(continue="  ")
require(yuima)


###################################################
### code chunk number 2: sim-mod4AB
###################################################
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


###################################################
### code chunk number 3: plot-mod4AB
###################################################
pdf("figures/plot-mod4AB.pdf",width=9,height=6)
par(mfrow=c(2,1))
par(mar=c(2,3,1,1))
plot(X1,main="H=0.3")
plot(X2,main="H=0.7")
dev.off()


###################################################
### code chunk number 4: chapter5.Rnw:56-57
###################################################
str(mod4A)


###################################################
### code chunk number 5: chapter5.Rnw:222-229
###################################################
set.seed(123)
samp <- setSampling(Terminal=100, n=10000)
mod <- setModel(drift="-lambda*x", diffusion="sigma", hurst=NA)
ou <- setYuima(model=mod, sampling=samp)
fou <- simulate(ou, xinit=1, 
  true.param=list(lambda=2, sigma=1), hurst=0.7)
fou


###################################################
### code chunk number 6: chapter5.Rnw:234-235
###################################################
qgv(fou)


###################################################
### code chunk number 7: chapter5.Rnw:246-247
###################################################
mmfrac(fou)


###################################################
### code chunk number 8: chapter5.Rnw:253-255
###################################################
data(MWK151)
str(MWK151)


###################################################
### code chunk number 9: MWK151
###################################################
par(mfrow=c(1,2))
plot(MWK151, main="Methuselah Walk ring widths", xlab="year")
plot(acf(MWK151))


###################################################
### code chunk number 10: plot-MWK151
###################################################
pdf("figures/plot-MWK151.pdf",width=9,height=6)
par(mar=c(4,4,2,0))
par(mfrow=c(1,2))
plot(MWK151, main="Methuselah Walk ring widths", xlab="year")
plot(acf(MWK151))
dev.off()


###################################################
### code chunk number 11: chapter5.Rnw:277-281
###################################################
mod <- setModel(drift="-lambda *x", diffusion="sigma", hurst=NA)
mwk  <- setYuima(model=mod, data=setData(MWK151))
mwk
mmfrac(mwk)


