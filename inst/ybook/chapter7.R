### R code from vignette source '/Users/jago/Dropbox (VOICES)/yuima-book/chapter7.Rnw'

###################################################
### code chunk number 1: chapter7.Rnw:3-6
###################################################
options(width=60)
options(continue="  ")
require(yuima)


###################################################
### code chunk number 2: setCogarch
###################################################
# COGARCH(1,1) driven by CP
Cog11 <- setCogarch(p = 1, q=1,  measure = list(intensity="1", 
  df="dnorm(z, 0, 1)"), measure.type = "CP", XinExpr = TRUE)
Cog11

# COGARCH(2,2) driven by CP
Cog22 <- setCogarch(p=2, q=2, measure = list(intensity="1", 
  df="dnorm(z, 0, 1)"), measure.type = "CP", XinExpr = TRUE)
Cog22


###################################################
### code chunk number 3: str-cog
###################################################
class(Cog11)
slotNames(Cog11)
str(Cog11@info,2)


###################################################
### code chunk number 4: chapter7.Rnw:215-226
###################################################
# Param of the COGARCH(1,1)
paramCP11 <- list(a1 = 0.038, b1=  0.053, a0 = 0.04/0.053, 
 y01 = 50.31)
check11 <- Diagnostic.Cogarch(Cog11, param=paramCP11)
str(check11)

# Param of the COGARCH(2,2)
paramCP22 <- list(a1 = 0.04, a2 = 0.001, b1 = 0.705, b2 = 0.1, 
 a0 = 0.1, y01=01, y02 = 0)
check22 <- Diagnostic.Cogarch(Cog22, param=paramCP22)
str(check22)


###################################################
### code chunk number 5: cogarch-euler
###################################################
model1 <- setCogarch(p = 1, q = 1, 
  measure=list("rvgamma(z, 1, sqrt(2), 0, 0)"), 
  measure.type = "code", Cogarch.var = "G", 
  V.var = "v", Latent.var="x", XinExpr=TRUE)


###################################################
### code chunk number 6: cogarch-euler-bad
###################################################
param1 <- list(a1 = 0.038, b1 = 301, a0 =0.01, x01 = 0)
Diagnostic.Cogarch(model1, param=param1)
Terminal1 <- 5
n1 <- 750
samp1 <- setSampling(Terminal=Terminal1, n=n1)
set.seed(123)
sim1 <- simulate(model1, sampling = samp1, true.parameter = param1,  
 method="euler")


###################################################
### code chunk number 7: plot-cogarch
###################################################
plot(sim1, main="VG-COGARCH(1,1) model with Euler scheme")


###################################################
### code chunk number 8: plot-cogarch1
###################################################
pdf("figures/plot-cogarch1.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(sim1, 
 main="VG-COGARCH(1,1) model with Euler scheme")
dev.off()


###################################################
### code chunk number 9: sim-cogarch2
###################################################
set.seed(123)
sim2 <- simulate(model1, sampling = samp1, true.parameter = param1,
        method="mixed")


###################################################
### code chunk number 10: plot-cogach2
###################################################
plot(sim2, main="VG-COGARCH(1,1) model with mixed scheme")


###################################################
### code chunk number 11: plot-cogarch2
###################################################
pdf("figures/plot-cogarch2.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(sim2, 
 main="VG-COGARCH(1,1) model with mixed scheme")
dev.off()


###################################################
### code chunk number 12: chapter7.Rnw:342-345
###################################################
sampCP <- setSampling(0, 1000, 5000)
simCog11 <- simulate(Cog11, true.par=paramCP11, sampling=sampCP)
simCog22 <- simulate(Cog22, true.par=paramCP22, sampling=sampCP)


###################################################
### code chunk number 13: plot-cogachs
###################################################
plot(simCog11, main="CP-COGARCH(1,1) with Gaussian noise")
plot(simCog22, main="CP-COGARCH(2,2) with Gaussian noise")


###################################################
### code chunk number 14: plot-cogarchs2
###################################################
pdf("figures/plot-cogarchs1.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(simCog11, main="CP-COGARCH(1,1) with Gaussian noise")
dev.off()
pdf("figures/plot-cogarchs2.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(simCog22, main="CP-COGARCH(2,2) with Gaussian noise")
dev.off()


###################################################
### code chunk number 15: gmm-cogarch
###################################################
set.seed(123)
sampCP <- setSampling(0, 5000, 15000)
simCog11 <- simulate(Cog11, true.par=paramCP11, sampling=sampCP)
fit11 <- gmm(simCog11, start=paramCP11)
summary(fit11)
mat <- rbind(coef(fit11), unlist(paramCP11[names(coef(fit11))]))
rownames(mat) <- c("gmm", "true")
mat


###################################################
### code chunk number 16: est-cogarch3
###################################################
param.VG <- list(a1 = 0.038,  b1 =  0.053, a0 = 0.04 / 0.053,
  y01 = 50.33)
cog.VG <- setCogarch(p = 1, q = 1, work = FALSE,
  measure=list("rvgamma(z, 1, sqrt(2), 0, 0)"),
  measure.type = "code", XinExpr = TRUE)
samp.VG <- setSampling(Terminal = 1000, n = 15000)
set.seed(123)
sim.VG <- simulate(cog.VG, true.parameter = param.VG,
  sampling = samp.VG, method = "mixed")
fit.gmm <- gmm(sim.VG, start=param.VG)
fit.qmle <- qmle(sim.VG, start=param.VG, grideq=TRUE)
nm <- names(coef(fit.gmm))
mat <- rbind(coef(fit.gmm), coef(fit.qmle)[nm],
 unlist(param.VG[nm]))
rownames(mat) <- c("gmm", "qmle", "true")
round(mat,5)


###################################################
### code chunk number 17: chapter7.Rnw:579-593
###################################################
require(quantmod)
getSymbols("NXT.L",  to="2016-12-31")
S <- NXT.L$NXT.L.Close
X <- na.omit(diff(log(S)))
mX <- mean(X)
X <- X - mX
plot(X,  main="Log-returns of NEXT Plc")
require(rugarch)
spec <- ugarchspec(variance.model = 
 list(model = "sGARCH", garchOrder = c(1, 1)),
 mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
fitGARCH <- ugarchfit(data = X, spec = spec)
GARCH11param <- coef(fitGARCH)
GARCH11param


###################################################
### code chunk number 18: plot-nextplc
###################################################
pdf("figures/plot-nextplc.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(X,  main="Log-returns of NEXT Plc")
dev.off()


###################################################
### code chunk number 19: chapter7.Rnw:608-622
###################################################
Delta <- 1/252
ParGarToCog<- function(GARCH11param, dt, names=NULL){
    if(is.null(names))
     names <- names(GARCH11param)
    my.omega <- GARCH11param["omega"]
    my.alpha <- GARCH11param["alpha1"]
    my.beta <- GARCH11param["beta1"]
    a1 <- my.alpha/dt
    b1 <- -log(my.beta)/dt
    a0 <- my.omega/(b1*dt^2)
    qmleparInGARCH <- c(a0, a1, b1)
    names(qmleparInGARCH) <- c("a0", "a1", "b1")
    return(qmleparInGARCH)
}


###################################################
### code chunk number 20: chapter7.Rnw:625-637
###################################################
ParGarToCog(GARCH11param, Delta)
start <- as.list(ParGarToCog(GARCH11param, Delta))

modCog11 <- setCogarch(p=1, q=1, measure =
 list(intensity="1", df=list("dnorm(z, 0, 1)")), measure.type="CP")
NXT.data <- setData(cumsum(X), delta = Delta)
Cog11 <- setYuima(data = NXT.data, model = modCog11)
Cog11.fit <- qmle(yuima = Cog11, grideq=TRUE, 
 start = c(start, y1 = 0.1),
 aggregation = FALSE, method = "Nelder-Mead")
COGARCH11par <- coef(Cog11.fit)
COGARCH11par


###################################################
### code chunk number 21: chapter7.Rnw:640-653
###################################################
ParCogToGar<- function(COGARCH11param, dt, names=NULL){
    a0 <- COGARCH11param["a0"]
    a1 <- COGARCH11param["a1"]
    b1 <- COGARCH11param["b1"]
    my.omega <- a0*b1*dt^2
    my.alpha <- a1*dt
    my.beta <- exp(-b1*dt)
    qmleparInGARCH <- c(my.omega, my.alpha, my.beta)
    names(qmleparInGARCH) <- c("omega", "alpha1", "beta1")
    return(qmleparInGARCH)
}
ParCogToGar(COGARCH11par, Delta)
GARCH11param


