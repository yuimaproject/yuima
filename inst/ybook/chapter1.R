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
Rver <- paste(version$major,version$minor, collapse="",sep=".")
YUIMAver <- as.character(read.dcf(file=system.file("DESCRIPTION", package="yuima"),
    fields="Version"))

## ----eval=FALSE----------------------------------------------------------
## ybook(1)

## ----eval=FALSE----------------------------------------------------------
## ybook(3)

## ------------------------------------------------------------------------
data(cars)
class(cars)

## ------------------------------------------------------------------------
summary(cars)
mod <- lm(dist~speed, data=cars)
summary(mod)

## ------------------------------------------------------------------------
class(mod)

## ------------------------------------------------------------------------
methods(summary)

## ------------------------------------------------------------------------
x <- 1:4
x
class(x)
class(x) <- "lm"
class(x)

## ------------------------------------------------------------------------
require(stats4)
set.seed(123)
y <- rnorm(100, mean=1.5)
f <- function(theta=0) -sum(dnorm(x=y, mean=theta,log=TRUE))
fit <- mle(f)
fit

## ----results='hide'------------------------------------------------------
str(fit)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(fit)),width=60))

## ------------------------------------------------------------------------
fit@coef

## ------------------------------------------------------------------------
showMethods(summary)

## ----eval=FALSE----------------------------------------------------------
## install.packages("yuima")

## ----eval=FALSE----------------------------------------------------------
## install.packages("yuima",repos="http://R-Forge.R-project.org")

## ----eval=FALSE----------------------------------------------------------
## install.packages("yuima",repos="http://R-Forge.R-project.org",
##   type="source")

## ----echo=TRUE,eval=TRUE,results='markup',showWarnings=TRUE--------------
library(yuima)

## ----echo=TRUE,results='hide'--------------------------------------------
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)")

## ----results='hide'------------------------------------------------------
str(mod1)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(mod1)),width=60))

## ------------------------------------------------------------------------
mod1

## ----echo=TRUE,results='hide'--------------------------------------------
mod1b <- setModel(drift = "-3*s*y", diffusion = "1/(1+y^2)", 
state.var="y", time.var="s")

## ------------------------------------------------------------------------
tmp <- capture.output(str(mod1b))
writeLines(strwrap(tmp[c(2,3,4,17,19,23)],width=60))

## ----echo=TRUE,results='hide'--------------------------------------------
mod2 <- setModel(drift = "-mu*x", diffusion = "1/(1+x^gamma)")

## ------------------------------------------------------------------------
tmp <- capture.output(str(mod2)) 
writeLines(strwrap(tmp[c(2,3,4,9:13,17,19,23)],width=60))

## ------------------------------------------------------------------------
mod2

## ----sim-mod1,echo=TRUE,results='hide'-----------------------------------
set.seed(123)
X <- simulate(mod1)

## ----plot-mod1,echo=TRUE,fig.keep='none',results='hide'------------------
plot(X)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ----plot-mod1bis,echo=TRUE,fig.keep='none',results='hide'---------------
x0 <- 1
set.seed(123)
X <- simulate(mod1, xinit=x0)
plot(X)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1bis.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ----results='hide'------------------------------------------------------
str(X@data,vec.len=2)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(X@data,vec.len=2)),width=60))

## ------------------------------------------------------------------------
str(X@sampling,vec.len=2)

## ------------------------------------------------------------------------
tmp <- capture.output(str(X))
writeLines(strwrap(tmp[c(14:16,29,31,35,36)],width=60))

## ------------------------------------------------------------------------
X

## ----plot-mod1ter,echo=TRUE,fig.keep='none',results='hide'---------------
x0 <- 1
set.seed(123)
X <- simulate(mod1, xinit=x0, Initial=0.5, Terminal=1.2)
X
plot(X)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1ter.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ------------------------------------------------------------------------
str(X@sampling,vec.len=2)

## ----plot-mod1b,echo=TRUE,fig.keep='none',results='hide'-----------------
set.seed(123)
X <- simulate(mod1b, xinit=x0)
X
plot(X)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-mod1b.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ----sim-mod2,echo=TRUE,fig.keep='none',results='hide'-------------------
set.seed(123)
X <- simulate(mod2,true.param=list(mu=1,gamma=3))
plot(X)

## ----plot-mod2,echo=FALSE,results='hide'---------------------------------
pdf("figures/plot-mod2.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X)
dev.off()

## ------------------------------------------------------------------------
sol <- c("x1","x2") # variable for numerical solution
b <- c("-theta*x1","-x1-gamma*x2") # drift vector 
s <- matrix(c("1","x1","0","beta","x2","0"),2,3) # diff. mat.
mymod <- setModel(drift = b, diffusion = s, solve.variable = sol)

## ------------------------------------------------------------------------
samp <- setSampling(Terminal=3, n=3000)

## ----results='hide'------------------------------------------------------
str(samp)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(samp)),width=60))

## ------------------------------------------------------------------------
set.seed(123)
X2 <- simulate(mymod, sampling=samp, 
 true.param=list(theta=1,gamma=1,beta=1))
X2

## ----results='hide'------------------------------------------------------
str(X2@sampling)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(X2@sampling)),width=60))

## ----sub1----------------------------------------------------------------
newsamp <- setSampling(
random=list(rdist=c( function(x) rexp(x, rate=10), 
function(x) rexp(x, rate=20))) )

## ----results='hide'------------------------------------------------------
str(newsamp)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(newsamp)),width=60))

## ----sub2,echo=TRUE, fig.keep='none'-------------------------------------
newdata <- subsampling(X2, sampling=newsamp)
newdata
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green",pch=18)

## ----plot-sub2,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-sub2.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green",pch=18)
dev.off()

## ----sub3,echo=TRUE, fig.keep='none'-------------------------------------
newsamp <- setSampling(Terminal=3, delta=c(0.1,0.2), n=NULL)
newsamp
newdata <- subsampling(X2, sampling=newsamp)
newdata
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green", pch=18)

## ----plot-sub3,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-sub3.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green", pch=18)
dev.off()

## ------------------------------------------------------------------------
str(newdata@sampling)

## ----sub4,fig.keep='none'------------------------------------------------
set.seed(123)
Y.sub <- simulate(mymod,sampling=setSampling(delta=0.001,n=1000),
 subsampling=setSampling(delta=0.01,n=100), 
 true.par=list(theta=1,beta=1,gamma=1))
set.seed(123)
Y <- simulate(mymod, sampling=setSampling(delta=0.001,n=1000),
  true.par=list(theta=1,beta=1,gamma=1))
plot(Y, plot.type="single")
points(get.zoo.data(Y.sub)[[1]],col="red")
points(get.zoo.data(Y.sub)[[2]],col="green",pch=18)

## ----plot-sub4,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-sub4.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(Y, plot.type="single")
points(get.zoo.data(Y.sub)[[1]],col="red")
points(get.zoo.data(Y.sub)[[2]],col="green",pch=18)
dev.off()

## ----sub5,fig.keep='none'------------------------------------------------
plot(Y.sub, plot.type="single")

## ----plot-sub5,echo=FALSE, fig.keep='none',results='hide'----------------
pdf("figures/plot-sub5.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
plot(Y.sub, plot.type="single")
dev.off()

## ------------------------------------------------------------------------
Y
Y.sub

## ----eval=FALSE----------------------------------------------------------
## my.yuima <- setYuima(data=setData(X), model=mod)

## ----echo=TRUE,results='hide',tidy=TRUE----------------------------------
require(quantmod)
getSymbols("IBM", to = "2017-07-31")

## ----echo=TRUE-----------------------------------------------------------
str(IBM)

## ----echo=TRUE-----------------------------------------------------------
head(IBM)

## ------------------------------------------------------------------------
x <- setYuima(data=setData(IBM$IBM.Close))
str(x@data)

## ------------------------------------------------------------------------
y <- setYuima(data=setData(IBM$IBM.Close, delta=1/252))
str(y@data)

## ----setData,fig.keep='none'---------------------------------------------
plot(x, main="data with the original time stamps")
plot(y, main="time stamps of data rescaled")

## ----plot-setData,echo=FALSE, fig.keep='none',results='hide'-------------
pdf("figures/plot-setData.pdf",width=9,height=6)
par(mar=c(4,4,2,1), mfrow=c(2,1))
plot(x, main="data with the original time stamps")
plot(y, main="time stamps of data rescaled")
dev.off()

## ------------------------------------------------------------------------
x
y

## ----quantmod,message=FALSE----------------------------------------------
library(quantmod)
getSymbols("IBM", to = "2017-07-31")
attr(IBM, "src")

## ------------------------------------------------------------------------
getSymbols("IBM", to = "2017-07-31",  src="google")
attr(IBM, "src")

## ------------------------------------------------------------------------
getSymbols("DEXUSEU",src="FRED")
attr(DEXUSEU, "src")
getSymbols("EUR/USD",src="oanda")
attr(EURUSD, "src")
str(EURUSD)

## ----results='hide'------------------------------------------------------
library(tseries)
x <- get.hist.quote("IBM")
str(x)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(x)),width=60))

## ------------------------------------------------------------------------
mydat <- get.zoo.data(y)[[1]]
str(mydat)

## ------------------------------------------------------------------------
head(y@data@original.data)

## ----results='hide'------------------------------------------------------
str(y@data@original.data)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(y@data@original.data)),width=60))

## ------------------------------------------------------------------------
set.seed(123)
some.data <- rnorm(12)
str(some.data)

## ------------------------------------------------------------------------
X <- ts(some.data, frequency = 4, start = c(1961, 2))
X

## ------------------------------------------------------------------------
set.seed(123)
X <- ts(some.data, start = c(1964, 2), frequency = 12)
X

## ------------------------------------------------------------------------
time(X)[1:12]
deltat(X)
start(X)
end(X)
frequency(X)

## ------------------------------------------------------------------------
window(X, frequency=4)

## ------------------------------------------------------------------------
require(zoo)
X <- zoo( some.data )
X
str(X)

## ------------------------------------------------------------------------
index(X)

## ------------------------------------------------------------------------
rtimes <- cumsum(rexp(12,rate=0.2))
rtimes

## ------------------------------------------------------------------------
X <- zoo( rnorm(12), order.by = rtimes)
X
str(X)

## ------------------------------------------------------------------------
Xreg <- zooreg(some.data, start = c(1964, 2),  frequency = 12)
time(Xreg)

## ------------------------------------------------------------------------
Y <- as.ts(X)
time(X)
time(Y)

## ------------------------------------------------------------------------
require(xts)
my.time.stamps <- as.Date(rtimes)
my.time.stamps
X <- xts( some.data , order.by = my.time.stamps)
X
str(X)

## ------------------------------------------------------------------------
X.ts <- ts(some.data, start = c(1964, 2), frequency = 12)
X.ts
X.zoo <- as.zoo(X.ts)
X.zoo
X.xts <- as.xts(X.ts)
X.xts

## ----xts,fig.keep='none'-------------------------------------------------
plot(X)

## ----echo=FALSE,results='hide'-------------------------------------------
pdf("figures/plot-xts.pdf",width=9,height=4)
par(mar=c(4,4,2,1))
plot(X)
dev.off()

## ------------------------------------------------------------------------
require(tseries)
X <- irts( rtimes, some.data)
X
str(X)

## ------------------------------------------------------------------------
require(timeSeries)
X <- timeSeries( some.data,  my.time.stamps)
X
str(X)

## ------------------------------------------------------------------------
d <- ISOdate(2008,7,3)
d

## ----results='hide'------------------------------------------------------
args(ISOdate)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(args(ISOdate)),width=60))

## ------------------------------------------------------------------------
class(d)

## ------------------------------------------------------------------------
names(as.POSIXlt(d))
unlist(as.POSIXlt(d))

## ------------------------------------------------------------------------
format(d,"%a") # week day
format(d,"%A")

format(d,"%b") # month
format(d,"%B")

format(d,"%c") # full date

format(d,"%D") # yy/dd/mm
format(d,"%T") # hh:mm:ss

format(d,"%A %B  %d %H:%M:%S %Y")
format(d,"%A   %d/%m/%Y")
format(d,"%d/%m/%Y (%A)")

## ------------------------------------------------------------------------
x <- c("10jan1962", "2feb1970", "11jul2011", "27jun1968")
strptime(x, "%d%b%Y")

## ------------------------------------------------------------------------
Sys.getlocale()
Sys.setlocale("LC_ALL", "it_it")
strptime(x, "%d%b%Y")
Sys.setlocale("LC_ALL", "en_GB")
strptime(x, "%d%b%Y")

## ------------------------------------------------------------------------
format(ISOdate(2006,6,9),"%H:%M:%S")
format(as.POSIXct("2006-06-09"),"%H:%M:%S")

## ------------------------------------------------------------------------
holidayNYSE()
holidayNERC()

## ------------------------------------------------------------------------
ISOdate(2006,7,10) - ISOdate(2005, 3, 1)

## ------------------------------------------------------------------------
my.dates <- timeDate(c("2001-01-09", "2001-02-25")) 
diff(my.dates)

## ------------------------------------------------------------------------
listFinCenter("America*")[1:50]

## ------------------------------------------------------------------------
dA <- timeDate("2011-02-05", Fin="Europe/Zurich") 
dB <- timeDate("2016-01-22", Fin="America/Chicago") 
dA
dB

## ------------------------------------------------------------------------
set.seed(123) 
mydata <- rnorm(9) 
chardata <- sprintf("2010-0%s-01", 9:1) 
chardata

## ------------------------------------------------------------------------
X1 <- zoo(mydata, as.Date(chardata))
X2 <- xts(mydata, as.Date(chardata))
X3 <- timeSeries(mydata, chardata)

## ------------------------------------------------------------------------
X1
X2
X3

## ------------------------------------------------------------------------
zA <- zoo(mydata, as.POSIXct(chardata))
zB <- zoo(mydata, ISOdatetime(2016, 9:1, 1, 0,0,0))
zC <- zoo(mydata, ISOdate(2016, 9:1, 1, 0))
zA
zB
zC

## ------------------------------------------------------------------------
set.seed(123) 
val1 <- rnorm(9) 
val2 <- rnorm(6) 
mydate1 <- ISOdate(2016,1:9,1)
mydate2 <- ISOdate(2015,6:11,1)
Z1 <- zoo(val1, mydate1)
Z2 <- zoo(val2, mydate2)
rbind(Z1,Z2)
X1 <- xts(val1, mydate1)
X2 <- xts(val2, mydate2)
rbind(X1,X2)
W1 <- timeSeries(val1, mydate1)
W2 <- timeSeries(val2, mydate2)
rbind(W1,W2)

## ------------------------------------------------------------------------
mydate2 <- ISOdate(2016,4:9,1)
Z2 <- zoo(val2, mydate2)

## ----eval=FALSE----------------------------------------------------------
## rbind(Z1,Z2)

## ----echo=FALSE----------------------------------------------------------
cat(unclass(try(rbind(Z1,Z2))))

## ------------------------------------------------------------------------
X2 <- xts(val2, mydate2)
rbind(X1,X2)
W2 <- timeSeries(val2, mydate2)
rbind(W1,W2)

## ------------------------------------------------------------------------
merge(Z1,Z2)
merge(X1,X2)

## ------------------------------------------------------------------------
merge(W1,W2)

## ------------------------------------------------------------------------
W2 <- timeSeries(val2, mydate2, units="MyData")
merge(W1,W2)

## ------------------------------------------------------------------------
mydate1 <- ISOdate(2016,1:9,1)
mydate2 <- ISOdate(2015,6:11,1)
W1 <- timeSeries(val1, mydate1)
W2 <- timeSeries(val2, mydate2)

## ------------------------------------------------------------------------
rbind(W1,W2)

## ------------------------------------------------------------------------
rbind(W2,W1)

## ------------------------------------------------------------------------
sort( rbind(W2,W1) )
sort( rbind(W2,W1), decr=TRUE)

## ------------------------------------------------------------------------
W2
rev(W2)

## ----results="hide",message=FALSE----------------------------------------
require(sde)
data(quotes)
str(quotes)

## ----echo=FALSE----------------------------------------------------------
writeLines(strwrap(capture.output(str(quotes)),width=60))

## ------------------------------------------------------------------------
quotes[2,2:4]
quotes[10:20,"INTEL"]

## ------------------------------------------------------------------------
quotes$INTEL[10:20]

## ------------------------------------------------------------------------
mydate <- as.Date(sprintf("2006-08-%.2d",20:10))
mydate
quotes[mydate, 5:9]

## ------------------------------------------------------------------------
initial <- as.Date("2007-05-15")
terminal <- as.Date("2007-05-21")
quotes[ (time(quotes) >= initial) & (time(quotes)<= terminal), 4:9]

## ------------------------------------------------------------------------
getSymbols("IBM", from="2015-01-01", to = "2016-12-31")
str(IBM)

## ------------------------------------------------------------------------
IBM["2015-01","IBM.Close"]

## ------------------------------------------------------------------------
IBM["2016-02-11/2016-03-05","IBM.Close"]

## ------------------------------------------------------------------------
IBM["/2015-02-11","IBM.Close"]

## ------------------------------------------------------------------------
mod2 <- setModel(drift = "-mu*x", diffusion = "1/(1+x^gamma)")
mod2

## ------------------------------------------------------------------------
toLatex(mod2)

## ----echo=FALSE, results='asis', include=TRUE----------------------------
toLatex(mod2)

## ------------------------------------------------------------------------
sol <- c("x1","x2") # variable for numerical solution
b <- c("-theta*x1","-x1-gamma*x2") # drift vector 
s <- matrix(c("1","x1","0","delta","x2","0"),2,3) # diff. mat.
mymod <- setModel(drift = b, diffusion = s, solve.variable = sol)

## ----echo=FALSE, results='asis', include=TRUE----------------------------
toLatex(mymod)

## ----eval=FALSE----------------------------------------------------------
## install.packages("yuimaGUI")

## ----eval=FALSE----------------------------------------------------------
## library(yuimaGUI)
## yuimaGUI()

