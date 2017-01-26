### R code from vignette source '/Users/jago/Dropbox (VOICES)/yuima-book/chapter1.Rnw'

###################################################
### code chunk number 1: chapter1.Rnw:3-6
###################################################
options(width=60)
options(continue="  ")
require(yuima)


###################################################
### code chunk number 2: chapter1.Rnw:110-112
###################################################
data(cars)
class(cars)


###################################################
### code chunk number 3: chapter1.Rnw:115-118
###################################################
summary(cars)
mod <- lm(dist~speed, data=cars)
summary(mod)


###################################################
### code chunk number 4: chapter1.Rnw:121-122
###################################################
class(mod)


###################################################
### code chunk number 5: chapter1.Rnw:125-126
###################################################
methods(summary)


###################################################
### code chunk number 6: chapter1.Rnw:130-135
###################################################
x <- 1:4
x
class(x)
class(x) <- "lm"
class(x)


###################################################
### code chunk number 7: chapter1.Rnw:140-146
###################################################
require(stats4)
set.seed(123)
y <- rnorm(100, mean=1.5)
f <- function(theta=0) -sum(dnorm(x=y, mean=theta,log=TRUE))
fit <- mle(f)
fit


###################################################
### code chunk number 8: chapter1.Rnw:149-150
###################################################
str(fit)


###################################################
### code chunk number 9: chapter1.Rnw:153-154
###################################################
fit@coef


###################################################
### code chunk number 10: chapter1.Rnw:157-158
###################################################
showMethods(summary)


###################################################
### code chunk number 11: chapter1.Rnw:165-166 (eval = FALSE)
###################################################
## install.packages("yuima")


###################################################
### code chunk number 12: chapter1.Rnw:176-177 (eval = FALSE)
###################################################
## install.packages("yuima",repos="http://R-Forge.R-project.org")


###################################################
### code chunk number 13: chapter1.Rnw:181-183 (eval = FALSE)
###################################################
## install.packages("yuima",repos="http://R-Forge.R-project.org", 
##   type="source")


###################################################
### code chunk number 14: chapter1.Rnw:188-189
###################################################
library(yuima)


###################################################
### code chunk number 15: chapter1.Rnw:288-289
###################################################
mod1 <- setModel(drift = "-3*x", diffusion = "1/(1+x^2)")


###################################################
### code chunk number 16: chapter1.Rnw:295-296
###################################################
str(mod1)


###################################################
### code chunk number 17: chapter1.Rnw:302-303
###################################################
mod1


###################################################
### code chunk number 18: chapter1.Rnw:315-317
###################################################
mod1b <- setModel(drift = "-3*s*y", diffusion = "1/(1+y^2)", 
state.var="y", time.var="s")


###################################################
### code chunk number 19: chapter1.Rnw:320-322
###################################################
tmp <- capture.output(str(mod1b))
cat(tmp[c(2,3,4,17,19,23)],sep="\n")


###################################################
### code chunk number 20: chapter1.Rnw:334-335
###################################################
mod2 <- setModel(drift = "-mu*x", diffusion = "1/(1+x^gamma)")


###################################################
### code chunk number 21: chapter1.Rnw:340-342
###################################################
capture.output(str(mod2)) -> tmp
cat(tmp[c(2,3,4,9:13,17,19,23)],sep="\n")


###################################################
### code chunk number 22: chapter1.Rnw:345-346
###################################################
mod2


###################################################
### code chunk number 23: sim-mod1
###################################################
set.seed(123)
X <- simulate(mod1)


###################################################
### code chunk number 24: plot-mod1
###################################################
plot(X)


###################################################
### code chunk number 25: chapter1.Rnw:370-374
###################################################
pdf("figures/plot-mod1.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X)
dev.off()


###################################################
### code chunk number 26: plot-mod1bis
###################################################
x0 <- 1
set.seed(123)
X <- simulate(mod1, xinit=x0)
plot(X)


###################################################
### code chunk number 27: chapter1.Rnw:390-394
###################################################
pdf("figures/plot-mod1bis.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X)
dev.off()


###################################################
### code chunk number 28: chapter1.Rnw:401-402
###################################################
str(X@data,vec.len=2)


###################################################
### code chunk number 29: chapter1.Rnw:405-406
###################################################
str(X@sampling,vec.len=2)


###################################################
### code chunk number 30: chapter1.Rnw:411-413
###################################################
tmp <- capture.output(str(X))
cat(tmp[c(14:16,29,31,35,36)],sep="\n")


###################################################
### code chunk number 31: chapter1.Rnw:416-417
###################################################
X


###################################################
### code chunk number 32: plot-mod1ter
###################################################
x0 <- 1
set.seed(123)
X <- simulate(mod1, xinit=x0, Initial=0.5, Terminal=1.2)
X
plot(X)


###################################################
### code chunk number 33: chapter1.Rnw:428-432
###################################################
pdf("figures/plot-mod1ter.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X)
dev.off()


###################################################
### code chunk number 34: chapter1.Rnw:439-440
###################################################
str(X@sampling,vec.len=2)


###################################################
### code chunk number 35: plot-mod1b
###################################################
set.seed(123)
X <- simulate(mod1b, xinit=x0)
X
plot(X)


###################################################
### code chunk number 36: chapter1.Rnw:458-462
###################################################
pdf("figures/plot-mod1b.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X)
dev.off()


###################################################
### code chunk number 37: sim-mod2
###################################################
set.seed(123)
X <- simulate(mod2,true.param=list(mu=1,gamma=3))
plot(X)


###################################################
### code chunk number 38: plot-mod2
###################################################
pdf("figures/plot-mod2.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X)
dev.off()


###################################################
### code chunk number 39: setModel1
###################################################
sol <- c("x1","x2") # variable for numerical solution
b <- c("-theta*x1","-x1-gamma*x2")   # drift vector 
s <- matrix(c("1","x1","0","delta","x2","0"),2,3)  #  diffusion matrix
mymod <- setModel(drift = b, diffusion = s, solve.variable = sol)


###################################################
### code chunk number 40: chapter1.Rnw:551-552
###################################################
samp <- setSampling(Terminal=3, n=3000)


###################################################
### code chunk number 41: chapter1.Rnw:555-556
###################################################
str(samp)


###################################################
### code chunk number 42: chapter1.Rnw:560-564
###################################################
set.seed(123)
X2 <- simulate(mymod, sampling=samp, 
 true.param=list(theta=1,gamma=1,delta=1))
X2


###################################################
### code chunk number 43: chapter1.Rnw:567-568
###################################################
str(X2@sampling)


###################################################
### code chunk number 44: sub1
###################################################
newsamp <- setSampling(
random=list(rdist=c( function(x) rexp(x, rate=10), 
function(x) rexp(x, rate=20))) )
str(newsamp)


###################################################
### code chunk number 45: sub2
###################################################
newdata <- subsampling(X2, sampling=newsamp)
newdata
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green",pch=18)


###################################################
### code chunk number 46: plot-sub2
###################################################
pdf("figures/plot-sub2.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green",pch=18)
dev.off()


###################################################
### code chunk number 47: sub3
###################################################
newsamp <- setSampling(Terminal=3, delta=c(0.1,0.2), n=NULL)
newsamp
newdata <- subsampling(X2, sampling=newsamp)
newdata
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green", pch=18)


###################################################
### code chunk number 48: plot-sub3
###################################################
pdf("figures/plot-sub3.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(X2,plot.type="single", lty=c(1,3),ylab="X2")
points(get.zoo.data(newdata)[[1]],col="red")
points(get.zoo.data(newdata)[[2]],col="green", pch=18)
dev.off()


###################################################
### code chunk number 49: chapter1.Rnw:629-630
###################################################
str(newdata@sampling)


###################################################
### code chunk number 50: sub4
###################################################
set.seed(123)
Y.sub <- simulate(mymod,sampling=setSampling(delta=0.001,n=1000),
 subsampling=setSampling(delta=0.01,n=100), 
 true.par=list(theta=1,delta=1,gamma=1))
set.seed(123)
Y <- simulate(mymod, sampling=setSampling(delta=0.001,n=1000),
  true.par=list(theta=1,delta=1,gamma=1))
plot(Y, plot.type="single")
points(get.zoo.data(Y.sub)[[1]],col="red")
points(get.zoo.data(Y.sub)[[2]],col="green",pch=18)


###################################################
### code chunk number 51: plot-sub4
###################################################
pdf("figures/plot-sub4.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(Y, plot.type="single")
points(get.zoo.data(Y.sub)[[1]],col="red")
points(get.zoo.data(Y.sub)[[2]],col="green",pch=18)
dev.off()


###################################################
### code chunk number 52: sub5
###################################################
plot(Y.sub, plot.type="single")


###################################################
### code chunk number 53: plot-sub5
###################################################
pdf("figures/plot-sub5.pdf",width=9,height=4)
par(mar=c(4,4,0,0))
plot(Y.sub, plot.type="single")
dev.off()


###################################################
### code chunk number 54: chapter1.Rnw:675-677
###################################################
Y
Y.sub


###################################################
### code chunk number 55: chapter1.Rnw:681-682 (eval = FALSE)
###################################################
## my.yuima <- setYuima(data=setData(X), model=mod)


###################################################
### code chunk number 56: chapter1.Rnw:688-690
###################################################
data <- read.csv("http://chart.yahoo.com/table.csv?s=IBM&g=d&x=.csv", 
 stringsAsFactor=FALSE)


###################################################
### code chunk number 57: chapter1.Rnw:694-695
###################################################
head(data)


###################################################
### code chunk number 58: chapter1.Rnw:698-700
###################################################
tmp <- zoo(data$Close,order.by=as.Date(data$Date, format="%Y-%m-%d"))
str(tmp)


###################################################
### code chunk number 59: chapter1.Rnw:704-706
###################################################
x <- setYuima(data=setData(tmp))
str(x@data)


###################################################
### code chunk number 60: chapter1.Rnw:713-715
###################################################
y <- setYuima(data=setData(tmp, delta=1/252))
str(y@data)


###################################################
### code chunk number 61: setData
###################################################
plot(x, main="data with the original time stamps")
plot(y, main="time stamps of data rescaled")


###################################################
### code chunk number 62: plot-setData
###################################################
pdf("figures/plot-setData.pdf",width=9,height=6)
par(mar=c(4,4,2,0), mfrow=c(2,1))
plot(x, main="data with the original time stamps")
plot(y, main="time stamps of data rescaled")
dev.off()


###################################################
### code chunk number 63: chapter1.Rnw:736-738
###################################################
x
y


###################################################
### code chunk number 64: quantmod
###################################################
library(quantmod)
getSymbols("IBM")
attr(IBM, "src")


###################################################
### code chunk number 65: chapter1.Rnw:753-755
###################################################
getSymbols("IBM",src="google")
attr(IBM, "src")


###################################################
### code chunk number 66: chapter1.Rnw:758-763
###################################################
getSymbols("DEXUSEU",src="FRED")
attr(DEXUSEU, "src")
getSymbols("EUR/USD",src="oanda")
attr(EURUSD, "src")
str(EURUSD)


###################################################
### code chunk number 67: chapter1.Rnw:768-771
###################################################
require(fImport)
X <- yahooSeries("IBM")
str(X)


###################################################
### code chunk number 68: chapter1.Rnw:774-776
###################################################
X <- yahooImport("IBM")
str(X)


###################################################
### code chunk number 69: chapter1.Rnw:783-789
###################################################
library(tseries)
x <- get.hist.quote("IBM")
str(x)
x <- get.hist.quote(instrument = "EUR/USD", provider = "oanda", 
 start = Sys.Date() - 300)
str(x)


###################################################
### code chunk number 70: chapter1.Rnw:801-803
###################################################
mydat <- get.zoo.data(y)[[1]]
str(mydat)


###################################################
### code chunk number 71: chapter1.Rnw:806-808
###################################################
str(y@data@original.data)
head(y@data@original.data)


###################################################
### code chunk number 72: chapter1.Rnw:826-829
###################################################
set.seed(123)
some.data <- rnorm(12)
str(some.data)


###################################################
### code chunk number 73: chapter1.Rnw:832-834
###################################################
X <- ts(some.data, frequency = 4, start = c(1961, 2))
X


###################################################
### code chunk number 74: chapter1.Rnw:837-840
###################################################
set.seed(123)
X <- ts(some.data, start = c(1964, 2), frequency = 12)
X


###################################################
### code chunk number 75: chapter1.Rnw:845-850
###################################################
time(X)[1:12]
deltat(X)
start(X)
end(X)
frequency(X)


###################################################
### code chunk number 76: chapter1.Rnw:853-854
###################################################
window(X, frequency=4)


###################################################
### code chunk number 77: chapter1.Rnw:860-864
###################################################
require(zoo)
X <- zoo( some.data )
X
str(X)


###################################################
### code chunk number 78: chapter1.Rnw:867-868
###################################################
index(X)


###################################################
### code chunk number 79: chapter1.Rnw:871-873
###################################################
rtimes <- cumsum(rexp(12,rate=0.2))
rtimes


###################################################
### code chunk number 80: chapter1.Rnw:876-879
###################################################
X <- zoo( rnorm(12), order.by = rtimes)
X
str(X)


###################################################
### code chunk number 81: chapter1.Rnw:882-884
###################################################
Xreg <- zooreg(some.data, start = c(1964, 2),  frequency = 12)
time(Xreg)[1:12]


###################################################
### code chunk number 82: chapter1.Rnw:887-890
###################################################
Y <- as.ts(X)
time(X)
time(Y)


###################################################
### code chunk number 83: chapter1.Rnw:894-900
###################################################
require(xts)
my.time.stamps <- as.Date(rtimes)
my.time.stamps
X <- xts( some.data , order.by = my.time.stamps)
X
str(X)


###################################################
### code chunk number 84: chapter1.Rnw:905-911
###################################################
X.ts <- ts(some.data, start = c(1964, 2), frequency = 12)
X.ts
X.zoo <- as.zoo(X.ts)
X.zoo
X.xts <- as.xts(X.ts)
X.xts


###################################################
### code chunk number 85: xts
###################################################
plot(X)


###################################################
### code chunk number 86: chapter1.Rnw:919-923
###################################################
pdf("figures/plot-xts.pdf",width=9,height=4)
par(mar=c(4,4,2,0))
plot(X)
dev.off()


###################################################
### code chunk number 87: chapter1.Rnw:936-940
###################################################
require(tseries)
X <- irts( rtimes, some.data)
X
str(X)


###################################################
### code chunk number 88: chapter1.Rnw:946-950
###################################################
require(timeSeries)
X <- timeSeries( some.data,  my.time.stamps)
X
str(X)


###################################################
### code chunk number 89: chapter1.Rnw:959-961
###################################################
d <- ISOdate(2008,7,3)
d


###################################################
### code chunk number 90: chapter1.Rnw:964-965
###################################################
args(ISOdate)


###################################################
### code chunk number 91: chapter1.Rnw:969-970
###################################################
class(d)


###################################################
### code chunk number 92: chapter1.Rnw:973-975
###################################################
names(as.POSIXlt(d))
unlist(as.POSIXlt(d))


###################################################
### code chunk number 93: chapter1.Rnw:979-993
###################################################
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


###################################################
### code chunk number 94: chapter1.Rnw:997-999
###################################################
x <- c("10jan1962", "2feb1970", "11jul2011", "27jun1968")
strptime(x, "%d%b%Y")


###################################################
### code chunk number 95: chapter1.Rnw:1002-1007
###################################################
Sys.getlocale()
Sys.setlocale("LC_ALL", "it_it")
strptime(x, "%d%b%Y")
Sys.setlocale("LC_ALL", "en_GB")
strptime(x, "%d%b%Y")


###################################################
### code chunk number 96: chapter1.Rnw:1010-1012
###################################################
format(ISOdate(2006,6,9),"%H:%M:%S")
format(as.POSIXct("2006-06-09"),"%H:%M:%S")


###################################################
### code chunk number 97: chapter1.Rnw:1016-1018
###################################################
holidayNYSE()
holidayNERC()


###################################################
### code chunk number 98: chapter1.Rnw:1021-1022
###################################################
ISOdate(2006,7,10) - ISOdate(2005, 3, 1)


###################################################
### code chunk number 99: chapter1.Rnw:1025-1027
###################################################
my.dates <- timeDate(c("2001-01-09", "2001-02-25")) 
diff(my.dates)


###################################################
### code chunk number 100: chapter1.Rnw:1030-1031
###################################################
listFinCenter("America*")


###################################################
### code chunk number 101: chapter1.Rnw:1034-1038
###################################################
dA <- timeDate("2011-02-05", Fin="Europe/Zurich") 
dB <- timeDate("2016-01-22", Fin="America/Chicago") 
dA
dB


###################################################
### code chunk number 102: chapter1.Rnw:1046-1050
###################################################
set.seed(123) 
mydata <- rnorm(9) 
chardata <- sprintf("2010-0%s-01", 9:1) 
chardata


###################################################
### code chunk number 103: chapter1.Rnw:1053-1056
###################################################
X1 <- zoo(mydata, as.Date(chardata))
X2 <- xts(mydata, as.Date(chardata))
X3 <- timeSeries(mydata, chardata)


###################################################
### code chunk number 104: chapter1.Rnw:1059-1062
###################################################
X1
X2
X3


###################################################
### code chunk number 105: chapter1.Rnw:1065-1071
###################################################
zA <- zoo(mydata, as.POSIXct(chardata))
zB <- zoo(mydata, ISOdatetime(2016, 9:1, 1, 0,0,0))
zC <- zoo(mydata, ISOdate(2016, 9:1, 1, 0))
zA
zB
zC


###################################################
### code chunk number 106: chapter1.Rnw:1076-1090
###################################################
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


###################################################
### code chunk number 107: chapter1.Rnw:1093-1095
###################################################
mydate2 <- ISOdate(2016,4:9,1)
Z2 <- zoo(val2, mydate2)


###################################################
### code chunk number 108: chapter1.Rnw:1098-1099 (eval = FALSE)
###################################################
## rbind(Z1,Z2)


###################################################
### code chunk number 109: chapter1.Rnw:1101-1102
###################################################
cat(unclass(try(rbind(Z1,Z2))))


###################################################
### code chunk number 110: chapter1.Rnw:1105-1109
###################################################
X2 <- xts(val2, mydate2)
rbind(X1,X2)
W2 <- timeSeries(val2, mydate2)
rbind(W1,W2)


###################################################
### code chunk number 111: chapter1.Rnw:1113-1115
###################################################
merge(Z1,Z2)
merge(X1,X2)


###################################################
### code chunk number 112: chapter1.Rnw:1118-1119
###################################################
merge(W1,W2)


###################################################
### code chunk number 113: chapter1.Rnw:1122-1124
###################################################
W2 <- timeSeries(val2, mydate2, units="MyData")
merge(W1,W2)


###################################################
### code chunk number 114: chapter1.Rnw:1127-1131
###################################################
mydate1 <- ISOdate(2016,1:9,1)
mydate2 <- ISOdate(2015,6:11,1)
W1 <- timeSeries(val1, mydate1)
W2 <- timeSeries(val2, mydate2)


###################################################
### code chunk number 115: chapter1.Rnw:1134-1135
###################################################
rbind(W1,W2)


###################################################
### code chunk number 116: chapter1.Rnw:1138-1139
###################################################
rbind(W2,W1)


###################################################
### code chunk number 117: chapter1.Rnw:1142-1144
###################################################
sort( rbind(W2,W1) )
sort( rbind(W2,W1), decr=TRUE)


###################################################
### code chunk number 118: chapter1.Rnw:1147-1149
###################################################
W2
rev(W2)


###################################################
### code chunk number 119: chapter1.Rnw:1155-1158
###################################################
require(sde)
data(quotes)
str(quotes)


###################################################
### code chunk number 120: chapter1.Rnw:1162-1164
###################################################
quotes[2,2:4]
quotes[10:20,"INTEL"]


###################################################
### code chunk number 121: chapter1.Rnw:1167-1168
###################################################
quotes$INTEL[10:20]


###################################################
### code chunk number 122: chapter1.Rnw:1171-1174
###################################################
mydate <- as.Date(sprintf("2006-08-%.2d",20:10))
mydate
quotes[mydate, 5:9]


###################################################
### code chunk number 123: chapter1.Rnw:1178-1181
###################################################
initial <- as.Date("2007-05-15")
terminal <- as.Date("2007-05-21")
quotes[ (time(quotes) >= initial) & (time(quotes)<= terminal), 4:9]


###################################################
### code chunk number 124: chapter1.Rnw:1188-1190
###################################################
mod2 <- setModel(drift = "-mu*x", diffusion = "1/(1+x^gamma)")
mod2


###################################################
### code chunk number 125: chapter1.Rnw:1193-1194
###################################################
toLatex(mod2)


###################################################
### code chunk number 126: chapter1.Rnw:1197-1198
###################################################
toLatex(mod2)


###################################################
### code chunk number 127: chapter1.Rnw:1201-1205
###################################################
sol <- c("x1","x2") # variable for numerical solution
b <- c("-theta*x1","-x1-gamma*x2")   # drift vector 
s <- matrix(c("1","x1","0","delta","x2","0"),2,3)  #  diff. matrix
mymod <- setModel(drift = b, diffusion = s, solve.variable = sol)


###################################################
### code chunk number 128: chapter1.Rnw:1208-1209
###################################################
toLatex(mymod)


###################################################
### code chunk number 129: chapter1.Rnw:1217-1218 (eval = FALSE)
###################################################
## install.packages("yuimaGUI")


###################################################
### code chunk number 130: chapter1.Rnw:1222-1224 (eval = FALSE)
###################################################
## library(yuimaGUI)
## yuimaGUI()


